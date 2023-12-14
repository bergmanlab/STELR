import os
import glob
telr_dir = config["telr_dir"]
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path


rule all:
    input:
        config["output"]
    localrule: True

"""

DOWNLOAD GENOME AND MODEL FILES

"""

def if_accession(wildcards): 
    #rudimentary check that the accession number given follows the right format eg starts with "GCA"
    #if no accession number was given, the config["input options"]["{reference} reference"] dict will not contain the "accession" key
    #in that case, this function fails at the getdict(), the exception is caught, and it returns []
    #this should only happen when the genome file was given and the accession # wasn't, so the download_genome rule will not be run either way
    #However, doing it this way prevents an error when snakemake is assessing this function while calculating the DAG.
    try:
        accession = getdict(config, ["input options", f"{wildcards.reference} reference","accession"])
        if accession[:3] == "GCA":
            return accession
    except: pass
    return []
rule download_genome:
    params:
        accession = if_accession #points to the above function, which snakemake runs when assessing the DAG
    output:
        "{reference}_reference.fasta"
    localrule: True # don't start a cluster job for this rule
    conda:
        config["conda"]["datasets"]
    shell:
        """
        mkdir ncbi_temp_{params.accession}
        cd ncbi_temp_{params.accession}
        datasets download genome accession {params.accession}
        unzip ncbi_dataset.zip
        cat ncbi_dataset/data/{params.accession}/*_genomic.fna | sed -e 's/.*chromosome />chr/g' | sed 's/, whole genome shotgun sequence//g' > raw_{output}
        samtools faidx raw_{output} chr2L chr2R chr3L chr3R chr4 chrX > {output}
        cd ..
        mv ncbi_temp_{params.accession}/{output} {output}
        rm -r ncbi_temp_{params.accession}        
        """
        #Note: I am really not sure if the sed script will work with other genome files; this may need to be amended if the eval pipeline is expanded to more genomes

rule download_model:
    #If the model is given without a file path, it downloads it from the pbsim2 creator's git repo
    params:
        model = config["simulation parameters"]["model"]
    output:
        "{model}.model"
    localrule: True
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        wget -L https://raw.githubusercontent.com/yukiteruono/pbsim2/master/data/{output} -O "{output}"
        """

"""

SIMULATE LONG READS

"""

def simulation_time(ideal_time):
    #this really still needs to be worked out
    config_time = config["slurm options"]["time limit"]
    time = str(min(int(config_time.replace(":","")), ideal_time))
    time = int("0"+time[:-4])*60 + int("0"+time[-4:-2]) + int("0"+time[-2:])/60
    print(f"config time: {config_time} ideal time: {ideal_time} returned time: {time}")
    return time
rule simulate_reads:
    #This is the rule that actually runs pbsim2 to simulate the reads
    #It does so on "sub-simulations" which are requested by the combine_simulations rule below
    #This rule is run 1-2 times in total for each simulation the pipeline performs
    input:
        reference = "{ref_type}_reference.fasta",
        model = config["simulation parameters"]["model"]["file"]
    output:
        out = "{simulation}/{cov}x_{ref_type}_simulated_reads.fq",
        log = "{simulation}/{cov}x_{ref_type}_simulated_reads.report"
    params:
        prefix = "{cov}x_{ref_type}_simulated_reads",
        sub_sim_dir = "{simulation}/subsim_{ref_type}"
    resources: 
        mem_mb = 100000,
        runtime = lambda wildcards: simulation_time(20000)
    group: "simulation"
    shell:
        """
        if [ ! -d {params.sub_sim_dir} ]; then
            mkdir {params.sub_sim_dir}
        fi
        cd {params.sub_sim_dir}
        pbsim --depth {wildcards.cov} --seed {params.seed} --prefix {params.prefix} --id-prefix '{wildcards.ref_type}' --hmm_model {input.model} ../../{input.reference} 2> {output.log}
        cd ../..
        cat {params.sub_sim_dir}/{params.prefix}*.fastq > {output.out}
 	    rm -r {params.sub_sim_dir}
        """

def sub_simulations(wildcards):
    #Determines the two subsimulations from the coverage and genotype of each simulation
    coverage = int(wildcards.cov.split("/")[-1])
    simulation = f"{wildcards.cov}x_{wildcards.proportion_genotype}"
    #Proportion dict is determined by the "proportion_genotype" wildcard eg "diploid_homozygous"
    #The "proportion" in question is the proportion of reads simulated from the mapping vs community reference
    #For a diploid heterozygous simulation eg, this is 1:1. 
    proportion = {
        "diploid_heterozygous":{"community":0.5,"mapping":0.5},
        "diploid_homozygous":{"community":1, "mapping":0},
        "tetraploid_simplex":{"community":0.25,"mapping":0.75},
        "tetraploid_duplex":{"community":0.5,"mapping":0.5},
        "tetraploid_triplex":{"community":0.75,"mapping":0.25},
        "tetraploid_quadruplex":{"community":1, "mapping":0},
    }[wildcards.proportion_genotype]
    #Multiplies the total coverage by the proportion for that reference genome, and that is the coverage for that reference
    #This is a list comprehension so that it returns only one subsimulation for homozygous simulations, rather than asking snakemake to make a "0x" simulation
    return [f"{simulation}/{proportion[reference]*coverage}x_{reference}_simulated_reads.fq" for reference in ("community","mapping") if not proportion[reference] == 0]
rule combine_simulation:
    #When making a heterozygous / mixed simulation, it simulates a proportional fraction of each genome
    #These two parts are just combined straightforwardly with cat
    input:
        sub_simulations
    output:
        "{cov}x_{proportion_genotype}/simulated_reads.fq"
    group: "simulation"
    resources: 
        mem_mb = 5000,
        runtime = lambda wildcards: simulation_time(3000)
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        cat {input} > {output}
        """


"""

RUN STELR OR TELR

"""

#run_stelr and run_telr are separate rules because the outputs have different names (.stelr vs .telr)
#and because this will make it easier to adapt this pipeline for new options which may be added to STELR in the future

rule get_mapping_annotation:
    #If running STELR, the evaluation pipeline will go ahead and generate a single repeatmasked reference annotation
    #This avoids having STELR run RepeatMasker on the reference genome individually for every simulation, which
    #is fairly resource intensive and completely redundant.
    input: "mapping_reference.fasta"
    output: "mapping_reference.fasta.te.bed"
    conda: config["telr_conda"]
    threads: config["resources"]["threads"]
    resources: 
        mem_mb = 60000,
        runtime = lambda wildcards: simulation_time(60000)#TODO: fix this
    conda: config["telr_conda"]
    shell:
        """
        python3 {config[telr]} -r {input} -t {threads} --make_annotation
        """

def stelr_command(wildcards):
    # this function returns the appropriate command for running stelr:
    # python3 + path to src; + the --resume tag if an intermediate file directory already exists for STELR
    # (this may happen if the evaluation pipeline was interrupted then resumed)
    if "telr version 1.x" in config["telr parameters"]:
        return "telr" #left in just to avoid errors during DAG calculation, but this will not actually be used.
    else:
        command = f"python3 {config['telr']}"
        stelr_dirs = glob.glob(f"{wildcards.simulation}/stelr_run_*")
        if len(stelr_dirs) == 1:
            run_number = stelr_dirs[0].split("stelr_run_")[1]
            return f"{command} --resume {run_number}"
        else: return command

checkpoint run_stelr:
    input:
        reference = "mapping_reference.fasta",
        reads = "{simulation}/simulated_reads.fq",
        library = "library.fa",
        reference_annotation = "mapping_reference.fasta.te.bed"
    output:
        "{simulation}/simulated_reads.stelr.{output}"
    params:
        polish_iterations = config["telr parameters"]["polishing iterations"],
        assembler = config["telr parameters"]["assembler"],
        polisher = config["telr parameters"]["polisher"],
        command = stelr_command,
        slurm_log = "{simulation}_stelr_slurm.log",
        slurm_err = "{simulation}_stelr_slurm.err"
    threads: config["resources"]["threads"]
    resources: 
        mem_mb = 60000,
        runtime = lambda wildcards: simulation_time(60000)
    conda: config["telr_conda"]
    shell:
        """
        {params.command} -i {input.reads} -r {input.reference} -l {input.library} -t {threads} -k -p {params.polish_iterations} --assembler {params.assembler} --polisher {params.polisher} -a {input.reference_annotation} -o {wildcards.simulation}
        """

checkpoint run_telr:
    input:
        reference = "mapping_reference.fasta",
        reads = "{simulation}/simulated_reads.fq",
        library = "library.fa"
    output:
        "{simulation}/simulated_reads.telr.{output}"
    params:
        polish_iterations = config["telr parameters"]["polishing iterations"],
        assembler = config["telr parameters"]["assembler"],
        polisher = config["telr parameters"]["polisher"],
        slurm_log = "{simulation}_stelr_slurm.log",
        slurm_err = "{simulation}_stelr_slurm.err"
    threads: config["resources"]["threads"]
    resources: 
        mem_mb = 60000,
        runtime = lambda wildcards: simulation_time(60000)
    conda: config["telr_conda"]
    shell:
        """
        telr -i {input.reads} -r {input.reference} -l {input.library} -t {threads} -k -p {params.polish_iterations} --assembler {params.assembler} --polisher {params.polisher} {params.annotation} -o {wildcards.simulation}
        """
    


"""

LIFTOVER EVALUATION

"""
#These rules are pulled from the script liftover_evaluation.py





'''
rule liftover_eval:
    input:
        telr_out = "{simulation}/simulated_reads.stelr.contig.fasta",
        region_mask = "regular_recomb.bed"
        #annotation = "liftover_nonref.bed"
    output:
        "{simulation}/liftover_eval/annotation.filter.bed"
    localrule: True
    shell:
        """
        if [ ! -d {wildcards.simulation}/liftover_eval ]; then
            mkdir {wildcards.simulation}/liftover_eval
        fi
        touch {output}
        """

rule af_eval:
    input:
        telr_out = "{simulation}/simulated_reads.stelr.contig.fasta"
        #annotation = "liftover_nonref.bed"
    output:
        "{simulation}/af_eval/telr_eval_af.json"
    localrule: True
    shell:
        """
        if [ ! -d {wildcards.simulation}/af_eval ]; then
            mkdir {wildcards.simulation}/af_eval
        fi
        touch {output}
        """

rule seq_eval:
    input:
        telr_out = "{simulation}/simulated_reads.stelr.contig.fasta"
        #annotation = "liftover_nonref.bed"
    output:
        "{simulation}/seq_eval/seq_eval.json"
    localrule: True
    shell:
        """
        if [ ! -d {wildcards.simulation}/seq_eval ]; then
            mkdir {wildcards.simulation}/seq_eval
        fi
        touch {output}
        """
'''
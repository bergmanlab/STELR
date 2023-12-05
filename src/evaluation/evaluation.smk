import os
import glob
telr_dir = config["telr_dir"]
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path


rule all:
    input:
        config["output"]
    localrule: True


def if_accession(wildcards):
    try:
        accession = getdict(config, ["input options", f"{wildcards.reference} reference","accession"])
        if accession[:3] == "GCA":
            return accession
    except: pass
    return []
rule download_genome:
    params:
        accession = if_accession
    output:
        "{reference}_reference.fasta"
    localrule: True
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

rule download_model:
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

def simulation_time(ideal_time):
    config_time = config["slurm options"]["time limit"]
    time = str(min(int(config_time.replace(":","")), ideal_time))
    time = int("0"+time[:-4])*60 + int("0"+time[-4:-2]) + int("0"+time[-2:])/60
    print(f"config time: {config_time} ideal time: {ideal_time} returned time: {time}")
    return time
rule simulate_reads:
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
    coverage = int(wildcards.cov.split("/")[-1])
    simulation = f"{wildcards.cov}x_{wildcards.proportion_genotype}"
    proportion = {
        "diploid_heterozygous":{"community":0.5,"mapping":0.5},
        "diploid_homozygous":{"community":1, "mapping":0},
        "tetraploid_simplex":{"community":0.25,"mapping":0.75},
        "tetraploid_duplex":{"community":0.5,"mapping":0.5},
        "tetraploid_triplex":{"community":0.75,"mapping":0.25},
        "tetraploid_quadruplex":{"community":1, "mapping":0},
    }[wildcards.proportion_genotype]
    return [f"{simulation}/{proportion[reference]*coverage}x_{reference}_simulated_reads.fq" for reference in ("community","mapping") if not proportion[reference] == 0]
rule combine_simulation:
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

def telr_command(wildcards):
    if "telr version 1.x" in config["telr parameters"]:
        return "telr"
    else:
        command = f"python3 {config['telr']}"
        telr_dirs = glob.glob(f"{wildcards.simulation}/telr_run_*")
        if len(telr_dirs) == 1:
            run_number = telr_dirs[0].split("telr_run_")[1]
            return f"{command} --resume {run_number}"
        else: return command

checkpoint run_telr:
    input:
        reference = "mapping_reference.fasta",
        reads = "{simulation}/simulated_reads.fq",
        library = "library.fa",
    output:
        "{simulation}/simulated_reads.telr.{output}"
    params:
        polish_iterations = config["telr parameters"]["polishing iterations"],
        assembler = config["telr parameters"]["assembler"],
        polisher = config["telr parameters"]["polisher"],
        command = telr_command,
        slurm_log = "{simulation}_stelr_slurm.log",
        slurm_err = "{simulation}_stelr_slurm.err"
    threads: config["resources"]["threads"]
    resources: 
        mem_mb = 60000,
        runtime = lambda wildcards: simulation_time(60000)
    conda: config["telr_conda"]
    shell:
        """
        {params.command} -i {input.reads} -r {input.reference} -l {input.library} -t {threads} -k -p {params.polish_iterations} --assembler {params.assembler} --polisher {params.polisher} -o {wildcards.simulation}
        """

rule liftover_annotation:
    input:
        community_reference = "community_reference.fasta",
        community_annotation = "community_annotation.bed",
        mapping_reference = "mapping_reference.fasta",
        mapping_annotation = "mapping_annotation.bed"
    output:
        "liftover_nonref.bed"
    threads: config["resources"]["threads"]
    localrule: True
    shell:
        """
        python3 {config[liftover]} --fasta1 {input.community_reference} --fasta2 {input.mapping_reference} -1 {input.community_annotation} -2 {input.mapping_reference} -o . -t {threads} -g 50 -p 50 -x "asm10" -k
        """

'''
rule liftover_eval:
    input:
        telr_out = "{simulation}/simulated_reads.telr.contig.fasta",
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
        telr_out = "{simulation}/simulated_reads.telr.contig.fasta"
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
        telr_out = "{simulation}/simulated_reads.telr.contig.fasta"
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
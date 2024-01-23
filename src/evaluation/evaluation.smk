import os
import glob
telr_dir = config["telr_dir"]
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path
import subprocess
import json
import traceback


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
        sub_sim_dir = "{simulation}/subsim_{ref_type}",
        seed = config["simulation parameters"]["seed"]
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
    input: 
        reads = "mapping_reference.fasta",
        library = "library.fa"
    output: "mapping_reference.fasta.te.bed"
    conda: config["telr_conda"]
    threads: config["resources"]["threads"]
    resources: 
        mem_mb = 60000,
        runtime = lambda wildcards: simulation_time(60000)#TODO: fix this
    conda: config["telr_conda"]
    shell:
        """
        python3 {config[telr]} -r {input.reads} -l {input.library} -t {threads} --make_annotation
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
        "{simulation}/simulated_reads.stelr.bed",
        "{simulation}/simulated_reads.stelr.te.json",
        "{simulation}/simulated_reads.stelr.loci.json",
        "{simulation}/simulated_reads.stelr.loci.fasta",
        "{simulation}/simulated_reads.stelr.te.fasta"
    params:
        polish_iterations = config["telr parameters"]["polishing iterations"],
        aligner = config["telr parameters"]["aligner"],
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
        {params.command} -i {input.reads} -r {input.reference} -l {input.library} -t {threads} -k -p {params.polish_iterations} --aligner {params.aligner} --assembler {params.assembler} --polisher {params.polisher} -a {input.reference_annotation} -o {wildcards.simulation}
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

EVALUATE STELR COORDINATE AND FAMILY PREDICTIONS

"""
#These rules are pulled from the script liftover_evaluation.py

rule region_family_filter:
    # Filter a bed file by intersection with the regular recombination region,
    # By excluding "nested" insertions ie insertions with multiple predicted families,
    # And by excluding insertions labelled as excluded families, ie "INE_1"
    input:
        input_bed = "{some_file}.bed",
        filter_region = "regular_recomb.bed"
    output:
        "{some_file}_filter.bed"
    params:
        exclude_families = ["INE_1"],
        exclude_nested_insertions = False
    run:
        # set up functions to filter file based on params
        exclude_families = set(params.exclude_families)
        filters = {}
        if params.exclude_nested_insertions:
            filters["nested_insertions"] = lambda te: "|" in te.split("\t")[3]
            if exclude_families:
                filters["exclude_families"] = lambda te: te.split("\t")[3] in exclude_families
        elif exclude_families:
            filters["exclude_families"] = lambda te: len(exclude_families.intersection(te.split("\t")[3].split("|"))) > 0
        
        # run bedtools intersect to find the intersection between the bed file and the regular recombination region
        command = ["bedtools","intersect","-a",input.input_bed,"-b",input.filter_region,"-u"]
        intersection = subprocess.run(command, capture_output=True, text=True).stdout.strip().split("\n")

        # filter the intersection using functions set up earlier based on params.
        try:
            filtered_annotation = [line for line in intersection if not any([filters[param](line) for param in filters])]
        except:
            print(intersection)
            traceback.raise_

        # write output
        with open(output[0], "w") as out:
            out.write("\n".join(filtered_annotation))
            traceback.print_exc()
            sys.exit(2)

rule count_liftover:
    # count the # of lines in the annotation liftover once instead of doing it for every simulation's evaluation
    input: "community_annotation_filter.bed"
    output: "liftover_count"
    shell: "cat {input} | wc -l > {output}"

rule prediction_quality_filter:
    # Filter STELR output further by a few quality control metrics.
    input:
        telr_json = "{simulation}/simulated_reads.{telr}.expanded.json",
        telr_bed_filtered = "{simulation}/simulated_reads.{telr}_filter.bed"
    output: 
        filtered_predictions = "{simulation}/{telr}_filtered_predictions.bed",
        prediction_count = "{simulation}/total_{telr}_predictions"
    run:
        # Create the filtered list of TE ids from filtered bed file
        with open(input.filtered_bed,"r") as input_file:
            filtered_te_ids = {"_".join(line.strip().split("\t")[:4]) for line in input_file}
        # Load the expanded output json from (s)telr
        with open(input.telr_json,"r") as input_file:
            telr_json = json.load(input_file)
        # Key each TE's output dict to its bed format string
        telr_json = {"\t".join([str(te_dict[key]) for key in ["chrom","start","end","family"]] + [".",te_dict["strand"]]):te_dict for te_dict in telr_json}

        # Define the quality checks TEs must pass to pass quality control
        quality_checks = {                                      # TEs pass each quality control step if:
            "ID": lambda value: value in filtered_te_ids,           # they are one of the filtered TEs from the previous step
            "allele_frequency": lambda value: value is not None,    # they have a predicted value for allele frequency
            "support": lambda value: value == "both_sides"          # they have both-sided support
        }

        # make a list of TEs which failed because of each quality check
        failed_tes = {value:[te for te in te_dict if not quality_checks[value](te_dict[te][value])] for value in quality_checks}
        # the filtered TEs are those which are not present in any of the quality check failure lists.
        filtered_predictions = [te for te in te_dict if not any([lambda value: te in failed_tes[value] for value in quality_checks])]

        # count and print the total number of TELR predictions that passed region, family, and quality control checks
        total_filtered_predictions = len(filtered_predictions)
        print(f"Total {wildcards.telr.upper()} Predictions (filtered): {total_filtered_predictions}")

        # write the # of total predictions to an output file
        with open(output.prediction_count,"w") as output_file:
            output_file.write(total_filtered_predictions)
        # write the list of filtered predictions in bed format to an output file
        with open(output.filtered_predictions,"w") as output_file:
            output.write("\n".join(filtered_predictions))

rule compare_with_annotation:
    # use bedtools to find the regions of overlap (and non-overlap) between TELR predictions and the annotation
    input:
        telr_filtered_predictions = "{simulation}/{telr}_filtered_predictions.bed",
        filtered_annotation = "community_annotation_filter.bed"
    output:
        overlap = "{simulation}/{telr}_annotation_overlap.bed",
        false_negatives = "{simulation}/{telr}_false_negatives.bed"
    params:
        window = 5
    shell: 
        """
        bedtools window -w {params.window} -a {input.telr_filtered_predictions} -b {input.filtered_annotation} > {output.overlap}
        bedtools window -w {params.window} -a {input.filtered_annotation} -b {input.telr_filtered_predictions} -v > {output.false_negatives}
        """

rule parse_overlap:
    # Calculate the number of true positives, false positives, and false negatives,
    # as well as the precision and recall from running STELR on this simulated dataset
    input:
        overlap = "{simulation}/{telr}_annotation_overlap.bed",
        false_negatives = "{simulation}/{telr}_false_negatives.bed",
        prediction_count = "{simulation}/total_{telr}_predictions",
        liftover_count = "liftover_count"
    output: "{simulation}/{telr}_eval_liftover.json"
    params:
        exclude_nested_insertions = False,
        relax_mode = False
    run:
        # set conditions by which two sets of families are considered a match based on params
        if params.exclude_nested_insertions: # only one family allowed per TE
            family_match = lambda overlap: overlap[3] == overlap[9]
        elif params.relax_mode: # match as long as any TE family label is present in both sets
            family_match = lambda overlap: len(set(overlap[3].split("|")).intersection(overlap[9].split("|"))) > 0
        else: # match only if all TE labels in each set are present in both
            family_match = lambda overlap: set(overlap[3].split("|")) == set(overlap[9].split("|"))
        
        # read in the telr/annotation overlap bed file
        with open(input.overlap_file,"r") as input_file:
            tes = input_file.read().split("\n")
        #read in the counts of telr filtered predictions and the annotation
        with open(input.prediction_count,"r") as input_file:
            total_filtered_predictions = int(input_file.read())
        with open(input.liftover_count,"r") as input_file:
            liftover_count = int(input_file.read())
        
        # create a list of true positives
        true_positives = [te for te in tes if family_match(te.split("\t"))]
        # calculate the # of true positives and false positives
        summary_dict = {"total predictions":total_filtered_predictions,"true positives":len(true_positives)}
        summary_dict["false positives"] = total_filtered_predictions - summary_dict["true positives"]
        # count the number of lines in the bed file containing the false negatives
        with open(input.false_negatives,"rb") as input_file:
            summary_dict["false negatives"] = sum(1 for _ in input_file)
        # calculate the precision and recall
        summary_dict["precision"] = round(summary_dict["true positives"]/total_filtered_predictions, 3)
        summary_dict["recall"] = round(summary_dict["true positives"]/liftover_count, 3)

        # print summary statistcs and write them to output file
        telr = wildcards.telr.upper()
        for stat in ["true positives","false positives","false negatives"]:
            print(f"Number of {telr} {stat}: {summary_dict[stat]}")        
        with open(output[0],"w") as output_file:
            json.dump(summary_dict,output_file)

"""

SEQUENCE EVALUATION

"""

rule index_contigs:
    input: "{simulation}/simulated_reads.{telr}.contig.fasta"
    output: "{simulation}/simulated_reads.{telr}.contig.fasta.fai"
    shell:
        """
        samtools faidx {input}
        """

# Shunhua's script sets up a list of contigs that match the criteria and then executes the sequence evaluation process on them in parallel.
# My plan is to rebuild this in the opposite direction, first setting up rules for the individual sequence evaluations and then defining which ones it runs on
# remember the filtering step, which will need to occur at the "end" rule that determines which TEs this process is run on.

rule make_stelr_te_fasta:
    # Shunhua's code calls this aligning the contig to reference; however I'm fairly certain it's only aligning the TE sequence.
    #TELR evaluation will require an extra rule ahead of this one to make a similar output to the 18_output.json one
    input: "{simulation}/{stelr_dir}/contigs/{contig}/{te}/18_output.json"
    output: "{simulation}/{stelr_dir}/contigs/{contig}/{te}/te.fasta"
    shell:
        """
        python3 {config[seq_eval]} make_stelr_te_fasta {input} {output}
        """

rule align_te_to_ref:
    input: 
        te = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/te.fasta",
        reference = "community_reference.fasta"
    output: "{simulation}/{stelr_dir}/contigs/{contig}/{te}/alignment_to_ref.paf"
    params: 
        preset = "asm10"
    shell:
        """
        minimap2 -cx {params.preset} -v 0 --secondary=no {input.reference} {input.te} > {output}
        """

rule find_corresponding_annotation:
    input:
        telr_json = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/18_output.json",
        paf_file = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/alignment_to_ref.paf",
        community_annotation = "community_annotation.bed",
        ref = "community_reference.fasta"
    output: 
        bed = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/corresponding_annotation.bed",
        fasta = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/corresponding_annotation.fasta"
    shell:
        """
        python3 {config[seq_eval]} find_corresponding_annotation {input} {output}
        """

rule get_corresponding_ref_fa:
    input: 
        telr_sequence = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/te.fasta"
        reference_sequence = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/corresponding_annotation.fasta"
    output: 
        unsorted = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/diff_unsorted.paf",
        sorted = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/diff_sorted.paf"
    shell:
        """
        if [ -s {input.reference_sequence} ]; then
            minimap2 -cx asm5 --cs {input.reference_sequence} {input.telr_sequence} > {output.unsorted}
            sort -k6,6 -k8,8n {output.unsorted} > {output.sorted}
        else
            touch {output}
        fi
        """

rule get_sequence_similarity:
    input:
        ref_paf = "",
        annotation_bed = "",
        annotation_paf = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/diff_sorted.paf",
        telr_json = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/18_output.json"
    output:
        "{simulation}/{stelr_dir}/contigs/{contig}/{te}/sequence_evaluation.json"
    params:
        flank_len = 500
    shell:
        """
        python3 {config[seq_eval]} sequence_eval_report {params.flank_len} {input} {output}
        """
        



"""

Evaluate Allele Frequency Predictions

"""


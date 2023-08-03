import os
telr_dir = config["telr_dir"]
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path


rule all:
    input:
        config["output"]


def if_get_genome(wildcards):
    genome = os.path.basename(wildcards.reference)
    genomes_in_config = [
        ("options for simulated reads","Reference Genomes","Mapping Reference"),
        ("options for simulated reads","Reference Genomes","Community Reference"),
        ("options for reads from file","Mapping Reference"),
        ("TELR parameters","options for different reference"),
    ]
    for genome_from in genomes_in_config:
        if getdict(config,genome_from+("name",)) == genome:
            return getdict(config,genome_from+("accession",))
rule download_genome:
    params:
        accession = if_get_genome
    output:
        "{reference}.fna"
    shell:
        """
        mkdir ncbi_temp_{params.accession}
        cd ncbi_temp_{params.accession}
        datasets download genome accession {params.accession}
        unzip ncbi_dataset.zip
        cd ..
        mv ncbi_temp_{params.accession}/ncbi_dataset/data/{params.accession}/*_genomic.fna {output}
        rm -r ncbi_temp_{params.accession}
        """

rule download_model:
    params:
        model = config["options for simulated reads"]["pbsim2"]["model name"]
    output:
        "{model}.model"
    shell:
        """
        wget -L https://raw.githubusercontent.com/yukiteruono/pbsim2/master/data/{wildcards.model}.model -O "{output}"
        """


rule simulate_reads:
    input:
        reference = lambda wildcards: getdict(config,("options for simulated reads","Reference Genomes",f"{wildcards.ref_type} Reference"))["file path"],
        model = config["options for simulated reads"]["pbsim2"]["model file"]
    output:
        out = "{cov}x_{ref_type}_pbsim2.fq",
        log = "{cov}x_{ref_type}_pbsim2.report"
    params:
        prefix = "'{cov}x_{ref_type}_pbsim2'"
    resources: 
        mem_mb = 100000
    shell:
        """
        mkdir {params.prefix}
        cd {params.prefix}
        pbsim --depth {wildcards.cov} --prefix {params.prefix} --id-prefix '{wildcards.ref_type}' --hmm_model '{input.model}' '{input.reference}' 2> '../{output.log}'
        cd ..
        cat {params.prefix}/*.fastq > {output.out}
        rm -r {params.prefix}
        """

def sub_simulations(wildcards):
    proportion = config["options for simulated reads"]["Simulation type params"]["proportion"]
    coverage = config["options for simulated reads"]["coverage"]
    simulator = config["options for simulated reads"]["Simulator"]
    subsims = {
        "Mapping":proportion[0],
        "Community":proportion[1]
    }
    return [f"{coverage*subsims[subsim]}x_{subsim}_{simulator}.fq" for subsim in subsims if subsims[subsim] > 0]
rule combine_simulation:
    input:
        sub_simulations
    output:
        "{cov}x_{proportion_genotype}.fastq"
    shell:
        """
        cat {input} > {output}
        """

def reads_for_telr(wildcards):
    source = config["read source"]
    if source == "reads from file":
        return config["options for reads from file"]["Long read sequencing data"]
    elif source == "simulated reads":
        cov = config["options for simulated reads"]["coverage"]
        simulation_type = config["options for simulated reads"]["Simulation type"]
        if simulation_type == "proportion":
            proportion = config["options for simulated reads"]["Simulation type params"]["proportion"]
            return f"{cov}x_{proportion[0]}-{proportion[1]}.fastq"
        else:
            genotype = config["options for simulated reads"]["Simulation type params"]["genotype"]
            return f"{cov}x_{simulation_type}_{genotype}.fastq"
def reference_for_telr(wildcards):
    if config["TELR parameters"]["Different mapping reference"]:
        return config["TELR parameters"]["options for different reference"]["file path"]
    else:
        if config["read source"] == "reads from file":
            return config["options for reads from file"]["Mapping Reference"]["file path"]
        elif config["read source"] == "simulated reads":
            return config["options for simulated reads"]["Reference Genomes"]["Mapping Reference"]["file path"]
checkpoint run_telr:
    input:
        reference = reference_for_telr,
        reads = reads_for_telr,
        library = config["TELR parameters"]["TE library"],
    output:
        "{sample_name}.telr.{output}"
    params:
        polish_iterations = config["TELR parameters"]["Polishing iterations"],
        assembler = config["TELR parameters"]["Assembler"],
        polisher = config["TELR parameters"]["Polisher"],
        command = config["TELR command"]
    threads: config["Resources"]["Threads"]
    conda: config["telr_conda"]
    shell:
        """
        {params.command} -i {input.reads} -r {input.reference} -l {input.library} -t {threads} -k -p {params.polish_iterations} --assembler {params.assembler} --polisher {params.polisher}
        """

#'''
rule telr_liftover_eval:
    input:
        telr_out = lambda wildcards: f"{config['reads_name']}.telr.contig.fasta",
        region_mask = config["Evaluation requirements"]["Regular recombination region"],
        annotation = config["Evaluation requirements"]["Mapping reference annotation"]
    output:
        "liftover_eval/annotation.filter.bed"
    shell:
        """
        if [ ! -d liftover_eval ]
        then
            mkdir liftover_eval
        fi
        python3 {config[liftover_eval]} -i . -o liftover_eval -r {input.region_mask} -b {input.annotation} --exclude_families "INE_1"
        """

def get_ploidy(wildcards): #adjust this for non tetraploid and diploid options
    if config["read source"] == "simulated reads":
        string = config["options for simulated reads"]["Simulation type"]
    else: string = os.path.basename(config["options for reads from file"]["Long read sequencing data"])
    string = string.lower()
    if "tetraploid" in string: return 4
    else: return 2
def get_genotype(wildcards):
    if config["read source"] == "simulated reads":
        string = config["options for simulated reads"]["Simulation type params"]["genotype"]
    else: string = os.path.basename(config["options for reads from file"]["Long read sequencing data"])
    string = string.lower()
    ploidy = get_ploidy(wildcards)
    if ploidy == 4:
        for genotype in ["simplex","duplex","triplex","quadruplex"]:
            if genotype in string:
                return genotype
    else:
        for genotype in ["homozygous","heterozygous"]:
            if genotype in string:
                return genotype
        return "heterozygous"
rule telr_af_eval:
    input:
        telr_out = lambda wildcards: f"{config['reads_name']}.telr.contig.fasta",
        region_mask = config["Evaluation requirements"]["Regular recombination region"]
    output:
        "af_eval/telr_eval_af.json"
    params:
        ploidy = get_ploidy,
        genotype = get_genotype
    shell:
        """
        if [ ! -d af_eval ]
        then
            mkdir af_eval
        fi
        python3 {config[af_eval]} -i . -o af_eval -r {input.region_mask} --exclude_families "INE_1" --ploidy {params.ploidy} --genotype {params.genotype}
        """

rule telr_seq_eval:
    input:
        telr_out = lambda wildcards: f"{config['reads_name']}.telr.contig.fasta",
        region_mask = config["Evaluation requirements"]["Regular recombination region"],
        annotation = config["Evaluation requirements"]["Community reference annotation"],
        community_reference = config["Evaluation requirements"]["Community reference genome"]
    output:
        "seq_eval/seq_eval.json"
    threads: config["Resources"]["Threads"]
    shell:
        """
        if [ ! -d seq_eval ]
        then
            mkdir seq_eval
        fi
        python3 {config[seq_eval]} -i . -r {input.community_reference} -o seq_eval --region1 {input.region_mask} -a {input.annotation} -t {threads} --exclude_families "INE_1"
        """


'''
LIFTOVER_EVALUATION
'''

rule rename_chroms:
    input:
        telr_out = lambda wildcards: f"{config['reads_name']}.telr.contig.fasta"

def need_to_rename(wildcards):
    with open(config["Evaluation requirements"]["Regular recombination region"], "r") as region_mask:
        list_of_chroms = {line.split("\t")[0] for line in region_mask}
    telr_output_file = checkpoints.run_telr.get(output="bed").output[0]
    with open(telr_output_file,"r") as telr_output:
        chrom_match = len({line for line in telr_output for chrom in list_of_chroms if chrom in line})
    if chrom_match == 0:
        return "renamed_out.telr.bed"#TODO: make this file
    else:
        return telr_output_file
def rename_json(wildcards):
    return need_to_rename(wildcards).replace(".telr.bed",".telr.json")
rule filter_annotation_region:
    input:
        region_mask = config["Evaluation requirements"]["Regular recombination region"],
        annotation = config["Evaluation requirements"]["Community reference annotation"]
    output:
        "liftover_eval/filter_annotation_region.bed"
    shell:
        """
        if [ ! -d liftover_eval ]
        then
            mkdir liftover_eval
        fi
        bedtools intersect -a {input.annotation} -b {input.region_mask} -u > {output}
        """

rule filter_annotation_family:
    input:
        region_mask = config["Evaluation requirements"]["Regular recombination region"],
        annotation = "liftover_eval/filter_annotation_region.bed"
    output:
        "liftover_eval/filter_annotation_family.bed"
    params:
        family_filter = "INE-1"
    shell:
        """
        python3 {config[eval_utility]} filter_family_bed {input.annotation} {params.family_filter} {output} "exclude"
        """
rule filter_telr_region:
    input:
        region_mask = config["Evaluation requirements"]["Regular recombination region"],
        telr_output = need_to_rename
    output:
        "liftover_eval/filter_telr_region.bed"
    shell:
        """
        if [ ! -d liftover_eval ]
        then
            mkdir liftover_eval
        fi
        bedtools intersect -a {input.telr_output} -b {input.region_mask} -u > {output}
        """

rule filter_telr_family:
    input:
        region_mask = config["Evaluation requirements"]["Regular recombination region"],
        telr_output = "liftover_eval/filter_telr_region.bed",
    output:
        "liftover_eval/filter_telr_family.bed"
    params:
        family_filter = "INE-1"
    shell:
        """
        python3 {config[eval_utility]} filter_family_bed {input.telr_output} {params.family_filter} {output} "exclude"
        """

rule filtered_output_with_info:#total
    input:
        "liftover_eval/filter_telr_family.bed",
        rename_json
    output:
        "liftover_eval/filtered_telr_info.bed"
    shell:
        """
        python3 {config[liftover_evaluation]} give_output_info {input} {output}
        """

rule compare_liftover:
    input:
        pred_filtered = "liftover_eval/filtered_telr_info.bed",
        annotation_filtered = "liftover_eval/filter_annotation_family.bed"
    output:
        "liftover_eval/telr_liftover_overlap.bed"
    params:
        window = 5
    shell:
        """
        bedtools window -w {params.window} -a {input.pred_filtered} -b {input.annotation_filtered} -v > {output}
        """

rule parse_overlap:#tp
    input:
        "liftover_eval/telr_liftover_overlap.bed"
    output:
        "liftover_eval/parsed_liftover_overlap.bed"
    params:
        relax_mode = "false"
    shell:
        """
        python3 {config[liftover_evaluation]} parse_overlap {input} {output} {params.relax_mode}
        """

rule liftover_only:#fn
    input:
        pred_filtered = "liftover_eval/filtered_telr_info.bed",
        annotation_filtered = "liftover_eval/filter_annotation_family.bed"
    output:
        "liftover_eval/liftover_only.bed"
    params:
        window = 5
    shell:
        """
        bedtools window -w {params.window} -a {input.annotation_filtered} -b {input.pred_filtered} -v > {output}
        """

rule parse_overlap_2:
    input:
        all_telr_predictions = "liftover_eval/filtered_telr_info.bed",
        true_positives = "liftover_eval/parsed_liftover_overlap.bed",
        false_negatives = "liftover_eval/liftover_only.bed",
        annotation_filtered = "liftover_eval/filter_annotation_family.bed"
    output:
        "liftover_eval/eval_liftover.json"
    params:
        prefix = config["reads_name"]
    shell:
        """
        python3 {config[liftover_evaluation]} parse_overlap {params.prefix} {input} {output}
        """


#echo "hello$(wc -l < reads_sort.telr.bed)"

#echo "hello$(($(wc -l < reads_sort.telr.bed) - $(wc -l < reads_sort.telr.json)))"
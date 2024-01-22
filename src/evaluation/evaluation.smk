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

rule make_te_fasta:
    # Shunhua's code calls this aligning the contig to reference; however I'm fairly certain it's only aligning the TE sequence.
    #TELR evaluation will require an extra rule ahead of this one to make a similar output to the 18_output.json one
    input: "{simulation}/{stelr_dir}/contigs/{contig}/{te}/18_output.json"
    output: "{simulation}/{stelr_dir}/contigs/{contig}/{te}/te.fasta"
    run:
        with open(input[0],"r") as contig_jsons:
            te_dict = json.load(contig_jsons)
        with open(output[0],"w") as output_file:
            output_file.write(f">{te_dict['contig_id']}\n{te_dict['te_sequence']}")

#Test script: cd to te dir, python3, and then execute the following:
'''
import json
with open("18_output.json","r") as contig_jsons:
    te_dict = json.load(contig_jsons)

with open("te.fasta","w") as output_file:
    output_file.write(f">{te_dict['contig_id']}\n{te_dict['te_sequence']}")
'''

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

#Test script:
'''
minimap2 -cx asm10 -v 0 --secondary=no ../../../../../../community_reference.fasta te.fasta > alignment_to_ref.paf
'''

class paf_file:
    def __init__(file_path):
        def format_paf_line(line):
            if line.strip():
                line = line.split("\t")
                for n in [1,7,8,9,10,11]:
                    line[n] = int(line[n])
                return line
            else:
                return None
        with open(input.paf_file,"r") as paf_file:
            paf_data = [format_paf_line(line) for line in paf_file]
            self.data = [line for line in paf_data if line]
        
        max_len = max([line[10] for line in self.data])
        self.longest = [line for line in self.data if line[10] == max_len][0]
    
    def count():
        return len(self.data)
    
    def get(key):
        if not "characteristics" in self.__dict__:
            header = ["","query_len","","","","chrom","","start","end","matches","align_len","map_qual"]
            self.characteristics = {key:self.longest[header.index(key)] for key in header if key}
            self.characteristics["map_prop"] = self.characteristics["matches"]/self.characteristics["query_len"],
            self.characteristics["blast_id"] = self.characteristics["matches"]/self.characteristics["align_len"]
        return self.characteristics[key]

class bed_file:
    def __init__(bed_file):
        if type(bed_file) is list:
            self.data = bed_file
        else:
            def format_bed_line(line):
                if line.strip():
                    line = line.split("\t")
                    for n in [1,2]:
                        line[n] = int(line[n])
                    return line
                else: return [0,0,0,0,0,0]
            with open(file_path,"r") as input_file:
                self.data = [format_bed_line(line) for line in input_file if line.strip()]
    
    def get(key):
        if not "characteristics" in self.__dict__:
            header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
            self.characteristics = {key:self.data[0][header.index(key)] for key in header if key}
        return self.characteristics[key]
    
    def write_out(file_path):
        if type(self.data[0]) is list:
            lines = ["\t".join([str(i) for i in line]) for line in self.data]
        else: lines = ["\t".join([str(i) for i in self.data])]
        with open(file_path,"w") as output_file:
            output_file.write("\n".join(lines))

rule find_corresponding_annotation:
    input:
        telr_json = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/18_output.json",
        paf_file = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/alignment_to_ref.paf",
        community_annotation = "community_annotation.bed",
        ref = "community_reference.fasta"
    output: 
        bed = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/corresponding_annotation.bed",
        fasta = "{simulation}/{stelr_dir}/contigs/{contig}/{te}/corresponding_annotation.fasta"
    run:
        with open(input.telr_json,"r") as telr_json:
            telr_data = json.load(telr_json)

        #paf_data = get_paf_data(input.paf_file)
        ref_paf = paf_file(input.paf_file)
        ref_align_start = ref_paf.get("start")
        ref_align_end = ref_paf.get("end")

        annotation_data = bed_file(input.community_annotation).data
        
        aligned_chrom = ref_paf.longest[5]
        candidate_annotations = [line for line in annotation_data if line[0] == aligned_chrom and line[3] == telr_data["family"]]
        location_filters = [
            lambda start, end: start >= ref_align_start and end <= ref_align_end,
            lambda start, end: start <= ref_align_end and end >= ref_align_end,
            lambda start, end: start <= ref_align_start and end >= ref_align_start
        ]
        overlapping_annotations = [line for line in candidate_annotations if any([f(line[1],line[2]) for f in location_filters])]

        if len(overlapping_annotations) == 1:
            corresponding_annotation = overlapping_annotations[0]
            bed_file(corresponding_annotation).write_out(output.bed)
            coord = f"{corresponding_annotation[0]}:{corresponding_annotation[1]+1}-{corresponding_annotation[2]}"
            with open(output.fasta,"w") as output_file:
                subprocess.run(["samtools","faidx",input.ref,coord], stdout=output_file)
        else:
            if len(overlapping_annotations) > 1:
                print(f"{wildcards.contig}: more than one TE in the aligned region")
            else: print(f"{wildcards.contig}: no ref TEs in the aligned region"):
            with open(output[0],"w") as output_file: output_file.write("")
            with open(output[1],"w") as output_file: output_file.write("")

#Test script (python3):
#working on this rn
'''
import json
with open("18_output.json","r") as telr_json:
    telr_data = json.load(telr_json)

def format_line(line):
    if line.strip():
        line = line.split("\t")
        for n in [7,8,10]:
            line[n] = int(line[n])
        return line
    else:
        return [0,0,0,0,0,0,0,0,0,0,0]

with open("alignment_to_ref.paf","r") as paf_file:
    paf_data = [format_line(line) for line in paf_file]

max_alignment_length = max([line[10] for line in paf_data])
longest_alignment = [line for line in paf_data if line[10] == max_alignment_length][0]
ref_align_start = longest_alignment[7]
ref_align_end = longest_alignment[8]

def format_bed_line(line):
    if line.strip():
        line = line.split("\t")
        for n in [1,2]:
            line[n] = int(line[n])
        return line
    else: return [0,0,0,0,0,0]

with open("../../../../../../community_annotation.bed","r") as community_annotation:
    annotation_data = [format_bed_line(line) for line in community_annotation if line.strip()]

aligned_chrom = longest_alignment[5]
candidate_annotations = [line for line in annotation_data if line[0] == aligned_chrom and line[3] == telr_data["family"]]
location_filters = [
    lambda start, end: start >= ref_align_start and end <= ref_align_end,
    lambda start, end: start <= ref_align_end and end >= ref_align_end,
    lambda start, end: start <= ref_align_start and end >= ref_align_start
]
overlapping_annotations = [line for line in candidate_annotations if any([f(line[1],line[2]) for f in location_filters])]

if len(overlapping_annotations) == 1:
    import subprocess
    corresponding_annotation = overlapping_annotations[0]
    with open("corresponding_annotation.bed","w") as output_file:
        output_file.write("\t".join([str(i) for i in corresponding_annotation]))
    coord = f"{corresponding_annotation[0]}:{corresponding_annotation[1]+1}-{corresponding_annotation[2]}"
    with open("corresponding_annotation.fasta","w") as output_file:
        subprocess.run(["samtools","faidx","../../../../../../community_reference.fasta",coord], stdout=output_file)
else:
    if len(overlapping_annotations) > 1:
        print(f"{wildcards.contig}: more than one TE in the aligned region")
    else: print(f"{wildcards.contig}: no ref TEs in the aligned region"):
    with open("corresponding_annotation.bed","w") as output_file: 
        output_file.write("")
    with open("corresponding_annotation.fasta","w") as output_file: 
        output_file.write("")

'''

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

#Test script
'''
if [ -s corresponding_annotation.fasta ]; then
    minimap2 -cx asm5 --cs corresponding_annotation.fasta te.fasta > diff_unsorted.paf
    sort -k6,6 -k8,8n diff_unsorted.paf > diff_sorted.paf
else
    touch {output}
fi
'''

# at this point I run into the get_paftools_var_summary function in Shunhua's script and I think I really 
# need to play with the intermediate output of paftools.js to understand fully what it is counting.

# currently mid implementation of compare_contig_ref_te_seqs which is nearly the last step of the single sequence evaluation process

def get_paftools_summary(paf_file):
    paftools_fields = {
        "insertions":{
            "1bp":0,
            "2bp":0,
            "[3,50)":0,
            "[50,1000)":0,
            "total":0
        },
        "deletions":{
            "1bp":0,
            "2bp":0,
            "[3,50)":0,
            "[50,1000)":0,
            "total":0
        }
        "reference bases covered":0,
        "substitutions":0
    }
    
    paftools_summary = subprocess.run(["paftools.js","call","-l","100","-L","100",paf_file], capture_output=True, text=True).stderr

    for line in paftools_summary.split("\n"):
        try:
            value = int(line.split()[0])
            field = [field for field in paftools_fields if field in line][0]
            if type(paftools_fields[field]) is int: paftools_fields[field] = value
            else:
                subfield = [subfield for subfield in paftools_fields[field] if subfield in line][0]
                paftools_fields[field][subfield] = value
                paftools_fields[field]["total"] += value
        except: pass
    
    return paftools_fields


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
    run:
        ref_paf = paf_file(input.ref_paf)
        annotation_bed = bed_file(input.annotation_bed)
        annotation_paf = paf_file(input.annotation_paf)
        paftools_summary = get_paftools_summary(input.paf_file)
        with open(input.telr_json,"r") as input_file:
            json_data = json.load(input_file)
        


        eval_report = {
            "contig_te_plus_flank_start": max(0,json_data["contig_te_start"] - params.flank_len),
            "contig_te_plus_flank_end": min(json_data["contig_length"], json_data["contig_te_end"] + params.flank_len),
            "contig_te_plus_flank_size": None,#calculated below
            "num_contig_ref_hits": ref_paf.count(),
            "ref_aligned_chrom": ref_paf.get("chrom"),
            "ref_aligned_start": ref_paf.get("start"),
            "ref_aligned_end": ref_paf.get("end"),
            "contig_max_base_mapped_prop": ref_paf.get("map_prop"),
            "contig_mapp_qual": ref_paf.get("map_qual"),
            "contig_num_residue_matches": ref_paf.get("matches"),
            "contig_alignment_block_length": ref_paf.get("align_len"),
            "contig_blast_identity": ref_paf.get("blast_id"),
            "ref_te_family": annotation_bed.get("name"),
            "ref_te_length": annotation_bed.get("end") - annotation_bed.get("start"),
            "num_contig_ref_te_hits": annotation_paf.count(),
            "ref_te_aligned_chrom": annotation_paf.get("chrom"),
            "ref_te_aligned_start": annotation_paf.get("start"),
            "ref_te_aligned_end": annotation_paf.get("end"),
            "contig_te_mapp_qual": annotation_paf.get("map_qual"),
            "contig_te_max_base_mapped_prop": annotation_paf.get("map_prop"),
            "contig_te_num_residue_matches": annotation_paf.get("matches"),
            "contig_te_alignment_block_length": annotation_paf.get("align_len"),
            "contig_te_blast_identity": annotation_paf.get("blast_id"),
            "ref_te_num_bases_covered": paftools_summary["reference bases covered"],
            "ref_te_num_snvs": paftools_summary["substitutions"],
            "ref_te_num_1bp_del": paftools_summary["deletions"]["1bp"],
            "ref_te_num_1bp_ins": paftools_summary["insertions"]["1bp"],
            "ref_te_num_2bp_del": paftools_summary["deletions"]["2bp"],
            "ref_te_num_2bp_ins": paftools_summary["insertions"]["2bp"],
            "ref_te_num_50bp_del": paftools_summary["deletions"]["[3,50)"],
            "ref_te_num_50bp_ins": paftools_summary["insertions"]["[3,50)"],
            "ref_te_num_1kb_del": paftools_summary["deletions"]["[50,1000)"],
            "ref_te_num_1kb_ins": paftools_summary["insertions"]["[50,1000)"],
            "ref_te_num_ins": paftools_summary["insertions"]["total"],
            "ref_te_num_del": paftools_summary["deletions"]["total"],
            "ref_te_num_indel": paftools_summary["deletions"]["total"] + paftools_summary["insertions"]["total"],
        }
        eval_report["contig_te_plus_flank_size"] = eval_report["contig_te_plus_flank_end"] - eval_report["contig_te_plus_flank_start"]
        with open(output[0],"w") as output_file:
            json.dump(eval_report, output_file)
        



"""

Evaluate Allele Frequency Predictions

"""


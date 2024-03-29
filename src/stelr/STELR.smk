import os
import subprocess
import json
import glob
from pathlib import Path
from STELR_utility import get_contig_length, progress_bar, abs_path

rule all:
    input: 
        config["output"]

def input_reads_if_in_bam_format(wildcards):
    #   this rule returns the bam format input if the input is given in bam format, or an empty list if not.
    #returning an empty list essentially gives snakemake permission to accept the input file in fasta 
    #format; otherwise it tries to check the next file back in the workflow, and causes a cyclic dependency
    #error if the input and output for the rule bam_input are the same.
    if ".bam" in config["reads"]: return config["reads"]
    else: return []
rule bam_input: #if input is given in bam format, convert it to fasta format.
    #todo -- check if this bam format is ACTUALLY aligned.
    input:
        input_reads_if_in_bam_format
    output:
        config["fasta_reads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_alignment]} bam2fasta '{input}' '{output}'"

"""
1st stage: identify TE insertion candidate loci
"""

'''1st stage
"Read alignment (NGMLR)"
^(or minimap2)
'''
#only run if reads are supplied in fasta format
def find_alignment(wildcards):
    bam_input = f"input/reads.bam"
    if(os.path.isfile(bam_input)): return bam_input
    else: 
        return f"{config['aligner']}_alignment.sam"
rule sort_index_bam:
    input:
        find_alignment#gives name of input file (this depends on whether the user supplied input was in bam format or it was aligned later)
    output:
        "reads_sort.bam"
    threads: config["thread"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        samtools index -@ {threads} {output}
        """
    
rule align_with_ngmlr:#TODO: add timers to these alignments
    priority: 10
    input:
        reads = config["fasta_reads"],
        reference = config["reference"]
    output:
        "ngmlr_alignment.sam"
    params:
        presets = config["presets"],
        label = lambda wildcards: {"ont":"ont","pacbio":"pb"}[config["presets"]]
    threads: config["thread"]
    conda:
        config["conda"]["ngmlr"]
    shell:
        """
        ngmlr -r {input.reference} -q {input.reads} -x {params.presets} -t {threads} --rg-id reads --rg-sm reads --rg-lb {params.label} --no-progress | python3 {config[fix_ngmlr]} > {output}
        """

rule align_with_minimap2:
    priority: 10
    input:
        reads = config["fasta_reads"],
        reference = config["reference"]
    output:
        "minimap2_alignment.sam"
    params:
        presets = lambda wildcards: {"ont":"map-ont","pacbio":"map-pb"}[config["presets"]]
    threads: config["thread"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        minimap2 --cs --MD -Y -L -ax {params.presets} {input.reference} {input.reads} > {output}
        """


'''1st stage
SV calling (Sniffles)
'''

rule sv_detection_sniffles1:
    input:
        "reads_sort.bam"
    output:
        "sv-reads_Sniffles1.vcf"
    threads: config["thread"]
    conda:
        config["conda"]["sniffles1"]
    shell:
        "sniffles -n -1 --threads {threads} -m {input} -v {output}"

rule sv_detection_sniffles2:
    input:
        "reads_sort.bam"
    output:
        "sv-reads_Sniffles2.vcf"
    threads: config["thread"]
    conda:
        config["conda"]["sniffles2"]
    shell:
        """
        sniffles --output-rnames -t 10 -i {input} -v {output}
        """
def sv_detector(wildcards):
    if "sv_detector" in config:
        return config["sv_detector"]
    else:
        return "Sniffles1"
def parse_vcf_input(wildcards):
    return f"sv-reads_{sv_detector(wildcards)}.vcf"
rule parse_vcf:
    input:
        parse_vcf_input
    output:
        "reads.vcf_parsed.tsv.tmp"
    params:
        lambda wildcards: {
            "Sniffles1":'%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n',
            "Sniffles2":'%CHROM\\t%POS\\t%END\\t%SVLEN\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n'
        }[sv_detector(wildcards)]
    conda:
        config["conda"]["stable_environment"]
    shell:
        'bcftools query -i \'SVTYPE="INS" & ALT!="<INS>"\' -f "{params}" "{input}" > "{output}"'

#temporary measure to handle Sniffles1 output formatting issue, which causes an error when calling bcftools on it
rule bcftools_bypass:
    input:
        parse_vcf_input
    output:
        ""#"reads.vcf_parsed.tsv.tmp"
    shell:
        """
        python3 {config[STELR_sv]} bcftools {input} {output}
        """

rule swap_vcf_coordinate:
    input:
        "reads.vcf_parsed.tsv.tmp"
    output:
        "reads.vcf_parsed.tsv.swap"
    params:
        lambda wildcards: sv_detector(wildcards)
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_sv]} swap_coordinate '{input}' '{output}' {params}"

rule rm_vcf_redundancy:
    input:
        "reads.vcf_parsed.tsv.swap"
    output:
        "reads.vcf_parsed.tsv"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_sv]} rm_vcf_redundancy '{input}' '{output}'"

rule write_ins_seqs:
    input:
        "reads.vcf_parsed.tsv"
    output:
        "reads.vcf_ins.fasta"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_sv]} write_ins_seqs '{input}' '{output}'"

'''1st stage
Filter for TE insertion candidate (RepeatMasker)
'''

rule sv_repeatmask:
    input:
        ins_seqs = "reads.vcf_ins.fasta",
        library = config["library"]
    output:
        "vcf_ins_repeatmask/{ins_seqs}.out.gff"
    params:
        repeatmasker_dir = "vcf_ins_repeatmask",
    threads: 20
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        RepeatMasker -dir '{params.repeatmasker_dir}' -gff -s -nolow -no_is -xsmall -e ncbi -lib '{input.library}' -pa '{threads}' '{input.ins_seqs}'
        """

rule sv_RM_sort:
    input:
        "vcf_ins_repeatmask/{ins_seqs}.out.gff"
    output:
        "vcf_ins_repeatmask/{ins_seqs}.out.sort.gff"
        #bedtools
    conda:
        config["conda"]["stable_environment"]
    shell:
        "bedtools sort -i '{input}' > '{output}'"

rule sv_RM_merge:
    input:
        "vcf_ins_repeatmask/{ins_seqs}.out.sort.gff"
    output:
        "vcf_ins_repeatmask/{ins_seqs}.out.merge.bed"
        #bedtools
    conda:
        config["conda"]["stable_environment"]
    shell:
        "bedtools merge -i '{input}' > '{output}'"

rule sv_TE_extract:
    input:
        parsed_vcf = "reads.vcf_parsed.tsv",
        ins_seqs = "reads.vcf_ins.fasta",
        ins_rm_merge = "vcf_ins_repeatmask/reads.vcf_ins.fasta.out.merge.bed"
    output:
        ins_filtered = "reads.vcf.filtered.tmp.tsv",
        loci_eval = "reads.loci_eval.tsv"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_sv]} te_extract '{input.parsed_vcf}' '{input.ins_seqs}' '{input.ins_rm_merge}' '{output.ins_filtered}' '{output.loci_eval}'"

rule seq_merge:
    input:
        "reads.vcf.filtered.tmp.tsv"
    output: 
        "reads.vcf.merged.tmp.tsv"
    params:
        window = 20
        #bedtools
    conda:
        config["conda"]["stable_environment"]
    shell:
        'bedtools merge -o collapse -c 2,3,4,5,6,7,8,9,10,11,12,13,14 -delim ";" -d "{params.window}" -i "{input}" > "{output}"'

rule merge_parsed_vcf:##### can we thread back to here? (probably not easily)
    input:
        "reads.vcf.merged.tmp.tsv"
    output:
        "reads.vcf_filtered.tsv"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_sv]} merge_vcf '{input}' '{output}'"

"""
2nd stage: assembly and polish local TE contig
"""
'''
rule initialize_contig_dir:
    input:
        "reads.vcf_filtered.tsv"
    output:
        "contigs/{contig}/00_vcf_parsed.tsv"
    shell:
        "python3 {config[STELR_assembly]} make_contig_dir '{input}' '{wildcards.contig}' '{output}'"
'''
checkpoint initialize_contig_dirs:
    input:
        "reads.vcf_filtered.tsv",
        "config.json"
    output:
        "contigs/.contig_dirs_made"
    threads: config["thread"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_assembly]} make_contig_dirs {input} {threads}
        touch {output}
        """

def contig_threads(wildcards):
    checkpoints.initialize_contig_dirs.get()
    contig_config = f"contigs/{wildcards.contig}/config.json"
    with open(contig_config,"r") as input:
        contig_config = json.load(input)
    return contig_config["threads"]
rule contig_smk:
    input:
        json = "contigs/{contig}/config.json",
        bam = "reads_sort.bam",
        ref_repeatmask = lambda wildcards: f"ref_repeatmask/{os.path.basename(config['reference'])}.te.bed"
    output:
        "contigs/{contig}/.complete"
    threads: contig_threads
    run:
        contig_dir = f"contigs/{wildcards.contig}"
        progress=progress_bar()
        print(f"Now processing contig {progress}.")
        command = [
            "snakemake","-s",config["STELR_contig"],
            "--use-conda","--conda-prefix",config["envs"],
            "--configfile",abs_path(input.json),
            "--cores",str(threads)
        ]
        try:
            with open(f"{contig_dir}/log", "a") as log:
                if len(glob.glob(f"{contig_dir}/*")) > 2:
                    subprocess.run(command + ["--unlock"], cwd=contig_dir, stdout=log, stderr=log)
                    subprocess.run(command + ["--rerun-incomplete"] + ["--rerun-triggers","mtime"], cwd=contig_dir, stdout=log, stderr=log)
                else:
                    subprocess.run(command, cwd=contig_dir, stdout=log, stderr=log)
        except: pass
        Path(output[0]).touch()
        print(f"Finished processing contig {progress}.")

# repeatmask reference genome using custom TE library
#   Not sure which step in the workflow this actually belongs to
rule ref_repeatmask:
    priority: 1
    input:
        ref = config["reference"],
        lib = config["library"]
    output:
        "ref_repeatmask/{reference}.masked",
        "ref_repeatmask/{reference}.out.gff",
        "ref_repeatmask/{reference}.out"
    params:
        ref_rm_dir = "ref_repeatmask"
    threads: 20
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ ! -d '{params.ref_rm_dir}' ]; then mkdir '{params.ref_rm_dir}'
        fi
        RepeatMasker -dir '{params.ref_rm_dir}' -gff -s -nolow -no_is -e ncbi -lib '{input.lib}' -pa '{threads}' '{input.ref}'
        touch {output}
        """

rule ref_rm_process:
    input:
        gff = "ref_repeatmask/{reference}.out.gff",
        out = "ref_repeatmask/{reference}.out"
    output:
        "ref_repeatmask/{reference}.out.gff3"
    shell:
        "python3 {config[STELR_te]} parse_rm_out '{input.gff}' '{input.out}' '{output}'"
        #left off here, config[STELR_te]} repeatmask()

rule ref_te_bed:
    input:
        "ref_repeatmask/{reference}.out.gff3"
    output:
        "ref_repeatmask/{reference}.te.bed.unsorted"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_te]} gff3tobed '{input}' '{output}'
        touch '{output}'
        """

rule sort_ref_rm:
    input:
        "ref_repeatmask/{reference}.te.bed.unsorted"
    output:
        "ref_repeatmask/{reference}.te.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input}' ]; then
            bedtools sort -i '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """


'''
Write Output
'''

def all_contigs_output(wildcards):
    checkpoints.initialize_contig_dirs.get(**wildcards)
    contig_list = [x.split("/")[1] for x in glob.glob("contigs/*/config.json")]
    return [f"contigs/{contig}/.complete" for contig in contig_list]
rule compile_contig_info:
    input: all_contigs_output
    output: "all_tes.json"
    threads: config["thread"]
    shell:
        """
        python3 {config[STELR_output]} compile_te_data {output} {threads}
        """
rule contig_fasta_output:
    input: "all_tes.json"
    output: "reads.stelr.loci.fasta"
    threads: config["thread"]
    params: 
        contig_pattern = "03_contig1.fa"
    shell:
        """
        python3 {config[STELR_output]} write_contig_fasta_output {input} {output} {threads} {params.contig_pattern}
        """
rule te_fasta_output:
    input: "all_tes.json"
    output: "reads.stelr.te.fasta"
    shell:
        """
        python3 {config[STELR_output]} write_te_fasta_output {input} {output}
        """
rule bed_output:
    input: "all_tes.json"
    output: "reads.stelr.bed"
    shell:
        """
        python3 {config[STELR_output]} write_bed_output {input} {output}
        """
rule te_json_output:
    input: "all_tes.json"
    output: "reads.stelr.te.json"
    shell:
        """
        python3 {config[STELR_output]} write_te_json {input} {output}
        """
rule contig_json_output:
    input: 
        parsed_SVs = "reads.vcf_parsed.tsv",
        checkpoint = "all_tes.json"
    output: "reads.stelr.loci.json"
    threads: config["thread"]
    shell:
        """
        python3 {config[STELR_output]} write_contig_json {input.parsed_SVs} {output} {threads}
        """
rule vcf_output:
    input: 
        stelr_json = "reads.stelr.te.json",
        reference = config["reference"],
        reference_index = lambda wildcards: f"{config['reference']}.fai"
    output: "reads.stelr.vcf"
    shell:
        """
        python3 {config[STELR_output]} write_vcf_output {input} {output}
        """

'''
Util Rules
'''

rule minimap2bed:
    input:
        "{minimap_output}.paf"
    output:
        "{minimap_output}.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_utility]} minimap2bed '{input}' '{output}'"

rule build_index:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        samtools faidx '{input}' || true
        touch {output}
        """

rule installation:
    output:
        ".installed_{env}"
    conda:
        lambda wildcards: config["conda"][wildcards.env]
    shell:
        "touch {output}"
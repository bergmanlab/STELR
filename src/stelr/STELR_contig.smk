import os
import json
import traceback


def all_tes(wildcards): #expects annotation file to be in contigs/{contig}/tes/
    tes = []
    annotation_file = checkpoints.annotate_contig.get(**wildcards).output[0]
    with open(annotation_file, "r") as input:
        for line in input:
            entry = line.replace("\n","").split("\t")
            if len(entry) == 6:
                tes.append(f"te_{entry[1]}_{entry[2]}")
    return tes
def get_all_outputs(wildcards):
    tes = all_tes(wildcards)
    return [f"tes/{te}/18_output.json" for te in tes]
rule all:
    input:
        get_all_outputs
    shell:
        "touch .complete"

'''2nd stage
Local contig assembly and polishing (wtdbg2/flye + minimap2)
'''


rule get_read_ids: # get a list of all the read IDs from the parsed vcf file
    input:
        "00_vcf_parsed.tsv"
    output:
        "00_reads.id"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_assembly]} write_read_IDs '{input}' '{config[contig_name]}' '{output}'"

rule unique_IDlist: # get a list of unique IDs from the readlist
    input:
        "{read_ids}.id"
    output:
        "{read_ids}.id.unique"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "cat '{input}' | sort | uniq > '{output}'"

def fasta_reads(wildcards):
    fasta_reads = config["fasta_reads"]
    if fasta_reads[0] == "/":
        return fasta_reads
    else: return f"../../{fasta_reads}"
rule filter_readlist: # use seqtk to get the fasta reads from the input reads file
    input:
        reads = fasta_reads,
        unique = "{read_ids}.id.unique"
    output:
        "{read_ids}.fa"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "seqtk subseq '{input.reads}' '{input.unique}' | seqtk seq -a > '{output}'"

rule run_assembly:
    input:
        "00_reads.fa"
    output:
        "01_initial_assembly.fa"
    threads: 1
    conda:
        config["conda"][config["assembler"]]
    shell:
        """
        python3 {config[STELR_assembly]} run_{config[assembler]}_assembly '{input}' '{config[contig_name]}' '{threads}' '{config[presets]}' '{output}'
        """

rule run_polishing:
    input:
        reads = "00_reads.fa",
        initial_assembly = "01_initial_assembly.fa"
    output:
        "02_polished_assembly.fa"
    threads: 1
    conda:
        config["conda"][config["polisher"]]
    shell:
        """
        python3 {config[STELR_assembly]} run_{config[polisher]}_polishing '{input.initial_assembly}' '{output}' '{input.reads}' '{config[contig_name]}' '{threads}' '{config[polish_iterations]}' '{config[presets]}'
        """

rule get_parsed_contigs:
    input:
        "02_polished_assembly.fa"
    output:
        "03_contig1.fa"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_assembly]} parse_assembled_contig '{input}' '{config[contig_name]}' '{output}'
        """

"""
3rd stage: annotate TE and predict location in reference genome
"""

'''3rd stage
Contig TE annotation (minimap2 + RepeatMasker)
'''

rule get_vcf_seq:
    input:
        vcf_parsed = "00_vcf_parsed.tsv"
    output:
        "04_vcf_seq.fa"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_te]} get_vcf_seq '{config[contig_name]}' '{input.vcf_parsed}' '{output}'
        """

rule map_contig:
    input:
        subject = "03_contig1.fa",
        query = "04_vcf_seq.fa"
    output:
        "05_vcf_mm2.paf"
    params:
        presets = lambda wildcards: {"pacbio":"map-pb","ont":"map-ont"}[config["presets"]]
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        minimap2 -cx '{params.presets}' --secondary=no -v 0 -t '{threads}' '{input.subject}' '{input.query}' > '{output}'
        """

rule te_contig_map:
    input:
        minimap_initial = "05_vcf_mm2.paf",
        subject = "03_contig1.fa",
        library = config["library"]
    output:
        "06_te_mm2.paf"
    params:
        presets = lambda wildcards: {"pacbio":"map-pb","ont":"map-ont"}[config["presets"]]
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        minimap2 -cx '{params.presets}' '{input.subject}' '{input.library}' -v 0 -t '{threads}' > '{output}'
        """

rule minimap2bed:
    input:
        "{minimap_output}_mm2.paf"
    output:
        "{minimap_output}_mm2.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_utility]} minimap2bed '{input}' '{output}'"

rule vcf_alignment_filter_intersect:
    input:
        vcf_seq_mm2 = "05_vcf_mm2.bed",
        te_mm2 = "06_te_mm2.bed"
    output:
        "07_te2contig_filter.tsv"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        bedtools intersect -a '{input.te_mm2}' -b '{input.vcf_seq_mm2}' -wao > '{output}'
        """

rule vcf_alignment_filter:
    input:
        "07_te2contig_filter.tsv"
    output:
        "08_te2contig_filtered.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_te]} vcf_alignment_filter '{input}' '{output}'
        """

rule te_annotation_sort:
    input:
        "08_te2contig_filtered.bed"
    output:
        "08_te2contig_sorted.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        bedtools sort -i {input} > {output}
        """

rule te_annotation_merge:
    input:
        "08_te2contig_sorted.bed"
    output:
        "09_te2contig_merged.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim "|" -i '{input}' > '{output}'
        """

checkpoint annotate_contig:
    input:
        "09_te2contig_merged.bed"
    output:
        "tes/annotation.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_te]} annotate_contig '{input}' '{output}'
        """

## use RM to annotate config

rule rm_te_fasta:
    input:
        bed_file = "tes/annotation.bed",
        sequence = "03_contig1.fa"
    output:
        "rm_01_annotated_tes.fasta"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        bedtools getfasta -fi '{input.sequence}' -bed '{input.bed_file}' > '{output}'
        """

rule rm_annotate:
    input:
        te_fasta = "rm_01_annotated_tes.fasta",
        te_library = config["library"]
    output:
        "rm_01_annotated_tes.fasta.out.gff"
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        RepeatMasker -gff -s -nolow -no_is -xsmall -e ncbi -lib '{input.te_library}' -pa '{threads}' '{input.te_fasta}'
        """

rule rm_annotation_sort:
    input:
        "rm_01_annotated_tes.fasta.out.gff"
    output:
        "rm_02_annotated_tes.out.sort.gff"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        bedtools sort -i '{input}' > '{output}'
        """

rule rm_annotation_parse_merge:
    input:
        "rm_02_annotated_tes.out.sort.gff"
    output:
        "rm_03_annotated_tes_parsed.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_te]} rm_parse_merge '{input}' '{output}'
        """

rule rm_annotation_bedtools_merge:
    input:
        "rm_03_annotated_tes_parsed.bed"
    output:
        "rm_04_annotated_tes_merged.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        bedtools merge -c 4,6 -o distinct -delim "|" -i '{input}' > '{output}'
        """

rule rm_reannotate:
    input:
        repeat_masker_out = "rm_04_annotated_tes_merged.bed",
        original_bed = "tes/annotation.bed"
    output:
        "rm_05_annotation.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_te]} rm_reannotate '{input.repeat_masker_out}' '{input.original_bed}' '{output}'
        """

'''3rd stage
Identify TE insertion breakpoint (minimap2)
'''

def te_annotation(wildcards):
    return {True:f"tes/annotation.bed",False:f"rm_05_annotation.bed"}[config["minimap2_family"]]
rule make_te_dirs:
    input:
        te_annotation
    output:
        "tes/{te}/00_annotation.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_te]} make_te_dir '{input}' '{output}'"

rule make_te_json:
    input:
        "tes/{te}/00_annotation.bed"
    output:
        "tes/{te}/00_annotation.json"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_liftover]} make_json '{input}' '{output}'"

rule build_index:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        samtools faidx '{input}' || true
        touch {output}
        """

rule get_genome_size:
    input:
        "{genome}.fai"
    output:
        "{genome}.size"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_liftover]} get_genome_size '{input}' '{output}'"

rule flank_bed:
    input:
        fasta = "03_contig1.fa",
        contig_size = "03_contig1.fa.size",
        te_dict = "tes/{te}/00_annotation.json"
    output:
        "tes/{te}/12_{flank}_flank.bed"
    params:
        flank_len = config["flank_len"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_liftover]} flank_bed '{input.fasta}' '{input.contig_size}' '{input.te_dict}' '{params.flank_len}' '{output}'
        touch '{output}'
        """

rule flank_fasta:
    input:
        fasta = "03_contig1.fa",
        bed = "tes/{te}/12_{flank}_flank.bed"
    output:
        "tes/{te}/12_{flank}_flank.fa"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input.bed}' ]; then
            bedtools getfasta -fi '{input.fasta}' -bed '{input.bed}' -fo '{output}'
        else
            touch '{output}'
        fi
        """

rule align_flank:
    input:
        flank_fa = "tes/{te}/12_{flank}_flank.fa",
        ref_fa = config["reference"]
    output:
        "tes/{te}/13_{flank}_flank.paf"
    params:
        preset = "asm10",
        num_secondary = 10
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input.flank_fa}' ]; then
            minimap2 -t {threads} -cx '{params.preset}' -v 0 -N '{params.num_secondary}' '{input.ref_fa}' '{input.flank_fa}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule get_flank_alignment_info:
    input:
        "tes/{te}/13_{flank}_flank.paf"
    output:
        "tes/{te}/14_{flank}_flank.info"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_liftover]} get_paf_info '{input}' '{output}'
        else
            touch '{output}'
        fi
        """

rule flank_paf_to_bed:
    input:
        "tes/{te}/13_{flank}_flank.paf"
    output:
        "tes/{te}/14_{flank}_flank.bed_unsorted"
    params:
        different_contig_name = config["different_contig_name"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_liftover]} paf_to_bed '{input}' '{output}' '{config[contig_name]}' '{params.different_contig_name}'
        else
            touch '{output}'
        fi
        """

rule sort_flank_bed:
    input:
        "tes/{te}/14_{flank}_flank.bed_unsorted"
    output:
        "tes/{te}/14_{flank}_flank.bed"
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

rule closest_flank_maps_to_ref:
    input:
        flank_5p = "tes/{te}/14_5p_flank.bed",
        flank_3p = "tes/{te}/14_3p_flank.bed"
    output:
        "tes/{te}/15_flank_overlap.bed"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input.flank_5p}' ] && [ -s '{input.flank_3p}' ]; then
            bedtools closest -a '{input.flank_5p}' -b '{input.flank_3p}' -s -d -t all > '{output}'
        else
            touch '{output}'
        fi
        """

checkpoint json_for_report:
    input:
        overlap = "tes/{te}/15_flank_overlap.bed",
        info_5p = "tes/{te}/14_5p_flank.info",
        info_3p = "tes/{te}/14_3p_flank.info"
    output:
        "tes/{te}/15_flank_overlap.json"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_liftover]} bed_to_json {input.overlap} {input.info_5p} {input.info_3p} {output} || true
        touch {output}
        """

rule make_report:
    input:
        overlap = "tes/{te}/15_flank_overlap.json",
        te_json = "tes/{te}/00_annotation.json",
        ref_bed = lambda wildcards: f"../../ref_repeatmask/{os.path.basename(config['reference'])}.te.bed",
        reference = config['reference']
    output:
        "tes/{te}/16_{overlap_id}_report.json"
    params: 
        flank_overlap_max = config["overlap"],
        flank_gap_max = config["gap"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_liftover]} make_report {input.overlap} {wildcards.overlap_id} {input.te_json} {input.ref_bed} {input.reference} {params.flank_overlap_max} {params.flank_gap_max} {output} || true
        touch {output}
        """

def overlap_ids_report(wildcards):
    overlap_file = checkpoints.json_for_report.get(**wildcards).output[0]
    try:
        with open(overlap_file, "r") as overlap:
            overlap_dict = json.load(overlap)
            overlap_ids = [key for key in overlap_dict]
        return [f"tes/{wildcards.te}/16_{overlap_id}_report.json" for overlap_id in overlap_ids]
    except:
        print(traceback.format_exc())
        return []
rule best_report:
    input:
        flanks = [
            "tes/{te}/14_5p_flank.info",
            "tes/{te}/14_3p_flank.info",
            "tes/{te}/14_5p_flank.bed",
            "tes/{te}/14_3p_flank.bed"
            ],
        ref_bed = lambda wildcards: f"../../ref_repeatmask/{os.path.basename(config['reference'])}.te.bed",
        te_json = "tes/{te}/00_annotation.json",
        overlap_reports = overlap_ids_report
    output:
        "tes/{te}/17_best_report.json"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_liftover]} choose_report {output} {input}"


'''4th stage
Read extraction (samtools)
'''

rule read_context:
    input:
        vcf_parsed = "00_vcf_parsed.tsv",
        bam = "../../reads_sort.bam"
    output:
        read_ids = "00_read_context.id",
        vcf_parsed_new = "00_parsed_vcf_with_readcount.tsv"
    params:
        window = 1000
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_assembly]} read_context '{config[contig_name]}' '{input.vcf_parsed}' '{input.bam}' '{output.read_ids}' '{output.vcf_parsed_new}' '{params.window}'"

"""
4th stage: estimate intra-sample TE allele frequency (TAF)
"""

'''4th stage
Read alignment to TE contig (minimap2)
'''

rule get_reverse_complement:
    input:
        "03_contig1.fa"
    output:
        "10_revcomp.fa"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_utility]} get_rev_comp_sequence '{input}' '{output}'
        """
        
rule realignment:
    input:
        contig = "{contig_revcomp}.fa",
        reads = "00_read_context.fa"#TODO check this
    output:
        "{contig_revcomp}_realign.sam"
    params:
        presets = lambda wildcards: {"pacbio":"map-pb","ont":"map-ont"}[config["presets"]]
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        minimap2 -t {threads} -a -x '{params.presets}' -v 0 '{input.contig}' '{input.reads}' > '{output}'
        """

rule realignment_to_bam:
    input:
        "{contig_revcomp}_realign.sam"
    output:
        "{contig_revcomp}_realign.bam"
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input}' ]; then
            samtools view -@ {threads} -bS '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule sort_index_realignment:
    input:
        "{contig_revcomp}_realign.bam"
    output:
        "{contig_revcomp}_realign.sort.bam"
    threads: config["threads"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input}' ]; then
            samtools sort -@ '{threads}' -o '{output}' '{input}'
            samtools index -@ '{threads}' '{output}'
        else
            touch '{output}'
        fi
        """


'''4th stage
Depth-based TAF estimation (SAMtools)
'''

rule estimate_te_depth:
    input:
        bam = "{contig_revcomp}_realign.sort.bam",
        contig = "{contig_revcomp}.fa"
    output:
        depth_5p = "tes/{te}/{contig_revcomp}_5p_te.depth",
        depth_3p = "tes/{te}/{contig_revcomp}_3p_te.depth"
    params:
        te_interval = config["af_te_interval"],
        te_offset = config["af_te_offset"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input.bam}' ]; then
            python3 {config[STELR_te]} estimate_te_depth '{input.bam}' '{input.contig}' '{wildcards.te}' '{params.te_interval}' '{params.te_offset}' '{output.depth_5p}' '{output.depth_3p}'
        else
            touch '{output}'
        fi
        """

rule estimate_flank_depth:
    input:
        bam = "{contig_revcomp}_realign.sort.bam",
        contig = "{contig_revcomp}.fa"
    output:
        depth_5p = "tes/{te}/{contig_revcomp}_5p_flank.depth",
        depth_3p = "tes/{te}/{contig_revcomp}_3p_flank.depth"
    params:
        flank_len = config["af_flank_interval"],
        flank_offset = config["af_flank_offset"]
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        if [ -s '{input.bam}' ]; then
            python3 {config[STELR_te]} estimate_flank_depth '{input.bam}' '{input.contig}' '{wildcards.te}' '{params.flank_len}' '{params.flank_offset}' '{output.depth_5p}' '{output.depth_3p}'
        else
            touch '{output}'
        fi
        """
    
rule estimate_coverage:
    input:
        te_5p = "tes/{te}/{contig_revcomp}_5p_te.depth",
        te_3p = "tes/{te}/{contig_revcomp}_3p_te.depth",
        flank_5p = "tes/{te}/{contig_revcomp}_5p_flank.depth",
        flank_3p = "tes/{te}/{contig_revcomp}_3p_flank.depth"
    output:
        "tes/{te}/{contig_revcomp}.freq"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_te]} estimate_coverage '{input.te_5p}' '{input.te_3p}' '{input.flank_5p}' '{input.flank_3p}' '{output}'"

rule get_allele_frequency:
    input:
        fwd = "tes/{te}/03_contig1.freq",
        rev = "tes/{te}/10_revcomp.freq"
    output:
        "tes/{te}/11_allele_frequency.json"
    conda:
        config["conda"]["stable_environment"]
    shell:
        "python3 {config[STELR_te]} get_af '{input.fwd}' '{input.rev}' '{output}'"

'''
Write Output
'''

rule individual_json:
    input:
        liftover_file = "tes/{te}/17_best_report.json",
        af_file = "tes/{te}/11_allele_frequency.json",
        vcf_parsed = "00_vcf_parsed.tsv",
        annotation_file = "tes/{te}/00_annotation.bed",
        contig_file = "03_contig1.fa"
    output:
        "tes/{te}/18_output.json"
    conda:
        config["conda"]["stable_environment"]
    shell:
        """
        python3 {config[STELR_output]} make_json_output {input} {output} {config[contig_name]}
        """

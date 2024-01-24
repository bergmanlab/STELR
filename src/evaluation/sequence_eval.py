import os
import glob
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path
import subprocess
import json
import traceback
from binf_util import paf_file, bed_file, minimap2

def make_stelr_te_fasta(input,output):
    with open(input,"r") as contig_jsons:
        te_dict = json.load(contig_jsons)
    with open(output,"w") as output_file:
        output_file.write(f">{te_dict['contig_id']}\n{te_dict['te_sequence']}")

def find_corresponding_annotation(telr_json, input_paf, community_annotation, ref, bed_out, fasta_out):
    with open(telr_json,"r") as telr_json:
        telr_data = json.load(telr_json)

    ref_paf = paf_file(input_paf)
    ref_align_start = ref_paf.get("start")
    ref_align_end = ref_paf.get("end")

    annotation_data = bed_file(community_annotation).data
    
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
        bed_file(corresponding_annotation).write_out(bed_out)
        coord = f"{corresponding_annotation[0]}:{corresponding_annotation[1]+1}-{corresponding_annotation[2]}"
        with open(fasta_out,"w") as output_file:
            subprocess.run(["samtools","faidx",ref,coord], stdout=output_file)
    else:
        if len(overlapping_annotations) > 1:
            print(f"{wildcards.contig}: more than one TE in the aligned region")
        else: print(f"{wildcards.contig}: no ref TEs in the aligned region"):
        with open(output[0],"w") as output_file: output_file.write("")
        with open(output[1],"w") as output_file: output_file.write("")

def sequence_eval_report(flank_len, ref_paf, annotation_bed, annotation_paf, telr_json, output_json)
    #parse params
    flank_len = int(flank_len)

    #parse input files
    ref_paf = paf_file(ref_paf)
    annotation_bed = bed_file(annotation_bed)
    annotation_paf = paf_file(annotation_paf)
    paftools_summary = annotation_paf.paftools_summary()
    with open(telr_json,"r") as input_file:
        json_data = json.load(input_file)

    #compile report
    eval_report = {
        "contig_te_plus_flank_start": max(0,json_data["contig_te_start"] - flank_len),
        "contig_te_plus_flank_end": min(json_data["contig_length"], json_data["contig_te_end"] + flank_len),
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

    #write output
    with open(output_json,"w") as output_file:
        json.dump(eval_report, output_file)

def evaluate_sequence(
        telr_json,
        ref_fasta,
        community_annotation,
        flank_len,
        report_file
        stelr="stelr"):

    keep_intermediates = False
    stelr = stelr.lower()

    try:
        telr_data = json.loads(telr_json)
        telr_data["family"]
    except:
        telr_data = json.load(telr_json)

    prediction_head = f"{telr_data['ID']}.{stelr}"
    predicted_te_fasta = f"{prediction_head}.fasta"
    with open(predicted_te_fasta,"w") as te_fasta_file:
        te_fasta_file.write(f">{telr_data['ID']}\n{telr_data['te_sequence']}")

    if keep_intermediates:
        prediction_paf = f"{prediction_head}-to-{'.'.join(ref_fasta.split('.')[:-1])}.paf"
        minimap2(ref_fasta,predicted_te_fasta,"cigar",output=prediction_paf)
        prediction_paf = paf_file(prediction_paf)
    else:
        prediction_paf = paf_file(minimap2(ref_fasta,predicted_te_fasta,"cigar"))

    ref_align_start = prediction_paf.get("start")
    ref_align_end = prediction_paf.get("end")

    annotation_data = bed_file(community_annotation).data
    
    aligned_chrom = prediction_paf.longest[5]
    candidate_annotations = [line for line in annotation_data if line[0] == aligned_chrom and line[3] == telr_data["family"]]
    location_filters = [
        lambda start, end: start >= ref_align_start and end <= ref_align_end,
        lambda start, end: start <= ref_align_end and end >= ref_align_end,
        lambda start, end: start <= ref_align_start and end >= ref_align_start
    ]
    overlapping_annotations = [line for line in candidate_annotations if any([f(line[1],line[2]) for f in location_filters])]

    if len(overlapping_annotations) != 1:
        with open(report_file,"w") as outfile:
            if len(overlapping_annotations) > 1:
                outfile.write(f'"{wildcards.contig}: more than one TE in the aligned region"')
            else: outfile.write(f'"{wildcards.contig}: no ref TEs in the aligned region"')
        quit()
    

    annotation = bed_file(corresponding_annotations[0])
    annotation_head = f"{prediction_head}-corresponding_annotation"
    annotation.faidx(f"{annotation_file_head}.fasta")

    minimap2(ref_fasta,annotation.fasta,"cs",preset="asm5",output=f"{annotation_head}_unsorted.paf")
    process(["sort","-k6,6","-k8,8n",f"{annotation_head}_unsorted.paf"]).write_to(f"{annotation_head}.paf")
    annotation_paf = paf_file(f"{annotation_head}.paf")
    paftools_summary = annotation_paf.paftools_summary()
    
    #parse params
    flank_len = int(flank_len)

    #compile report
    eval_report = {
        "contig_te_plus_flank_start": max(0,telr_data["contig_te_start"] - flank_len),
        "contig_te_plus_flank_end": min(telr_data["contig_length"], telr_data["contig_te_end"] + flank_len),
        "contig_te_plus_flank_size": None,#calculated below
        "num_contig_ref_hits": prediction_paf.count(),
        "ref_aligned_chrom": prediction_paf.get("chrom"),
        "ref_aligned_start": prediction_paf.get("start"),
        "ref_aligned_end": prediction_paf.get("end"),
        "contig_max_base_mapped_prop": prediction_paf.get("map_prop"),
        "contig_mapp_qual": prediction_paf.get("map_qual"),
        "contig_num_residue_matches": prediction_paf.get("matches"),
        "contig_alignment_block_length": prediction_paf.get("align_len"),
        "contig_blast_identity": prediction_paf.get("blast_id"),
        "ref_te_family": annotation.get("name"),
        "ref_te_length": annotation.get("end") - annotation.get("start"),
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

    #write output
    with open(output_json,"w") as output_file:
        json.dump(eval_report, output_file)


if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
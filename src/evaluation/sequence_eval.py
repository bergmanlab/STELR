import os
import sys
import json
from binf_util import paf_file, bed_file, minimap2, process

def evaluate_sequence(
        telr_json,
        ref_fasta,
        community_annotation,
        report_file,
        flank_len=500,
        stelr="stelr",
        keep_intermediates=False):

    intermediate_files = []
    keep_intermediates = keep_intermediates and not str(keep_intermediates).lower() in ["false","f"]
    flank_len = int(flank_len)
    stelr = stelr.lower()

    # load TELR or STELR json file
    try:
        telr_data = json.loads(telr_json)
        telr_data["family"]
    except:
        with open(telr_json,"r") as telr_json:
            telr_data = json.load(telr_json)

    # Make a fasta file containing the TE sequence predicted by STELR
    prediction_head = f"{telr_data['ID']}.{stelr}"
    predicted_te_fasta = f"{prediction_head}.fasta"
    with open(predicted_te_fasta,"w") as te_fasta_file:
        te_fasta_file.write(f">{telr_data['ID']}\n{telr_data['te_sequence']}")
        intermediate_files.append(predicted_te_fasta)

    # align the STELR-predicted TE sequence to the community reference genome
    if keep_intermediates:
        prediction_paf = f"{prediction_head}-to-{'.'.join(ref_fasta.split('.')[:-1])}.paf"
        minimap2(ref_fasta,predicted_te_fasta,"cigar",output=prediction_paf)
        prediction_paf = paf_file(prediction_paf)
    else:
        prediction_paf = paf_file(minimap2(ref_fasta,predicted_te_fasta,"cigar"))

    # extract the start and end of the longest alignment
    ref_align_start = prediction_paf.get("start")
    ref_align_end = prediction_paf.get("end")

    # read the community reference annotation into a table
    annotation_data = bed_file(community_annotation).data
    
    # filter all of the annotated TEs in the reference genome by whether they overlap with the position where the predicted TE aligned to the reference
    aligned_chrom = prediction_paf.longest[5]
    candidate_annotations = [line for line in annotation_data if line[0] == aligned_chrom and line[3] == telr_data["family"]]
    location_filters = [
        lambda start, end: start >= ref_align_start and end <= ref_align_end,
        lambda start, end: start <= ref_align_end and end >= ref_align_end,
        lambda start, end: start <= ref_align_start and end >= ref_align_start
    ]
    overlapping_annotations = [line for line in candidate_annotations if any([f(line[1],line[2]) for f in location_filters])]

    # check if there was exactly 1 TE annotation at the site of alignment
    if len(overlapping_annotations) != 1:
        with open(report_file,"w") as outfile:
            if len(overlapping_annotations) > 1:
                outfile.write(f'"{wildcards.contig}: more than one TE in the aligned region"')
            else: outfile.write(f'"{wildcards.contig}: no ref TEs in the aligned region"')
        quit()
    
    # extract the reference genome sequence at the site of alignment
    annotation = bed_file(overlapping_annotations[0])
    annotation_head = f"{prediction_head}-corresponding_annotation"
    annotation.faidx(ref_fasta,f"{annotation_head}.fasta")

    # use minimap2 to align STELR's predicted sequence for the TE to its actual reference sequence
    annotation_paf, unsorted_annotation_paf = f"{annotation_head}.paf", f"{annotation_head}_unsorted.paf"
    minimap2(ref_fasta,annotation.fasta,"cs",preset="asm5",output=unsorted_annotation_paf)
    intermediate_files.append(unsorted_annotation_paf)
    process(["sort","-k6,6","-k8,8n",unsorted_annotation_paf]).write_to(annotation_paf)
    intermediate_files.append(annotation_paf)
    annotation_paf = paf_file(annotation_paf)
    paftools_summary = annotation_paf.paftools_summary()
    
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

    # write output
    with open(report_file,"w") as output_file:
        json.dump(eval_report, output_file, indent=4)
    
    # remove intermediate files
    if not keep_intermediates:
        for file in intermediate_files:
            os.remove(file)


if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
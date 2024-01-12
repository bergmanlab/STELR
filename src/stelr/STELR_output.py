import os
import sys
import pandas as pd
import logging
import json
from Bio import SeqIO
from datetime import date
import subprocess
import glob
from STELR_utility import check_exist, abs_path
from multiprocessing import Pool


def get_contig_info(reference_index):
    contig_info = []
    with open(reference_index, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_info.append("##contig=<ID={},length={}>".format(entry[0], entry[1]))
    return contig_info

def make_json_output(liftover_file, af_file, vcf_parsed_file, annotation_file, contig_file, json_output, contig_name):
    with open(liftover_file, "r") as data:
        liftover_report = json.load(data)
    if not "type" in liftover_report:
        quit()
    if not liftover_report["type"] == "non-reference":
        quit()
    with open(af_file, "r") as data:
        af_dict = json.load(data)
    with open(vcf_parsed_file, "r") as data:
        sniffles_info = [line for line in data][0].replace("\n","").split("\t")
        sniffles_info = {
            "genotype":sniffles_info[10],
            "alt_count":sniffles_info[12],
            "ref_count":sniffles_info[11]
        }
    sequence = subprocess.run(f"bedtools getfasta -fi '{contig_file}' -bed '{annotation_file}' -s", shell=True, capture_output=True, text=True).stdout.split("\n")[1]
    te_name = liftover_file.split("tes/")[1].split("/")[0]
    te_start, te_end = [int(index) for index in te_name.replace("te_","").split("_")]
    contig_length = 0
    with open(contig_file, "r") as data:
        next(data)
        for line in data:
            contig_length += len(line.replace("\n",""))

    full_report = {
        "type": "non-reference",
        "ID": f"{liftover_report['chrom']}_{liftover_report['start']}_{liftover_report['end']}_{liftover_report['family']}",
        "chrom": liftover_report["chrom"],
        "start": liftover_report["start"],
        "end": liftover_report["end"],
        "family": liftover_report["family"],
        "strand": liftover_report["strand"],
        "support": "single_side",#will be processed further
        "tsd_length": liftover_report["TSD_length"],
        "tsd_sequence": liftover_report["TSD_sequence"],#will be processed further
        "te_sequence": sequence,
        "genotype": sniffles_info["genotype"],
        "num_sv_reads": sniffles_info["alt_count"],
        "num_ref_reads": sniffles_info["ref_count"],
        "allele_frequency": af_dict["freq"],
        "gap_between_flank": liftover_report["gap"],
        "te_length": len(sequence),
        "contig_id": contig_name,
        "contig_length": contig_length,
        "contig_te_start": te_start,
        "contig_te_end": te_end,
        "5p_flank_align_coord": liftover_report["5p_flank_align_coord"],
        "5p_flank_mapping_quality": liftover_report["5p_flank_mapping_quality"],
        "5p_flank_num_residue_matches": liftover_report["5p_flank_num_residue_matches"],
        "5p_flank_alignment_block_length": liftover_report["5p_flank_alignment_block_length"],
        "5p_flank_sequence_identity": liftover_report["5p_flank_sequence_identity"],
        "3p_flank_align_coord": liftover_report["3p_flank_align_coord"],
        "3p_flank_mapping_quality": liftover_report["3p_flank_mapping_quality"],
        "3p_flank_num_residue_matches": liftover_report["3p_flank_num_residue_matches"],
        "3p_flank_alignment_block_length": liftover_report["3p_flank_alignment_block_length"],
        "3p_flank_sequence_identity": liftover_report["3p_flank_sequence_identity"],
        "te_5p_cov": af_dict["fwd"]["5p"]["te"],
        "te_3p_cov": af_dict["fwd"]["3p"]["te"],
        "flank_5p_cov": af_dict["fwd"]["5p"]["flank"],
        "flank_3p_cov": af_dict["fwd"]["3p"]["flank"],
        "te_5p_cov_rc": af_dict["rev"]["5p"]["te"],
        "te_3p_cov_rc": af_dict["rev"]["3p"]["te"],
        "flank_5p_cov_rc": af_dict["rev"]["5p"]["flank"],
        "flank_3p_cov_rc": af_dict["rev"]["3p"]["flank"]
    }

    if full_report["tsd_sequence"]:
        full_report["tsd_sequence"] = full_report["tsd_sequence"].upper()
    
    if full_report["5p_flank_align_coord"] and full_report["3p_flank_align_coord"]:
        full_report["support"] = "both_sides"

    with open(json_output, "w") as output:
        json.dump(full_report, output)


def get_te_data(contig):
    #basic_report_keys = ["type","ID","chrom","start","end","family","strand","support","tsd_length","tsd_sequence","te_sequence","genotype","num_sv_reads","num_ref_reads","allele_frequency"]
    try:
        tes = next(os.walk(f"contigs/{contig}/tes"))[1]
        te_data = []
        contig_path = ""
        for te in tes:
            try:
                with open(f"contigs/{contig}/tes/{te}/18_output.json","r") as data:
                    te_dict = json.load(data)
                    #te_dict = {
                    #    "expanded_json":te_info,
                    #    "json":{key:te_info[key] for key in basic_report_keys},
                    #    "te_fasta":f">{te_info['ID']}\n{te_info['te_sequence']}\n",
                    #    "bed_out":f"{te_info['chrom']}\t{te_info['start']}\t{te_info['end']}\t{te_info['family']}\t.\t{te_info['strand']}\n"
                    #}
                    te_data.append(te_dict)
            except:pass
        if len(te_data) == 0: return []
        return [contig,te_data]
    except: return []

'''
python3 $STELR_output compile_te_data all_tes.json 10
'''
def compile_te_data(outfile, threads=1):
    threads = int(threads)
    contigs = next(os.walk("contigs"))[1]
    with Pool(threads) as p:
        compiled_data = [item for item in p.map(get_te_data,contigs) if item]
        compiled_data.sort()
    with open(outfile,"w") as out:
        json.dump(compiled_data,out)

'''
time python3 $STELR_output write_contig_fasta_output all_tes.json reads.stelr.loci.fasta 10
'''
def return_content(file):
    with open(file,"r") as input_file:
        return [file,input_file.read()]
def write_contig_fasta_output(all_tes,out_fasta,threads=1,contig_pattern="03_contig1.fa"):
    threads = int(threads)

    with open(all_tes,"r") as data:
        te_data = json.load(data)
    
    contig_fastas = [f"contigs/{contig[0]}/{contig_pattern}" for contig in te_data]

    with Pool(threads) as p:
        contig_fastas = p.map(return_content,contig_fastas)
        contig_fastas.sort()

    with open(out_fasta,"w") as output:
        output.write("".join([item[1] for item in contig_fastas]))

'''
time python3 $STELR_output write_te_fasta_output all_tes.json reads.stelr.te.fasta
'''
def write_te_fasta_output(all_tes,out_fasta):
    with open(all_tes,"r") as data:
        te_data = json.load(data)
    
    te_fastas = "\n".join(["\n".join([f">{te_info['ID']}\n{te_info['te_sequence']}" for te_info in te_list[1]]) for te_list in te_data])

    with open(out_fasta,"w") as output:
        output.write(te_fastas)

def write_bed_output(all_tes,out_bed):
    with open(all_tes,"r") as data:
        te_data = json.load(data)
    
    te_beds = "\n".join(["\n".join([f"{te_info['chrom']}\t{te_info['start']}\t{te_info['end']}\t{te_info['family']}\t.\t{te_info['strand']}" for te_info in te_list[1]]) for te_list in te_data])
    
    with open(out_bed,"w") as output:
        output.write(te_beds)

'''
time python3 $STELR_output write_te_json all_tes.json reads.stelr.te.json
'''
def write_te_json(all_tes,output_file):
    with open(all_tes,"r") as data:
        te_data = json.load(data)
    
    expanded_json = [te_dict for te_list in te_data for te_dict in [te_dict for te_dict in te_list[1]]]

    with open(output_file,"w") as output:
        json.dump(expanded_json,output,indent=4)

    
    #basic_report_keys = ["type","ID","chrom","start","end","family","strand","support","tsd_length","tsd_sequence","te_sequence","genotype","num_sv_reads","num_ref_reads","allele_frequency"]

    #basic_json = [{key:report[key] for key in basic_report_keys} for report in expanded_json]
    
    #with open(basic_json_file,"w") as output:
    #    json.dump(basic_json,output)


def get_locus_info(locus):
    locus = locus.split("\t")
    contig_id = "_".join(locus[:3])
    locus[1], locus[2] = int(locus[1]), int(locus[2])
    info = {
        "ID":contig_id,
        "chrom":locus[0],
        "start":locus[1],
        "end":locus[2],
        "num_supporting_reads":0,
        "supporting_reads":locus[8].split(","),
        "passed_te_candidate_filter":False,
        "contig_assembled":False,
        "TE_predicted":False,
        "assembled_contig_length":None,
        "TEs":[]
    }
    info["num_supporting_reads"] = len(info["supporting_reads"])

    if os.path.isdir(f"contigs/{contig_id}"):
        info["passed_te_candidate_filter"] = True
    else: return info
    
    if check_exist(f"contigs/{contig_id}/03_contig1.fa"):
        info["contig_assembled"] = True
        with open(f"contigs/{contig_id}/03_contig1.fa","r") as contig_fa:
            info["assembled_contig_length"] = int(contig_fa.readline().strip().split("len=")[1])
    else: return info
    
    tes = glob.glob(f"contigs/{contig_id}/tes/*/18_output.json")
    for te in tes:
        try:
            with open(te,"r") as input_file:
                info["TEs"].append(input_file.read().split('"ID":')[1].split('"')[1])
        except: pass
    
    num_tes = len(info["TEs"])
    if num_tes > 0: info["TE_predicted"] = True
    
    return info
        

'''
time python3 $STELR_output write_contig_json reads.vcf_parsed.tsv reads.stelr.loci.json 10
'''
def write_contig_json(parsed_svs, output_file, threads=1):
    threads = int(threads)
    with open(parsed_svs,"r") as input_file:
        parsed_svs_list = input_file.read().strip().split("\n")
    
    with Pool(threads) as p:
        contig_info = p.map(get_locus_info,parsed_svs_list)
        contig_info.sort(key=lambda x:x["ID"])
    
    categorized_contig_info = {
        "no_TE_predicted":[c for c in contig_info if not c["passed_te_candidate_filter"]],
        "TE predicted but local assembly failed":[c for c in contig_info if c["passed_te_candidate_filter"] and not c["contig assembled?"]],
        "contig assembled, no TE found":[c for c in contig_info if c["contig assembled?"] and not c["TEs"]],
        "TE predicted at locus":[c for c in contig_info if c["TEs"]]
    }
    
    with open(output_file,"w") as output:
        json.dump(categorized_contig_info,output,indent=4)


'''
time python3 $STELR_output write_vcf_output reads.stelr.json input/reference.fasta input/reference.fasta.fai reads.stelr.vcf
'''
def write_vcf_output(basic_json, ref, ref_index, out_vcf):
    with open(basic_json,"r") as data:
        data = json.load(data)

    ref_info = get_contig_info(ref_index)
    df = pd.DataFrame(data)
    if not df.empty:
        df["ID"] = df.index
        df["start"] = df["start"] + 1
        df["REF"] = "N"
        df["QUAL"] = "."
        df["FILTER"] = "PASS"
        df["FORMAT"] = "GT:DR:DV"
        df["gt"] = df["genotype"] + ":" + df["num_sv_reads"] + ":" + df["num_ref_reads"]
        df["INFO"] = df.apply(
            lambda x: f"SVTYPE=INS;END={x.end};FAMILY={x.family};STRANDS={x.strand};SUPPORT_TYPE={x.support};RE={x.num_sv_reads};AF={x.allele_frequency};TSD_LEN={x.tsd_length};TSD_SEQ={x.tsd_sequence}",
            axis=1,
        )

        df = df[
            [
                "chrom",
                "start",
                "ID",
                "REF",
                "te_sequence",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "gt",
            ]
        ]
        df = df.fillna("NA")
    with open(out_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.1\n")
        vcf.write("##fileDate={}".format(date.today()) + "\n")
        vcf.write("##source=STELR\n")
        vcf.write(f"##reference={ref}\n")
        vcf.write("\n".join(ref_info) + "\n")
        vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structure variant">\n')
        vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structure variant">\n')
        vcf.write('##INFO=<ID=STRANDS,Number=A,Type=String,Description="Strand orientation">\n')
        vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        vcf.write('##INFO=<ID=FAMILY,Number=1,Type=String,Description="TE family">\n')
        vcf.write('##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">\n')
        vcf.write('##INFO=<ID=SUPPORT_TYPE,Number=1,Type=String,Description="single_side or both_sides">\n')
        vcf.write('##INFO=<ID=TSD_LEN,Number=1,Type=String,Description="Length of the TSD sequence if available">\n')
        vcf.write('##INFO=<ID=TSD_SEQ,Number=1,Type=String,Description="TSD sequence if available">\n')
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">\n')
        vcf.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        
    if not df.empty:
        df.to_csv(out_vcf, sep="\t", mode="a", index=False, header=False)


def write_output(contig_fa_outfile, te_fa_outfile, bed_outfile, json_outfile, expanded_json_outfile, vcf_outfile, reference, reference_index, output_pattern, threads=1):
    threads = int(threads)
    contigs = next(os.walk("contigs"))[1]
    with Pool(threads) as p:
        out = [item for item in p.map(get_te_data,contigs) if item]
        print(out)
    quit()
    json_files = [file for file in glob.glob(f"contigs/*/tes/*/{output_pattern}") if check_exist(file)]
    for file in json_files:
        try:
            with open(file, "r") as data:
                te_info = json.load(data)
                contig_name = te_info["contig_name"]
                te_name = te_info["json"]["ID"]
                if not contig_name in contigs:
                    contigs[contig_name] = {"contig_path":te_info["contig_path"],"te_list":{}}
                te_dict = {
                    "expanded_json":te_info["expanded_json"],
                    "json":te_info["json"],
                    "te_fasta":te_info["te_fasta"],
                    "bed_out":te_info["bed_out"]
                }
                contigs[contig_name]["te_list"][te_name] = te_dict
        except: pass
    
    json_output = []
    json_expanded_output = []
    if check_exist(contig_fa_outfile):
        os.remove(contig_fa_outfile)
    with open(te_fa_outfile, "w") as te_fa, open(bed_outfile, "w") as bed_out:
        for contig in contigs:
            subprocess.run(f"cat '{contigs[contig]['contig_path']}' >> {contig_fa_outfile}", shell=True)
            te_list = contigs[contig]["te_list"]
            for te in te_list:
                outputs = te_list[te]
                json_output.append(outputs["json"])
                json_expanded_output.append(outputs["expanded_json"])
                te_fa.write(outputs["te_fasta"])
                bed_out.write(outputs["bed_out"])
    with open(json_outfile, "w") as output:
        json.dump(json_output, output, indent=4, sort_keys=False)
    with open(expanded_json_outfile, "w") as output:
        json.dump(json_expanded_output, output, indent=4, sort_keys=False)

    write_vcf(json_output, reference, reference_index, vcf_outfile)


if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
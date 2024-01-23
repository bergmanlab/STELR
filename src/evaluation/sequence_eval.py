import os
import glob
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path
import subprocess
import json
import traceback

class paf_file:
    def __init__(file_path):
        self.file_path = file_path

        def format_paf_line(line):
            if line.strip():
                line = line.split("\t")
                for n in [1,7,8,9,10,11]:
                    line[n] = int(line[n])
                return line
            else:
                return None
        with open(file_path,"r") as paf_file:
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
    
    def paftools_summary():
        if not "paftools_fields" in self.__dict__:
            self.paftools_fields = {
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
            
            paftools_summary = subprocess.run(["paftools.js","call","-l","100","-L","100",self.file_path], capture_output=True, text=True).stderr

            for line in paftools_summary.split("\n"):
                try:
                    value = int(line.split()[0])
                    field = [field for field in self.paftools_fields if field in line][0]
                    if type(self.paftools_fields[field]) is int: self.paftools_fields[field] = value
                    else:
                        subfield = [subfield for subfield in self.paftools_fields[field] if subfield in line][0]
                        self.paftools_fields[field][subfield] = value
                        self.paftools_fields[field]["total"] += value
                except: pass
        
        return self.paftools_fields

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

def make_stelr_te_fasta(input,output):
    with open(input,"r") as contig_jsons:
        te_dict = json.load(contig_jsons)
    with open(output,"w") as output_file:
        output_file.write(f">{te_dict['contig_id']}\n{te_dict['te_sequence']}")

def find_corresponding_annotation(telr_json, input_paf, community_annotation, ref, bed_out, fasta_out):
    with open(telr_json,"r") as telr_json:
        telr_data = json.load(telr_json)

    #paf_data = get_paf_data(input_paf)
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


if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
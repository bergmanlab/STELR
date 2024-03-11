import os
import sys
import json
from binf_util import paf_file, bed_file, minimap2, process
from multiprocessing import Pool
from time import perf_counter
from shutil import rmtree
import traceback

def region_family_filter(filter_region, unfiltered_bed, filtered_bed = None, unfiltered_json = None, filtered_json = None, exclude_families = ["INE_1"], exclude_nested_insertions = False):
    # Filter a bed file by intersection with the regular recombination region,
    # By excluding "nested" insertions ie insertions with multiple predicted families,
    # And by excluding insertions labelled as excluded families, ie "INE_1"
    
    # set up functions to filter file based on params
    exclude_families = set(exclude_families)
    filters = {}
    if exclude_nested_insertions:
        filters["nested_insertions"] = lambda te: "|" in te.split("\t")[3]
        if exclude_families:
            filters["exclude_families"] = lambda te: te.split("\t")[3] in exclude_families
    elif exclude_families:
        filters["exclude_families"] = lambda te: len(exclude_families.intersection(te.split("\t")[3].split("|"))) > 0
    
    # run bedtools intersect to find the intersection between the bed file and the regular recombination region
    intersection = process(["bedtools","intersect","-a",unfiltered_bed,"-b",filter_region,"-u"]).output(lines=True)

    # filter the intersection using functions set up earlier based on params.
    try:
        filtered_annotation = [line for line in intersection if not any([filters[param](line) for param in filters])]
    except:
        print(intersection)
        sys.exit(1)#TODO

    filtered_annotation = "\n".join(filtered_annotation)

    # write output
    if filtered_bed:
        with open(filtered_bed, "w") as out:
            out.write(filtered_annotation)
    
    if filtered_json:
        with open(unfiltered_json,"r") as input_file:
            json_data = json.load(input_file)
        
        filtered_json_data = [entry for entry in json_data if entry["ID"].replace("_","\t") in filtered_annotation]
        
        with open(filtered_json,"w") as output_file:
            json.dump(filtered_json_data,output_file,indent=4)

def evaluate_family_and_position(
        stelr_json,
        filtered_annotation,
        summary_file,
        exclude_nested_insertions = False,
        relax_mode = False,
        family_match_required=True,
        window=5,
        stelr="stelr",
        quality_checks=["allele_frequency","support"]):
    stelr = stelr.upper()

    # Load the expanded output json from (s)stelr
    with open(stelr_json,"r") as input_file:
        stelr_json = json.load(input_file)
    # Key each TE's output dict to its bed format string
    te_dict = {te_dict["ID"]:te_dict for te_dict in stelr_json}

    # Define the quality checks TEs must pass to pass quality control
    defined_quality_checks = {                                  # TEs pass each quality control step if:
        "allele_frequency": lambda value: value is not None,    # they have a predicted value for allele frequency
        "support": lambda value: value == "both_sides"          # they have both-sided support
    }
    quality_checks = {key:defined_quality_checks[key] for key in quality_checks}
    # set conditions by which two sets of families are considered a match based on params
    if family_match_required:
        if exclude_nested_insertions: # only one family allowed per TE
            family_match = lambda overlap: overlap[3] == overlap[9]
        elif relax_mode: # match as long as any TE family label is present in both sets
            family_match = lambda overlap: len(set(overlap[3].split("|")).intersection(overlap[9].split("|"))) > 0
        else: # match only if all TE labels in each set are present in both
            family_match = lambda overlap: set(overlap[3].split("|")) == set(overlap[9].split("|"))
    else: family_match = lambda _: True

    # make a list of TEs which failed because of each quality check
    failed_tes = {value:[te for te in te_dict if not quality_checks[value](te_dict[te][value])] for value in quality_checks}
    #print(failed_tes)
    # the filtered TEs are those which are not present in any of the quality check failure lists.

    filtered_predictions = bed_file([te_dict[te] for te in te_dict if not any([(lambda x: te in failed_tes[x])(value) for value in quality_checks])])
    filtered_predictions.write_out("test_intermediate.bed")
    filtered_annotation = bed_file(filtered_annotation)

    # create a list of true positives
    true_positives = bed_file([te[:6] for te in filtered_predictions.intersect(filtered_annotation,window=window) if family_match(te)])
    false_negatives = filtered_annotation.exclude(filtered_predictions,window=window)

    # count and print the total number of TELR predictions that passed region, family, and quality control checks
    total_filtered_predictions = len(filtered_predictions)
    # calculate the # of true positives and false positives
    summary_dict = {"total predictions":total_filtered_predictions,"true positives":len(true_positives)}
    summary_dict["false positives"] = total_filtered_predictions - summary_dict["true positives"]
    # count the number of lines in the bed file containing the false negatives
    summary_dict["false negatives"] = len(false_negatives)
    # calculate the precision and recall
    summary_dict["precision"] = round(summary_dict["true positives"]/total_filtered_predictions, 3)
    summary_dict["recall"] = round(summary_dict["true positives"]/len(filtered_annotation), 3)

    # print summary statistcs and write them to output file
    print(f"Total {stelr} Predictions (filtered): {total_filtered_predictions}")    
    for stat in ["true positives","false positives","false negatives"]:
        print(f"Number of {stelr} {stat}: {summary_dict[stat]}")        
    with open(summary_file,"w") as output_file:
        json.dump(summary_dict,output_file, indent=4)

def single_seq_eval(args):
    single_seq_report = args[4]
    if not os.path.isfile(single_seq_report):
        completed = evaluate_sequence(*args)
    else: completed = True
    if completed: 
        with open(single_seq_report,"r") as input:
            try:
                return json.load(input)
            except:
                print(f"sequence eval failed for {args[0]['ID']}", flush=True)
                quit(1)
    else: return None

def eval_i(stelr_json, stelr_contig_fa, ref_fasta, community_annotation, i=0, stelr="stelr"):
    i = int(i)
    with open(stelr_json,"r") as input:
        stelr_json = json.load(input)
    print(json.dumps([evaluate_sequence(stelr_json[i],stelr_contig_fa,ref_fasta,community_annotation,stelr=stelr,keep_intermediates=True)]))
#python3 $stelr_src/evaluation/eval_scripts.py eval_i 50x_diploid_homozygous/simulated_reads.telr.te_filtered.json community_reference.fasta annotation_liftover_filtered.bed 0 telr

def evaluate_sequence(
        stelr_json,
        stelr_contig_fa,
        ref_fasta,
        community_annotation,
        report_file = None,
        flank_len=500,
        stelr="stelr",
        keep_intermediates=False,
        same_family = True):
    start = perf_counter()
    intermediate_files = []
    stelr = stelr.lower()

    if report_file:
        out_dir = "/".join(report_file.split("/")[:-1])
    else: 
        out_dir = "."
    
    prediction_paf = paf_file()
    annotation = bed_file()
    annotation_paf = paf_file()

    def output(te_id, keep_intermediates, intermediate_files, output, report_file, start):
        if not keep_intermediates:
            for file in intermediate_files:
                os.remove(file)
        print(f"{te_id} processed in {perf_counter()-start} seconds",file=sys.stderr, flush=True)
        if not report_file: return output
        with open(report_file,"w") as outfile:
            json.dump(output,outfile,indent=4)
        return True
            

    # load TELR or STELR json file
    if type(stelr_json) is dict:
        stelr_data = stelr_json
    else:
        try:
            stelr_data = json.loads(stelr_json)
            stelr_data["family"]
        except:
            with open(stelr_json,"r") as stelr_json:
                stelr_data = json.load(stelr_json)
    te_id = stelr_data['ID']

    # Filter out records with single-sided support
    if not stelr_data["support"] == "both_sides": 
        print(f"{te_id} not processed due to single-sided support",file=sys.stderr, flush=True)
        return None

    # Make a fasta file containing the TE sequence predicted by STELR
    prediction_head = f"{out_dir}/{te_id}-flanks.{stelr}"
    te_flank_fasta = f"{prediction_head}.fasta"
    te_fasta = f"{out_dir}/{te_id}.{stelr}.fasta"

    if stelr_contig_fa.split("/")[-1] == "03_contig1.fa":
        with open(stelr_contig_fa,"r") as i:
            contig_seq = "".join(i.read().split("\n")[1:])
    else:
        with open(stelr_contig_fa,"r") as i:
            contig_seq = "".join(i.read().split(stelr_data["contig_id"])[1].split(">")[0].split("\n")[1:])    
    contig_seq = contig_seq[max(0,stelr_data["contig_te_start"]-flank_len):min(stelr_data["contig_length"],stelr_data["contig_te_end"]+flank_len)]

    with open(te_flank_fasta,"w") as o:
        o.write(f">{te_id}\n{contig_seq}")
        intermediate_files.append(te_flank_fasta)
    with open(te_fasta,"w") as o:
        o.write(f">{te_id}\n{stelr_data['te_sequence']}")
        intermediate_files.append(te_fasta)

    # align the STELR-predicted TE sequence to the community reference genome
    if keep_intermediates:
        prediction_paf = f"{prediction_head}-to-{'.'.join(ref_fasta.split('.')[:-1])}.paf"
        minimap2(ref_fasta,te_flank_fasta,"cigar",output=prediction_paf)
        prediction_paf = paf_file(prediction_paf)
    else:
        prediction_paf = paf_file(minimap2(ref_fasta,te_flank_fasta,"cigar"))

    if not prediction_paf:
        return output(te_id, keep_intermediates, intermediate_files, f'{te_id}: TE did not align to reference', report_file, start)
        
    # extract the start and end of the longest alignment
    ref_align_start = prediction_paf.get("start")
    ref_align_end = prediction_paf.get("end")

    # read the community reference annotation into a table
    annotation_data = bed_file(community_annotation).data
    
    # filter all of the annotated TEs in the reference genome by whether they overlap with the position where the predicted TE aligned to the reference
    if same_family: family_match = lambda x: x[3] == stelr_data["family"]
    else: family_match = lambda _: True
    aligned_chrom = prediction_paf.get("chrom")
    candidate_annotations = [line for line in annotation_data if line[0] == aligned_chrom and family_match(line)]
    location_filters = [
        lambda start, end: start >= ref_align_start and end <= ref_align_end,
        lambda start, end: start <= ref_align_end and end >= ref_align_end,
        lambda start, end: start <= ref_align_start and end >= ref_align_start
    ]
    overlapping_annotations = [line for line in candidate_annotations if any([f(line[1],line[2]) for f in location_filters])]

    # check if there was exactly 1 TE annotation at the site of alignment
    if len(overlapping_annotations) == 1:
        contig_cover_te = 1
    
        # extract the reference genome sequence at the site of alignment
        annotation = bed_file(overlapping_annotations[0])
        annotation_head = f"{prediction_head}-corresponding_annotation"
        annotation.faidx(ref_fasta,f"{annotation_head}.fasta")

        # use minimap2 to align STELR's predicted sequence for the TE to its actual reference sequence
        annotation_paf, unsorted_annotation_paf = f"{annotation_head}.paf", f"{annotation_head}_unsorted.paf"
        minimap2(annotation.fasta,te_fasta,"cs","cigar","secondary",preset="asm5",output=unsorted_annotation_paf)#changed, pay attention, TODO figure out if correct
        intermediate_files.append(unsorted_annotation_paf)
        process(["sort","-k6,6","-k8,8n",unsorted_annotation_paf]).write_to(annotation_paf)
        intermediate_files.append(annotation_paf)
        annotation_paf = paf_file(annotation_paf)    
    else:
        contig_cover_te = 0
        if len(overlapping_annotations) > 1:
            print(f'{te_id}: more than one TE in the aligned region',file=sys.stderr, flush=True)
        else:
            print(f'{te_id}: no ref TEs in the aligned region',file=sys.stderr, flush=True)
    
    paftools_summary = annotation_paf.paftools_summary()

    #compile report
    try:
        eval_report = {
            "ID":                               te_id,
            "contig_te_plus_flank_start":       max(0,stelr_data["contig_te_start"] - flank_len),
            "contig_te_plus_flank_end":         min(stelr_data["contig_length"], stelr_data["contig_te_end"] + flank_len),
            "contig_te_plus_flank_size":        None,#calculated below
            "num_contig_ref_hits":              prediction_paf.count(),
            "ref_aligned_chrom":                prediction_paf.get("chrom"),
            "ref_aligned_start":                prediction_paf.get("start"),
            "ref_aligned_end":                  prediction_paf.get("end"),
            "contig_max_base_mapped_prop":      prediction_paf.get("map_prop"),
            "contig_mapp_qual":                 prediction_paf.get("map_qual"),
            "contig_num_residue_matches":       prediction_paf.get("matches"),
            "contig_alignment_block_length":    prediction_paf.get("align_len"),
            "contig_blast_identity":            prediction_paf.get("blast_id"),
            "ref_te_family":                    annotation.get("name"),
            "ref_te_length":                    annotation.get("length"),
            "ref_te_num_1bp_del":               paftools_summary["deletions"]["1bp"],
            "ref_te_num_1bp_ins":               paftools_summary["insertions"]["1bp"],
            "ref_te_num_2bp_del":               paftools_summary["deletions"]["2bp"],
            "ref_te_num_2bp_ins":               paftools_summary["insertions"]["2bp"],
            "ref_te_num_50bp_del":              paftools_summary["deletions"]["[3,50)"],
            "ref_te_num_50bp_ins":              paftools_summary["insertions"]["[3,50)"],
            "ref_te_num_1kb_del":               paftools_summary["deletions"]["[50,1000)"],
            "ref_te_num_1kb_ins":               paftools_summary["insertions"]["[50,1000)"],
            "ref_te_num_ins":                   paftools_summary["insertions"]["total"],
            "ref_te_num_del":                   paftools_summary["deletions"]["total"],
            "ref_te_num_indel":                 annotation_paf.indels(),
            "contig_cover_te":                  contig_cover_te
        }
        eval_report["contig_te_plus_flank_size"] = eval_report["contig_te_plus_flank_end"] - eval_report["contig_te_plus_flank_start"]
    except:
        print(te_id)
        traceback.print_exc()
        sys.exit(1)

    return output(te_id, keep_intermediates, intermediate_files, eval_report, report_file, start)

def evaluate_all_sequences(
        stelr_json,
        stelr_fasta,
        ref_fasta,
        community_annotation,
        report_file,
        flank_len=500,
        stelr="stelr",
        threads=1,
        keep_intermediates = False):
    
    start = perf_counter()
    
    with open(stelr_json,"r") as input:
        stelr_json = json.load(input)
    
    destination_folder = "/".join(report_file.split("/")[:-1])
    tmp = destination_folder + "/seq_eval_tmp"
    if not os.path.exists(tmp): os.mkdir(tmp)
    else:
        import glob
        for file in glob.glob(f"{tmp}/*"):
            try:
                with open(file,"r") as input:
                    json.load(input)
            except:
                os.remove(file)
    
    eval_args = [[stelr_json[x], stelr_fasta, ref_fasta, community_annotation, f"{tmp}/{stelr_json[x]['ID']}", flank_len, stelr] for x in range(len(stelr_json))]

    with Pool(threads) as p:
        compiled_data = [item for item in p.map(single_seq_eval,eval_args) if item]
    
    with open(report_file,"w") as output:
        json.dump(compiled_data,output,indent=4)
    
    if not keep_intermediates:
        rmtree(tmp)
    
    print(f"STELR sequence evaluation completed in {perf_counter()-start}")
    
    
#For testing and debugging purposes, this rule can also be run on the command line as follows:
'''
python3 $stelr_src/evaluation/eval_scripts.py evaluate_all_sequences '{
    "stelr_json":"50x_diploid_homozygous/simulated_reads.telr.te_filtered.json",
    "ref_fasta":"community_reference.fasta",
    "community_annotation":"annotation_liftover_filtered.bed",
    "report_file":"50x_diploid_homozygous/all_seq_reports.telr.json",
    "flank_len":500,
    "stelr":"telr",
    "threads":10
}'
'''

def evaluate_zygosity(simulation, stelr_json, output_file):
    coverage, ploidy, genotype = simulation.split("_")

    target_bucket = {
        "diploid":{
            "homozygous":3,
            "heterozygous":1
        },
        "tetraploid":{
            "simplex":0,
            "duplex":1,
            "triplex":2,
            "quadruplex":3
        }
    }[ploidy][genotype]

    buckets = [[0,0.375],[0.375,0.625],[0.625,0.875],[0.875,1]]

    with open(stelr_json,"r") as input_file:
        stelr_json = json.load(input_file)
    
    zygosities = [d["allele_frequency"] for d in stelr_json if not d["allele_frequency"] is None]
    num_predicted = len(zygosities)
    num_not_predicted = len(stelr_json) - num_predicted
    
    buckets = [len([item for item in zygosities if item > bucket[0] and item <= bucket[1]]) for bucket in buckets]

    precision = buckets[target_bucket]/num_predicted

    if ploidy == "diploid":
        summary_dict = {
            "simulation": simulation,
            "coverage": coverage,
            "ploidy": ploidy,
            "actual_zygosity": genotype,
            "num_pred_total": num_predicted,
            "num_pred_hom": buckets[3],
            "num_pred_het": buckets[1],
            "num_pred_unclassified": buckets[0] + buckets[2],
            "precison": precision,
        }
    else:
        summary_dict = {
            "simulation": simulation,
            "coverage": coverage,
            "ploidy": ploidy,
            "actual_zygosity": genotype,
            "num_pred_total": num_predicted,
            "num_simplex": buckets[0],
            "num_duplex": buckets[1],
            "num_triplex": buckets[2],
            "num_quadruplex": buckets[3],
            "precison": precision,
        }
    
    with open(output_file,"w") as output:
        json.dump(summary_dict, output, indent=4)

def compare_seq_eval_results(oliver,shunhua, output_file=None):
    with open(oliver,"r") as input:
        odata = json.load(input)
    with open(shunhua,"r") as input:
        sdata = json.load(input)
    
    mismatches = {}
    num_not_mismatch = 0
    
    for n in range(len(odata)):
        if type(odata[n]) is str:
            pass
            #if sdata[n]["ref_te_family"] is not None:
            #    mismatches.append([odata[n],sdata[n]])
        else:
            match = [item for item in sdata if item["ID"] == odata[n]["ID"]]
            if match:
                match = {key:match[0][key] for key in match[0] if key in odata[n]}
                if not all([odata[n][x] == match[x] for x in odata[n]]):
                    mismatches[odata[n]["ID"]] = {key:[odata[n][key],match[key]] for key in odata[n] if odata[n][key] != match[key]}
                else: num_not_mismatch += 1
            else:
                print(f"No match for {odata[n]['ID']}")
    
    for n in range(len(sdata)):
        match = [item for item in odata if item["ID"] == sdata[n]["ID"]]
        if not match:
            print(f"Missed record for {sdata[n]['ID']}")
    
    for mismatch in mismatches:
        print(mismatch)
        print(mismatches[mismatch])
    print(num_not_mismatch)

    if output_file:
        with open(output_file,"w")as outfile:
            json.dump(mismatches,outfile,indent=4)

if __name__ == '__main__':
    args = []
    if len(sys.argv[2:]) == 1:
        globals()[sys.argv[1]](**json.loads(sys.argv[2]))
        quit()
    for arg in sys.argv[2:]:
        try:
            args.append(json.loads(arg))
        except:
            args.append(arg)
    globals()[sys.argv[1]](*args)
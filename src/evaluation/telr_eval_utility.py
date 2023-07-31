import os
import sys
import subprocess
from datetime import datetime, timedelta


def string2set(input_string, delimiter):
    string_set = set(input_string.split(delimiter))
    return string_set


def filter_family_bed(bed_in, family_filter, bed_out, method):
    """
    Filter BED file by including or excluding given TE families
    """
    with open(bed_out, "w") as output, open(bed_in, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            family_input_set = string2set(entry[3], delimiter="|")
            family_filter_set = string2set(family_filter, delimiter=",")
            if method == "include":
                if family_input_set.intersection(family_filter_set) == family_input_set:
                    output.write(line)
            else:
                if len(family_input_set.intersection(family_filter_set)) == 0:
                    output.write(line)


def filter_region_bed(bed_in, region_filter, bed_out):
    """
    Filter BED file by given regions
    """
    #with open(bed_out, "w") as output:
    output = subprocess.check_output(
        ["bedtools", "intersect", "-a", bed_in, "-b", region_filter, "-u"],
        text=True,
    )


def filter_annotation(
    bed_in,
    bed_out,
    filter_region,
    include_families,
    exclude_families,
    exclude_nested,
):
    """
    Filter BED file by given TE families and regions
    """
    # filter by region
    bed_filtered = bed_out + ".tmp"
    if filter_region is not None:
        filter_region_bed(
            bed_in=bed_in,
            region_filter=filter_region,
            bed_out=bed_filtered,
        )
    else:
        bed_filtered = bed_in

    # filter by included families
    if include_families is not None:
        bed_filtered_tmp = bed_out + ".tmp"
        filter_family_bed(
            bed_in=bed_filtered,
            family_filter=include_families,
            bed_out=bed_filtered_tmp,
            method="include",
        )
        os.rename(bed_filtered_tmp, bed_filtered)

    # filter by exluded families
    if exclude_families is not None:
        bed_filtered_tmp = bed_filtered + ".tmp"
        filter_family_bed(
            bed_in=bed_filtered,
            family_filter=exclude_families,
            bed_out=bed_filtered_tmp,
            method="exclude",
        )
        os.rename(bed_filtered_tmp, bed_filtered)

    # filter by nested families
    if exclude_nested is not None:
        bed_filtered_tmp = bed_filtered + ".tmp"
        filter_nested_bed(
            bed_in=bed_filtered,
            bed_out=bed_filtered_tmp,
        )
        os.rename(bed_filtered_tmp, bed_filtered)

    os.rename(bed_filtered, bed_out)


def filter_nested_bed(bed_in, bed_out):
    """
    Filter BED file by excluding nested TEs
    """
    with open(bed_out, "w") as output, open(bed_in, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if "|" not in entry[3]:
                output.write(line)


def count_lines(file):
    k = 0
    with open(file, "r") as input:
        for line in input:
            k = k + 1
    return k


def create_soft_link(input, out_dir):
    link = os.path.join(out_dir, os.path.basename(input))
    if os.path.islink(link):
        os.unlink(link)
    if not os.path.isabs(input):
        input = os.path.abspath(input)
    try:
        os.symlink(input, link)
    except Exception as e:
        print(e)
        sys.exit(1)
    return link


def check_exist(file):
    if os.path.isfile(file) and os.stat(file).st_size != 0:
        return True
    else:
        return False


def format_time(time):
    d = datetime(1, 1, 1) + timedelta(seconds=time)
    if d.hour == 0 and d.minute == 0:
        return "%d seconds" % (d.second)
    elif d.hour == 0 and d.minute != 0:
        return "%d minutes %d seconds" % (d.minute, d.second)
    else:
        return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)

def rename_chromosome(infile, output, chrom_dict="A4"):
    if chrom_dict == "A4":
        chrom_dict = {
            "CM010541.1":"chrX",
            "CM010542.1":"chr2L",
            "CM010543.1":"chr2R",
            "CM010544.1":"chr3L",
            "CM010545.1":"chr3R",
            "CM010546.1":"chr4"
        }
    elif not type(chrom_dict) is dict:
        with open(chrom_dict,"r") as chrom_json:
            chrom_dict = json.load(chrom_dict)
    f = open(infile,"r")
    bed_text = f.read()
    f.close()

    for key in chrom_dict:
        bed_text = bed_text.replace(key,chrom_dict[key])
    
    with open(output,"w") as out:
        out.write(bed_text)

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
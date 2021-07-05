import os
import sys
import subprocess
from datetime import datetime, timedelta


def filter_family_bed(bed_in, family_filter, bed_out, method):
    """
    Filter BED file by including or excluding given TE families
    """
    with open(bed_out, "w") as output, open(bed_in, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            family = entry[3].split("|")[0]
            if method == "include":
                if family in family_filter:
                    output.write(line)
            else:
                if family not in family_filter:
                    output.write(line)


def filter_region_bed(bed_in, region_filter, bed_out):
    """
    Filter BED file by given regions
    """
    with open(bed_out, "w") as output:
        subprocess.call(
            ["bedtools", "intersect", "-a", bed_in, "-b", region_filter, "-u"],
            stdout=output,
        )


def filter_annotation(
    bed_in, bed_out, filter_region, include_families, exclude_families
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

    os.rename(bed_filtered, bed_out)


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

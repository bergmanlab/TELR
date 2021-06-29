#!/usr/bin/env python3

__version__ = "1.1"
__author__ = "Shunhua Han"

import sys
import argparse
import os
import json
import subprocess
import shutil
from multiprocessing import Pool


def get_args(program_version, arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description="Program for lifting TE annotations from genome 1 to genome 2, this program requires SAMtools, bedtools and minimap2."
    )

    ## required ##
    parser.add_argument("--fasta1", type=str, help="genome 1", required=True)
    parser.add_argument("--fasta2", type=str, help="genome 2", required=True)
    parser.add_argument(
        "-1",
        "--bed1",
        type=str,
        help="annotation file for genome 1 in bed format",
        required=True,
    )
    parser.add_argument(
        "-2",
        "--bed2",
        type=str,
        help="annotation file for genome 2 in bed format",
        required=True,
    )

    ## optional ##
    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        help="minimap2 preset (default = 'asm10')",
        required=False,
    )
    parser.add_argument(
        "-l",
        "--flank_len",
        type=int,
        help="flanking sequence length (default = '500bp')",
        required=False,
    )
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        help="maximum distance between 5p and 3p flanking sequence alignments (default = '50bp')",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        help="output directory (default = 'Current working directory')",
        required=False,
    )
    parser.add_argument(
        "-k",
        "--keep_files",
        action="store_true",
        help="If provided then all intermediate files will be kept (default: remove intermediate files)",
        required=False,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="number of cores (default = '1')",
        required=False,
    )

    ## info ##
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )

    args = parser.parse_args(arguments)

    # checks if required files exist
    try:
        test = open(args.fasta1, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.fasta2, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.bed1, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.bed2, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    # set default values for optional arguments
    if args.window is None:
        args.window = 50

    if args.flank_len is None:
        args.flank_len = 500

    if args.preset is None:
        args.preset = "asm10"

    if args.threads is None:
        args.threads = 1

    return args


def get_ref_fa(ref, chrom, start, end, prefix, out_dir):
    """
    Extract subsequence from fasta based on coordinates
    """
    prefix_bed = "_".join([chrom, str(start), str(end)])
    bed = os.path.join(out_dir, prefix_bed + ".fa")
    subject = os.path.join(out_dir, prefix + ".fa")
    with open(bed, "w") as output:
        out_line = "\t".join([chrom, str(start), str(end)])
        output.write(out_line + "\n")
    with open(subject, "w") as output:
        subprocess.call(
            ["bedtools", "getfasta", "-fi", ref, "-bed", bed],
            stdout=output,
        )
    os.remove(bed)
    return subject


def paf_to_bed(paf, bed, filter=None):
    """
    convert PAF to sorted BED format
    """
    bed_tmp = bed + ".tmp"
    with open(paf, "r") as input, open(bed_tmp, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[5]
            if filter is not None:
                if chrom == filter:
                    start = entry[7]
                    end = entry[8]
                    name = entry[0]
                    score = entry[11]
                    strand = entry[4]
                    out_line = "\t".join([chrom, start, end, name, score, strand])
                    output.write(out_line + "\n")
            else:
                start = entry[7]
                end = entry[8]
                name = entry[0]
                score = entry[11]
                strand = entry[4]
                out_line = "\t".join([chrom, start, end, name, score, strand])
                output.write(out_line + "\n")

    if check_file(bed_tmp):
        with open(bed, "w") as output:
            subprocess.call(["bedtools", "sort", "-i", bed_tmp], stdout=output)
            os.remove(bed_tmp)


def align_flank(flank_fa, ref_fa, paf_out, preset, num_secondary):
    """
    Align flanking sequences to genome 2
    """
    with open(paf_out, "w") as output:
        subprocess.call(
            [
                "minimap2",
                "-cx",
                preset,
                "-v",
                "0",
                "-N",
                str(num_secondary),
                ref_fa,
                flank_fa,
            ],
            stdout=output,
        )


def get_coord(start_3p, end_3p, start_5p, end_5p, strand):
    if strand == "+":
        start = end_3p
        end = start_5p
    else:
        start = start_3p
        end = end_5p
    gap = end - start
    if start > end:
        start, end = end, start
    return start, end, gap


def create_bed(chrom, start, end, family, strand, filename):
    with open(filename, "w") as output:
        out_line = "\t".join([chrom, str(start), str(end), family, ".", strand])
        output.write(out_line)


def check_nearby_ref(
    chrom, start_query, end_query, family, strand, ref_bed, out_dir, threshold=5000
):
    """
    Check if flanking seqeunce alignments have nearby TE annotations in genome 2
    """
    prefix_query = "_".join([chrom, str(start_query), str(end_query), "query"])
    query_bed = os.path.join(out_dir, prefix_query + ".bed")
    create_bed(chrom, start_query, end_query, family, strand, query_bed)

    prefix_overlap = "_".join([chrom, str(start_query), str(end_query), "nearby"])
    overlap = os.path.join(out_dir, prefix_overlap + ".bed")
    with open(overlap, "w") as output:
        subprocess.call(
            [
                "bedtools",
                "closest",
                "-a",
                query_bed,
                "-b",
                ref_bed,
                "-d",
                "-D",
                "ref",
                "-k",
                "5"
                # "-s",
            ],
            stdout=output,
        )

    distance = None
    with open(overlap, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom1 = entry[0]
            chrom2 = entry[6]
            family1 = entry[3]
            family2 = entry[9]
            strand1 = entry[5]
            strand2 = entry[11]
            if chrom1 == chrom2 and family1 == family2 and strand1 == strand2:
                distance_new = int(entry[12])
                if distance is None:
                    distance = distance_new
                else:
                    distance = absmin(distance, distance_new)
    if distance is not None:
        if abs(distance) > threshold:
            distance = None
    return distance


def absmin(num1, num2):
    """
    Return the value with min aboslute
    """
    num1_abs = abs(num1)
    num2_abs = abs(num2)
    num_abs_min = min(num1_abs, num2_abs)
    if num_abs_min == num1_abs:
        return num1
    else:
        return num2


def run_liftover(input_json):
    """
    Main liftover function for a single annotation
    """
    # initiate liftover report dictionary
    lift_entries = dict()

    # load input data from json
    with open(input_json) as f:
        input_values = json.load(f)
    chrom = input_values["chrom"]
    start = input_values["start"]
    end = input_values["end"]
    family = input_values["family"]
    strand = input_values["strand"]
    fasta1 = input_values["fasta1"]
    fasta2 = input_values["fasta2"]
    out_dir = input_values["out_dir"]
    flank_len = input_values["flank_len"]
    window = input_values["window"]
    bed2 = input_values["bed2"]
    preset = input_values["preset"]

    prefix = "_".join([chrom, str(start), str(end), family])
    lift_entries["ID"] = prefix
    ref_coord = chrom + ":" + str(start) + "-" + str(end)
    lift_entries["ref_coord"] = ref_coord

    te_length = int(end) - int(start)
    lift_entries["te_length"] = te_length

    # extract 5 prime flank
    start_5p_flank = int(start) - flank_len + 1
    end_5p_flank = int(start)
    prefix_5p = prefix + "_5p"
    fa_5p = get_ref_fa(fasta1, chrom, start_5p_flank, end_5p_flank, prefix_5p, out_dir)
    # extract 3 prime flank
    start_3p_flank = int(end) + 1
    end_3p_flank = int(end) + 1 + flank_len
    prefix_3p = prefix + "_3p"
    fa_3p = get_ref_fa(fasta1, chrom, start_3p_flank, end_3p_flank, prefix_3p, out_dir)

    # intiate minimap2 parameters
    num_secondary = 10  # save 10 secondary alignments

    # map 5 prime flanks to ref2
    align_5p = out_dir + "/" + prefix_5p + ".paf"
    align_flank(fa_5p, fasta2, align_5p, preset, num_secondary)
    align_5p_bed = out_dir + "/" + prefix_5p + ".bed"
    paf_to_bed(
        align_5p, align_5p_bed, filter=chrom
    )  # only keep alignments on the same chromosome as original TE

    # map 3 prime flanks to ref2
    align_3p = out_dir + "/" + prefix_3p + ".paf"
    align_flank(fa_3p, fasta2, align_3p, preset, num_secondary)
    align_3p_bed = out_dir + "/" + prefix_3p + ".bed"
    paf_to_bed(
        align_3p, align_3p_bed, filter=chrom
    )  # only keep alignments on the same chromosome as original TE

    # clean intermediate files
    os.remove(align_3p)
    os.remove(align_5p)
    os.remove(fa_5p)
    os.remove(fa_3p)

    # find closest entries between 3 prime and 5 prime flank alignments, require two alignments to be on the same strand, also report distance
    overlap = out_dir + "/" + prefix + ".overlap"
    if check_file(align_5p_bed) and check_file(align_3p_bed):
        with open(overlap, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "closest",
                    "-a",
                    align_5p_bed,
                    "-b",
                    align_3p_bed,
                    "-s",
                    "-d",
                    "-t",
                    "all",
                ],
                stdout=output,
            )

    # identify liftover coordinates and report in json format
    lift_entries["report"] = []
    lift_start = 0
    lift_end = 0
    num_hits = 0  # number of non-reference liftover hits
    reported = False
    if check_file(overlap):
        with open(overlap, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                chrom_ref2_5p = entry[0]
                chrom_ref2_3p = entry[6]
                flank_strand = entry[5]
                mapp_quality_5p = int(entry[4])
                mapp_quality_3p = int(entry[10])
                if (
                    chrom_ref2_5p == chrom_ref2_3p == chrom
                ):  # to make sure the entry exists
                    start_5p = int(entry[1])
                    end_5p = int(entry[2])
                    start_3p = int(entry[7])
                    end_3p = int(entry[8])
                    # get final liftover coordiantes
                    lift_start, lift_end, lift_gap = get_coord(
                        start_5p, end_5p, start_3p, end_3p, flank_strand
                    )
                    # figure out the strand of the lifted annotation
                    if flank_strand == strand:
                        lift_strand = "+"
                    else:
                        lift_strand = "-"
                    # report the flanking sequence alignments
                    align_5p = chrom_ref2_5p + ":" + str(start_5p) + "-" + str(end_5p)
                    align_3p = chrom_ref2_3p + ":" + str(start_3p) + "-" + str(end_3p)
                    lift_entry = {
                        "type": None,
                        "family": family,
                        "chrom": chrom,
                        "start": lift_start,
                        "end": lift_end,
                        "mapp_quality_5p": mapp_quality_5p,
                        "mapp_quality_3p": mapp_quality_3p,
                        "strand": lift_strand,
                        "gap": lift_gap,
                        "align_5p": align_5p,
                        "align_3p": align_3p,
                        "distance_ref_5p": None,
                        "distance_ref_3p": None,
                        "comment": None,
                    }
                    # check the distance between flank alignment and nearast TE annotation
                    distance_5p = check_nearby_ref(
                        chrom,
                        start_5p,
                        end_5p,
                        family,
                        lift_strand,
                        bed2,
                        out_dir,
                    )
                    distance_3p = check_nearby_ref(
                        chrom,
                        start_3p,
                        end_3p,
                        family,
                        lift_strand,
                        bed2,
                        out_dir,
                    )
                    if distance_5p is not None:
                        lift_entry["distance_ref_5p"] = str(distance_5p)
                    if distance_3p is not None:
                        lift_entry["distance_ref_3p"] = str(distance_3p)

                    # if the overlap between flank alignments is greater than 50bp, don't report
                    if lift_gap < -window and lift_gap > -500:
                        if not reported:
                            reported = False  # redundant
                    elif abs(lift_gap) <= window:
                        # if the gap between flank alignments is smaller than 50bp
                        # report as reference if 1) there is a reference 2 TE in between two flanks or 2) the gap size and original TE size is similar or 3) the gap size is bigger than the size of the TE.
                        # report as non-reference otherwise.
                        if (
                            (
                                distance_5p is not None
                                and distance_5p >= 0
                                and distance_5p <= lift_gap
                                and distance_3p is not None
                                and distance_3p <= 0
                                and -distance_3p <= lift_gap
                            )
                            or check_nums_similar(lift_gap, te_length)
                            or (lift_gap >= te_length)
                        ):
                            lift_entry["type"] = "reference"
                            lift_entry[
                                "comment"
                            ] = "flanks overlap/gap size within threshold, include ref2 TE in between"
                        else:
                            lift_entry["type"] = "non-reference"
                            lift_entry[
                                "comment"
                            ] = "flanks overlap/gap size within threshold"
                            num_hits = num_hits + 1
                        lift_entries["report"].append(lift_entry)
                        reported = True
                    else:
                        # if the gap between flanks is greater than 50bp
                        if lift_gap <= 0.5 * te_length:
                            # if the gap size is smaller than half of the original TE size
                            # report as reference if 1) there is a reference 2 TE in between two flanks
                            # report as non-reference otherwise.
                            if (
                                distance_5p is not None
                                and distance_5p >= 0
                                and distance_5p <= lift_gap
                                and distance_3p is not None
                                and distance_3p <= 0
                                and -distance_3p <= lift_gap
                            ):
                                lift_entry["type"] = "reference"
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size less than half of TE annotation, include ref2 TE in between"
                            else:
                                lift_entry["type"] = "non-reference"
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size exceeds threshold but less than half of TE annotation, no ref2 TE in between"
                                num_hits = num_hits + 1
                            lift_entries["report"].append(lift_entry)
                            reported = True
                        elif lift_gap >= 0.5 * te_length and lift_gap <= 20000:
                            # if the gap size is greater than half of the original TE size and smaller than 20kb
                            lift_entry["type"] = "reference"
                            if (
                                distance_5p is not None
                                and distance_5p >= 0
                                and distance_5p <= lift_gap
                                and distance_3p is not None
                                and distance_3p <= 0
                                and -distance_3p <= lift_gap
                            ):
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size greater than half of TE annotation, include ref2 TE in between"
                            else:
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size greater than half of TE annotation, no ref2 TE in between"
                            lift_entries["report"].append(lift_entry)
                            reported = True
                        else:
                            # if the gap size is greater than 20kb
                            if not reported:
                                reported = False  # redundant

    # more than one report per annotation, need to filter
    # first step is to choose the best reference annotation by comparing gap size with TE length
    if len(lift_entries["report"]) > 1:
        best_ref_entry = dict()
        best_nonref_entry = dict()
        for report in lift_entries["report"]:
            if report["type"] == "reference":
                if not best_ref_entry:
                    best_ref_entry = report
                else:
                    current_gap_size = best_ref_entry["gap"]
                    new_gap_size = report["gap"]
                    if choose_new_size(te_length, current_gap_size, new_gap_size):
                        best_ref_entry = report
            if report["type"] == "non-reference":
                if not best_nonref_entry:
                    best_nonref_entry = report
                else:
                    reported = False
        lift_entries["report"] = []
        if reported:
            if best_ref_entry and best_nonref_entry:
                # if both reference and non-reference liftover can be found, report the non-reference one
                # lift_entries["report"].append(best_ref_entry)
                lift_entries["report"].append(best_nonref_entry)
            elif best_ref_entry:
                lift_entries["report"].append(best_ref_entry)
            elif best_nonref_entry:
                lift_entries["report"].append(best_nonref_entry)
            else:
                lift_entries["report"] = []
                reported = False

    if not reported:
        # TODO: if only one flank can be lifted, check to see if there is a reference TE nearby (same family, same strand, similar size?)
        align_5p = []
        align_3p = []
        if check_file(align_5p_bed):
            with open(align_5p_bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    align_5p.append(entry[0] + ":" + entry[1] + "-" + entry[2])
        if check_file(align_3p_bed):
            with open(align_3p_bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    align_3p.append(entry[0] + ":" + entry[1] + "-" + entry[2])
        if len(align_5p) == 0:
            align_5p = None
        if len(align_3p) == 0:
            align_3p = None
        lift_entry = {
            "type": "unlifted",
            "family": family,
            "chrom": None,
            "start": None,
            "end": None,
            "mapp_quality_5p": None,
            "mapp_quality_3p": None,
            "strand": None,
            "gap": None,
            "align_5p": align_5p,
            "align_3p": align_3p,
            "comment": "flank alignments not nearby each other / only one flank aligned",
        }
        lift_entries["report"].append(lift_entry)
    lift_entries["num_hits"] = num_hits

    # write
    liftover_report = out_dir + "/" + prefix + "_report.json"
    with open(liftover_report, "w") as output:
        json.dump(lift_entries, output, indent=4, sort_keys=False)

    return liftover_report


def choose_new_size(size_ref, size_old, size_new):
    if size_ref - size_old > size_ref - size_new:
        return True
    else:
        return False


def check_nums_similar(num1, num2):
    normalize_diff = abs(num1 - num2) / num2
    if normalize_diff <= 0.1:
        return True
    else:
        return False


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        print("Creation of the directory %s failed" % dir)


def check_file(filename):
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return True
    else:
        return False


def build_index(fa):
    subprocess.call(["samtools", "faidx", fa])


def main():
    args = get_args(program_version=__version__)

    # TODO: merge overlapped/nested insertions?

    # generate ref1 index file if not present
    ref1_index = args.fasta1 + ".fai"
    if not check_file(ref1_index):
        build_index(args.fasta1)

    # loop through TE annotations, prepare data for parallel liftover jobs
    input_json_dir = os.path.join(args.out, "input_json")
    mkdir(input_json_dir)
    tmp_dir = os.path.join(args.out, "tmp")
    mkdir(tmp_dir)

    liftover_pa_list = []
    with open(args.bed1, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[0]
            start = entry[1]
            end = entry[2]
            family = entry[3]
            strand = entry[5]

            prefix = "_".join([chrom, str(start), str(end), family])
            liftover_entry = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "family": family,
                "strand": strand,
                "fasta1": args.fasta1,
                "fasta2": args.fasta2,
                "out_dir": tmp_dir,
                "flank_len": args.flank_len,
                "window": args.window,
                "bed2": args.bed2,
                "preset": args.preset,
            }
            liftover_input = input_json_dir + "/" + prefix + "_input.json"
            with open(liftover_input, "w") as output:
                json.dump(liftover_entry, output, indent=4, sort_keys=False)
            liftover_pa_list.append(liftover_input)

    # run liftover jobs for all annotations in parallel
    print("Perform liftover...")
    try:
        pool = Pool(processes=args.threads)
        liftover_report_list = pool.map(run_liftover, liftover_pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("liftover failed, exiting...")
        sys.exit(1)

    # final report
    print("Generate final report...")
    data = []
    for liftover_report in liftover_report_list:
        with open(liftover_report) as f:
            report = json.load(f)
            data.append(report)

    json_report = args.out + "/" + "liftover_report.json"
    with open(json_report, "w") as output:
        json.dump(data, output, indent=4, sort_keys=False)

    # add a step to merge close annotations if inserions are : 1) supported by only one flank 2) same family 3) same strand

    # write non-reference lifted annotations in BED format
    bed_report = args.out + "/" + "liftover_nonref.bed"
    with open(bed_report, "w") as output:
        for item in data:
            if item["num_hits"] == 1:
                info = item["report"][0]
                chrom = info["chrom"]
                start = info["start"]
                end = info["end"]
                family = info["family"]
                score = "."
                strand = info["strand"]
                out_line = "\t".join(
                    [chrom, str(start), str(end), family, score, strand]
                )
                output.write(out_line + "\n")

    # generate summary
    nonref_total = 0
    ref_total = 0
    unlift_total = 0
    nonref_comments = dict()
    ref_comments = dict()
    unlift_comments = dict()
    for item in data:
        infos = item["report"]
        for info in infos:
            if info["type"] == "non-reference":
                nonref_total = nonref_total + 1
                if "comment" in info.keys():
                    if info["comment"] not in nonref_comments.keys():
                        nonref_comments[info["comment"]] = 1
                    else:
                        nonref_comments[info["comment"]] = (
                            nonref_comments[info["comment"]] + 1
                        )
            if info["type"] == "reference":
                ref_total = ref_total + 1
                if "comment" in info.keys():
                    if info["comment"] not in ref_comments.keys():
                        ref_comments[info["comment"]] = 1
                    else:
                        ref_comments[info["comment"]] = (
                            ref_comments[info["comment"]] + 1
                        )
            if info["type"] == "unlifted":
                unlift_total = unlift_total + 1
                if "comment" in info.keys():
                    if info["comment"] not in unlift_comments.keys():
                        unlift_comments[info["comment"]] = 1
                    else:
                        unlift_comments[info["comment"]] = (
                            unlift_comments[info["comment"]] + 1
                        )
    summary_data = {
        "non-reference": {"total": nonref_total, "comments": nonref_comments},
        "reference": {"total": ref_total, "comments": ref_comments},
        "unlifted": {"total": unlift_total, "comments": unlift_comments},
    }

    summary_report = args.out + "/" + "liftover_summary.json"
    with open(summary_report, "w") as output:
        json.dump(summary_data, output, indent=4, sort_keys=False)

    # clean tmp files
    if not args.keep_files:
        print("Clean intermediate files...")
        shutil.rmtree(tmp_dir)
        shutil.rmtree(input_json_dir)

    print("Liftover finished!")


main()

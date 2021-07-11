#!/usr/bin/env python3
import sys
import argparse
import os
import subprocess
import glob
import json
from eval_utility import filter_annotation, count_lines


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to evaluate TELR predictions in terms of coordinate, family and AF"
    )

    ## required ##
    parser.add_argument(
        "-i",
        "--telr_out_dir",
        type=str,
        help="directory that include TELR runs",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        help="Regions in reference genome assembly to include in the analysis in BED format",
        required=False,
    )
    parser.add_argument(
        "--include_families",
        type=str,
        help="TE families to use in the analysis (separated by comma)",
        required=False,
    )
    parser.add_argument(
        "--exclude_families",
        type=str,
        help="TE families to exclude in the analysis (separated by comma)",
        required=False,
    )
    parser.add_argument(
        "--exclude_nested",
        action="store_true",
        help="Exclude nested TEs in the analysis",
        required=False,
    )
    parser.add_argument(
        "-b",
        "--bed",
        type=str,
        help="TE annotation for reference genome assembly in BED format",
        required=True,
    ),
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        help="window size when comparing TELR prediction with TE annotation (default = '5bp')",
        required=False,
    )
    parser.add_argument("-o", "--out", type=str, help="output directory", required=True)

    args = parser.parse_args()

    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.window is None:
        args.window = 5

    return args


def main():
    args = get_args()

    # filter TE annotation file
    annotation_filtered = os.path.join(args.out, "annotation.filter.bed")
    filter_annotation(
        bed_in=args.bed,
        bed_out=annotation_filtered,
        filter_region=args.region,
        include_families=args.include_families,
        exclude_families=args.exclude_families,
    )

    # find and read BED files from TELR runs
    pattern = "/**/*telr.bed"  # TODO: update here
    pred_file_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    print(pred_file_list)

    # for each file, do evaluations
    summary_dict = dict()
    for pred_file in pred_file_list:
        prefix = os.path.basename(pred_file).replace(".telr.bed", "")
        print("prefix: " + prefix)

        # filter TELR predictions by region and family
        pred_filtered = args.out + "/" + prefix + ".parse.filter.bed"
        filter_annotation(
            bed_in=pred_file,
            bed_out=pred_filtered,
            filter_region=args.region,
            include_families=args.include_families,
            exclude_families=args.exclude_families,
        )

        # get the number of predictions
        num_pred_total = count_lines(pred_filtered)
        print("Total TELR predictions:" + str(num_pred_total))

        # compare with lift over set
        overlap = args.out + "/" + prefix + ".overlap.bed"
        with open(overlap, "w") as output:  # TODO: strand?
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(args.window),
                    "-a",
                    pred_filtered,
                    "-b",
                    annotation_filtered,
                ],
                stdout=output,
            )

        # parse overlap file and get the number of matched insertions
        match_set = set()
        with open(overlap, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                family1 = entry[3]  # telr
                family2 = entry[9]  # liftover
                if (
                    family1 == family2
                ):  # TODO: think about nested insertions, should I remove them in both set?
                    te_locus = "_".join([entry[0], entry[1], entry[2], entry[3]])
                    match_set.add(te_locus)
                else:
                    print(line)
        num_tp = len(match_set)  # TEs predicted by TELR that are true positives
        print("Number of TELR true positives:" + str(num_tp))

        num_fp = num_pred_total - num_tp
        print(
            "Number of TELR False positives:" + str(num_fp)
        )  # TEs predicted by TELR that are false positives

        # get telr only set # TODO: do I need this?
        telr_only = args.out + "/" + prefix + ".telr_only.bed"
        with open(telr_only, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(args.window),
                    "-a",
                    pred_filtered,
                    "-b",
                    annotation_filtered,
                    "-v",
                ],
                stdout=output,
            )
        num_fp2 = count_lines(
            telr_only
        )  # TEs predicted by TELR that are not in liftover set
        # print("num_fp2:" + str(num_fp2))

        # compare with lift over set, get things that aren't predicted
        lift_only = args.out + "/" + prefix + ".lift_only.bed"
        with open(lift_only, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(args.window),
                    "-a",
                    annotation_filtered,
                    "-b",
                    pred_filtered,
                    "-v",
                ],
                stdout=output,
            )
        num_fn = count_lines(lift_only)
        print("Number of false negatives:" + str(num_fn))

        precision = round(num_tp / num_pred_total, 3)
        num_liftover = count_lines(annotation_filtered)
        recall = round(num_tp / num_liftover, 3)

        coverage = prefix.split("_")[0]
        coverage = int(coverage.replace("x", ""))
        ploidy = prefix.split("_")[1]
        data_zygosity = prefix.split("_")[2]

        # get af eval
        n_hom = 0
        n_het = 0
        n_unknown = 0

        # document this
        if ploidy == "diploid":
            hom_min = 0.75
            het_min = 0.25
            het_max = 0.75
        elif ploidy == "tetraploid":
            hom_min = 0.9
            het_min = 0.1
            het_max = 0.4
        else:
            print("unrecognized ploidy, please provide diploid/tetraploid, exiting...")
            sys.exit(1)

        with open(pred_filtered, "r") as input:  # TODO: mention this in the doc
            for line in input:
                entry = line.replace("\n", "").split("\t")
                af = float(entry[4])
                if af > hom_min:
                    n_hom = n_hom + 1
                elif af > het_min and af < het_max:
                    n_het = n_het + 1
                else:
                    n_unknown = n_unknown + 1

        if "het" in data_zygosity:
            false_af_rate = round((n_hom + n_unknown) / num_pred_total, 3)
        else:
            false_af_rate = round((n_het + n_unknown) / num_pred_total, 3)

        summary_dict = {
            "prefix": prefix,
            "coverage": coverage,
            "ploidy": ploidy,
            "data_zygosity": data_zygosity,
            "num_pred_total": num_pred_total,
            "num_tp": num_tp,
            "num_fp": num_fp,
            # "num_fp2": num_fp2,
            "num_fn": num_fn,
            "precision": precision,
            "recall": recall,
            "num_hom": n_hom,
            "num_het": n_het,
            "num_unclassified": n_unknown,
            "false_af_rate": false_af_rate,
        }

    # write
    eval_out = args.out + "/eval_liftover.json"
    with open(eval_out, "w") as output:
        json.dump(summary_dict, output, indent=4, sort_keys=False)


if __name__ == "__main__":
    main()

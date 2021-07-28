#!/usr/bin/env python3
import argparse
import os
import glob
import json
from eval_utility import filter_annotation, count_lines


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to evaluate TELR TE allele frequency predictions"
    )

    ## required ##
    parser.add_argument(
        "-i",
        "--telr_out_dir",
        type=str,
        help="TELR output directory",
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
        "--ploidy",
        type=str,
        help="expected ploidy, please choose '2' or '4' (default = '2')",
        required=False,
    )
    parser.add_argument(
        "--genotype",
        type=str,
        help="expected genotype, please choose 'homozygous' or 'heterozygous' (default = 'heterozygous')",
        required=False,
    )
    parser.add_argument("-o", "--out", type=str, help="output directory", required=True)

    args = parser.parse_args()

    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.ploidy is None:
        args.ploidy = "2"
    elif args.ploidy not in ["2", "4"]:
        raise Exception("ploidy must be '2' or '4'")

    if args.genotype is None:
        args.genotype = "heterozygous"
    elif args.genotype not in ["homozygous", "heterozygous"]:
        raise Exception("genotype must be 'homozygous' or 'heterozygous'")

    return args


def main():
    args = get_args()

    # find and read BED and JSON files from TELR runs
    pattern = "/**/*telr.bed"  # TODO: update here
    pred_bed_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    pred_bed = pred_bed_list[0]
    print("TELR output BED file: " + pred_bed)

    pattern = "/**/*telr.json"  # TODO: update here
    pred_json_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    pred_json = pred_json_list[0]
    print("TELR output JSON file: " + pred_json)

    # for each file, do evaluations
    summary_dict = dict()
    prefix = os.path.basename(pred_bed).replace(".telr.bed", "")
    print("prefix: " + prefix)

    # filter TELR predictions by region and family
    pred_filtered = args.out + "/" + prefix + ".parse.filter.bed"
    filter_annotation(
        bed_in=pred_bed,
        bed_out=pred_filtered,
        filter_region=args.region,
        include_families=args.include_families,
        exclude_families=args.exclude_families,
    )

    # get the number of predictions
    num_pred_total = count_lines(pred_filtered)
    print("Total TELR predictions:" + str(num_pred_total))

    coverage = prefix.split("_")[0]
    coverage = int(coverage.replace("x", ""))
    ploidy = args.ploidy
    data_zygosity = args.genotype

    # get af eval
    n_hom = 0
    n_het = 0
    n_unknown = 0

    # document this
    if ploidy == "2":
        hom_min = 0.9
        het_min = 0.375
        het_max = 0.625
    else:
        hom_min = 0.9
        het_min = 0.1
        het_max = 0.4

    # read data into set
    te_ids = set()
    with open(pred_filtered) as f:
        for line in f:
            line = line.strip().split("\t")
            te_id = "_".join(line[0:4])
            te_ids.add(te_id)

    with open(pred_json, "r") as input:  # TODO: mention this in the doc
        telr_json = json.load(input)
        for item in telr_json:
            te_id = item["ID"]
            if te_id in te_ids:
                af = item["allele_frequency"]
                if af is None:
                    n_unknown += 1
                elif af >= het_min and af <= het_max:
                    n_het += 1
                elif af >= hom_min:
                    n_hom += 1
                else:
                    n_unknown += 1

    if data_zygosity == "heterozygous":
        false_af_rate = round((n_hom + n_unknown) / num_pred_total, 3)
    else:
        false_af_rate = round((n_het + n_unknown) / num_pred_total, 3)

    summary_dict = {
        "prefix": prefix,
        "coverage": coverage,
        "ploidy": ploidy,
        "data_zygosity": data_zygosity,
        "num_pred_total": num_pred_total,
        "num_hom": n_hom,
        "num_het": n_het,
        "num_unclassified": n_unknown,
        "false_af_rate": false_af_rate,
    }

    # write
    eval_out = args.out + "/telr_eval_af.json"
    with open(eval_out, "w") as output:
        json.dump(summary_dict, output, indent=4, sort_keys=False)


if __name__ == "__main__":
    main()

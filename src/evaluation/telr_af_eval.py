#!/usr/bin/env python3
import argparse
import os
import glob
import json
from telr_eval_utility import filter_annotation, count_lines


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
        "--exclude_nested",
        action="store_true",
        help="Exclude nested TEs in the analysis",
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
        args.genotype = "homozygous"
    elif args.genotype not in [
        "homozygous",
        "heterozygous",
        "simplex",
        "duplex",
        "triplex",
        "quadruplex",
    ]:
        raise Exception(
            "genotype must be 'homozygous/heterozygous/simplex/duplex/triplex/quadruplex'"
        )

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
        exclude_nested=args.exclude_nested,
    )

    coverage = prefix.split("_")[0]
    coverage = int(coverage.replace("x", ""))
    ploidy = args.ploidy
    data_zygosity = args.genotype

    # TAF cutoff
    ## under tetraploid
    simplex_min = 0
    simplex_max = 0.375
    duplex_min = 0.375
    duplex_max = 0.625
    triplex_min = 0.625
    triplex_max = 0.875
    quadruplex_min = 0.875
    quadruplex_max = 1
    ## under diploid
    homozygous_min = 0.875
    homozygous_max = 1
    heterozygous_min = 0.375
    heterozygous_max = 0.625

    n_simplex = 0
    n_duplex = 0
    n_triplex = 0
    n_quadruplex = 0
    n_hom = 0
    n_het = 0
    n_unknown = 0
    n_total = 0

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
                if af is not None:
                    n_total += 1
                    if ploidy == "2":
                        if af > homozygous_min and af <= homozygous_max:
                            n_hom += 1
                        elif af > heterozygous_min and af <= heterozygous_max:
                            n_het += 1
                        else:
                            n_unknown += 1
                    if ploidy == "4":
                        if af > simplex_min and af <= simplex_max:
                            n_simplex += 1
                        elif af > duplex_min and af <= duplex_max:
                            n_duplex += 1
                        elif af > triplex_min and af <= triplex_max:
                            n_triplex += 1
                        elif af > quadruplex_min and af <= quadruplex_max:
                            n_quadruplex += 1

    if ploidy == "2":
        if data_zygosity == "homozygous":
            precison = n_hom / n_total
        elif data_zygosity == "heterozygous":
            precison = n_het / n_total
        else:
            raise Exception("data zygosity must be 'homozygous' or 'heterozygous'")
            sys.exit(1)
    else:
        if data_zygosity == "simplex":
            precison = n_simplex / n_total
        elif data_zygosity == "duplex":
            precison = n_duplex / n_total
        elif data_zygosity == "triplex":
            precison = n_triplex / n_total
        elif data_zygosity == "quadruplex":
            precison = n_quadruplex / n_total
        else:
            raise Exception(
                "data zygosity must be 'simplex', 'duplex', 'triplex', or 'quadruplex'"
            )
            sys.exit(1)

    if ploidy == "2":
        summary_dict = {
            "prefix": prefix,
            "coverage": coverage,
            "ploidy": ploidy,
            "data_zygosity": data_zygosity,
            "num_pred_total": n_total,
            "num_pred_hom": n_hom,
            "num_pred_het": n_het,
            "num_pred_unclassified": n_unknown,
            "precison": precison,
        }
    else:
        summary_dict = {
            "prefix": prefix,
            "coverage": coverage,
            "ploidy": ploidy,
            "data_zygosity": data_zygosity,
            "num_pred_total": n_total,
            "num_simplex": n_simplex,
            "num_duplex": n_duplex,
            "num_triplex": n_triplex,
            "num_quadruplex": n_quadruplex,
            "precison": precison,
        }

    # write
    eval_out = args.out + "/telr_eval_af.json"
    with open(eval_out, "w") as output:
        json.dump(summary_dict, output, indent=4, sort_keys=False)


if __name__ == "__main__":
    main()

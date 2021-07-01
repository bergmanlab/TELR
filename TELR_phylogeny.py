#!/usr/bin/env python3
import argparse
import sys
import os
import glob
import json
import random
import string
import editdistance
from multiprocessing import Pool
from TELR_te import create_fa

# import time
import logging
import subprocess
from Bio import SeqIO
import shutil
from TELR_utility import mkdir

"""
This script builds maximum likelihood tree using TE sequences from TELR run
author: Shunhua Han
"""


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to detect build phylogeny from TE sequences"
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    # required
    optional.add_argument(
        "--family",
        type=str,
        help="TE families (separated by comma)",
        required=True,
    )

    required.add_argument(
        "--telr_dirs",
        type=str,
        help="list of TELR output directories",
        nargs="+",
        required=True,
    )

    required.add_argument(
        "--consensus",
        type=str,
        help="TE consensus sequence",
        required=True,
    )

    # optional
    optional.add_argument(
        "--out",
        type=str,
        help="directory to output data (default = '.')",
        required=False,
    )
    optional.add_argument(
        "--thread",
        type=int,
        help="max cpu threads to use (default = '1')",
        required=False,
    )
    optional.add_argument(
        "--bootstrap",
        type=int,
        help="bootstrap number (only apply when raxml is used for creating phylogeny)",
        required=False,
    )
    optional.add_argument(
        "--method",
        type=str,
        help="method to create phylogeny, raxml/iqtree/both (default: raxml)",
        required=False,
    )
    optional.add_argument(
        "--add_consensus",  # TODO
        action="store_true",
        help="If provided then add consensus sequence to the phylogeny (default: don't add consensus)",
        required=False,
    )
    optional.add_argument(
        "--allow_nested",  # TODO
        action="store_true",
        help="If provided then allow nested/composite sequences in the phylogeny (default: don't allow)",
        required=False,
    )
    optional.add_argument(
        "--exclude_contigs",
        type=str,
        help="Exclude contigs can be provided, multiple contigs should be separated by comma",
        required=False,
    )
    optional.add_argument(
        "--length_filter",
        type=float,
        help="percentage of TE sequence longer or shorter than consensus sequence (default: 100%%)",
        required=False,
    )
    optional.add_argument(
        "--divergence_filter",
        type=float,
        help="percentage of TE sequence divergent from consensus sequence (default: 10%%)",
        required=False,
    )
    parser.add_argument(
        "-t", "--threads", type=int, help="number of cores", required=False
    )
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # checks if in files exist
    try:
        test = open(args.consensus, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.consensus)
        sys.exit(1)

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    mkdir(args.out)

    # set up default value for optional arguments
    if args.thread is None:
        args.thread = 1

    if args.method is None:
        args.method = "iqtree"
    elif args.method != "raxml" and args.method != "iqtree" and args.method != "both":
        print("method not recognized, please check help page")
        sys.exit(1)

    if args.length_filter is None:
        args.length_filter = 1

    if args.divergence_filter is None:
        args.divergence_filter = 0.1

    if args.bootstrap is None:
        args.bootstrap = 5

    return args


def get_te_seq(te_seqs, family, consensus, add_consensus, out):
    with open(out, "w") as output:
        if add_consensus:
            record = consensus[family]
            output.write(">" + record.id + "\n" + str(record.seq) + "\n")
        for item in te_seqs:
            item_family = item["family"]
            item_quality = item["quality"]
            if item_family == family and item_quality == "pass":
                item_id = item["ID"]
                item_seq = item["sequence"]
                output.write(">" + item_id + "\n" + item_seq + "\n")


def run_raxml(alignment, tmp_dir, out_dir, prefix, bootstrap, thread):
    commands = (
        "raxmlHPC-PTHREADS-AVX2 --silent -f a -x 10350 -p 12345 -m GTRGAMMA --no-bfgs"
        + " -s "
        + alignment
        + " -n "
        + prefix
        + " -T "
        + str(thread)
        + " -w "
        + tmp_dir
        + " -# "
        + str(bootstrap)
    )
    subprocess.call(commands, shell=True)
    tree_file = os.path.join(tmp_dir, "RAxML_bipartitions." + prefix)
    if os.path.isfile(tree_file):
        new_tree_file = os.path.join(out_dir, prefix + ".raxml.tree")
        os.rename(tree_file, new_tree_file)


def run_iqtree(alignment, tmp_dir, out_dir, prefix):
    commands = (
        "iqtree -m GTR+G -T AUTO -bb 1000 --redo"
        + " -s "
        + alignment
        + " --prefix "
        + os.path.join(tmp_dir, prefix)
    )
    subprocess.call(commands, shell=True)
    tree_file = os.path.join(tmp_dir, prefix + ".treefile")
    if os.path.isfile(tree_file):
        new_tree_file = os.path.join(out_dir, prefix + ".iqtree.tree")
        os.rename(tree_file, new_tree_file)
    else:
        print("tree failed")
        sys.exit(1)


def build_tree(alignment, family, out_dir, tmp_dir, method, bootstrap, thread):
    if method == "raxml" or method == "both":
        run_raxml(
            alignment=alignment,
            prefix=family,
            tmp_dir=tmp_dir,
            out_dir=out_dir,
            bootstrap=bootstrap,
            thread=thread,
        )
    if method == "iqtree" or method == "both":
        run_iqtree(alignment=alignment, prefix=family, tmp_dir=tmp_dir, out_dir=out_dir)


def get_telr_seqs_phylogeny(
    family,
    consensus,
    data_json,
    out_dir,
    method,
    bootstrap,
    thread,
    add_consensus=False,
):
    # read parsed data
    with open(data_json) as f:
        data = json.load(f)
    # create dir for intermediate files
    tree_tmp_dir = os.path.join(out_dir, "tree_tmp_files", family)
    if os.path.isdir(tree_tmp_dir):
        shutil.rmtree(tree_tmp_dir)
    os.makedirs(tree_tmp_dir)

    # build consensus sequence from TELR TE output
    telr_seqs = os.path.join(tree_tmp_dir, family + ".telr.fa")
    get_te_seq(data, family, consensus, add_consensus, telr_seqs)

    # align TE sequences
    align_fa = os.path.join(tree_tmp_dir, family + ".align.fa")
    with open(align_fa, "w") as output:
        if thread > 8:
            mafft_thread = 8
        else:
            mafft_thread = thread
        subprocess.call(
            ["mafft", "--quiet", "--auto", "--thread", str(mafft_thread), telr_seqs],
            stdout=output,
        )

    # build tree
    build_tree(
        alignment=align_fa,
        tmp_dir=tree_tmp_dir,
        out_dir=out_dir,
        family=family,
        thread=thread,
        method=method,
        bootstrap=bootstrap,
    )


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return "".join(random.choice(chars) for _ in range(size))


def per_base_len_dist(consensus_len, seq_len):
    return (seq_len - consensus_len) / consensus_len


def diff_telr_consensus(consensus_seq, telr_seq):
    # create fasta file
    consensus_id = id_generator(3, string.ascii_uppercase)
    os.path.join(out_dir, consensus_id + ".fa")
    create_fa(consensus_id, consensus_seq, out_dir)
    consensus_length = len(consensus_seq)
    telr_len = len(telr_seq)
    len_diff = per_base_len_dist(consensus_length, telr_len)
    div_diff = editdistance.eval(consensus_seq, telr_seq) / consensus_length
    return len_diff, div_diff


def parse_telr_out(args):
    json_in = args[0]
    json_out = args[1]
    consensus = args[2]
    len_filter = args[3]
    div_filter = args[4]
    exclude_contigs = args[5]
    contigs_fail = set()
    if exclude_contigs:
        contigs = exclude_contigs.split(",")
        for contig in contigs:
            contigs_fail.add(contig)
    with open(json_in) as f:
        item = json.load(f)
    sample = item["sample"]
    family = item["family"]
    telr_seq = item["sequence"]
    len_diff, div_diff = diff_telr_consensus(consensus, telr_seq)
    contig = item["chr"]
    start = item["start"]
    new_te_id = "_".join([sample, contig, str(start), family])
    item["ID"] = new_te_id
    item["length_diff"] = len_diff
    item["edit_dist"] = div_diff
    if (
        abs(len_diff) <= len_filter
        and div_diff <= div_filter
        and contig not in contigs_fail
    ):
        item["quality"] = "pass"
    else:
        item["quality"] = "fail"
    with open(json_out, "w") as c:
        json.dump(item, c, indent=4)
    os.remove(json_in)


def main():
    args = get_args()

    # TODO process consensus sequence, remove stuff after #
    consensus_dict = SeqIO.index(args.consensus, "fasta")

    families = args.family.replace(" ", "").replace("_", "-").split(",")
    # parse TELR output directories and load consensus
    json_dir = os.path.join(args.out, "json_files")
    mkdir(json_dir)
    print("parse TELR output")
    meta_data = []
    for telr_out_dir in args.telr_dirs:
        if telr_out_dir[-1] == "/":
            telr_out_dir = telr_out_dir[:-1]
        telr_json_files = glob.glob(telr_out_dir + "/**/*telr.json", recursive=True)
    for telr_json in telr_json_files:
        sample = os.path.basename(telr_json).replace(".telr.json", "").replace("_", "-")
        with open(telr_json) as f:
            json_data = json.load(f)
        for item in json_data:
            family = item["family"]
            if family in families:
                item["sample"] = sample
                prefix = "_".join([sample, item["ID"], family])
                item_json_in = os.path.join(json_dir, prefix + ".in.json")
                with open(item_json_in, "w") as c:
                    json.dump(item, c, indent=4)
                item_json_out = os.path.join(json_dir, prefix + ".out.json")
                consensus_seq = str(consensus_dict[family].seq)
                item_out = {
                    "json_in": item_json_in,
                    "json_out": item_json_out,
                    "consensus": consensus_seq,
                    "len_filter": args.length_filter,
                    "div_filter": args.divergence_filter,
                    "exclude_contigs": args.exclude_contigs,
                }
                meta_data.append(item_out)

    # run parsing
    # read table
    pa_list = []
    for meta_item in meta_data:
        pa_list.append(
            [
                meta_item["json_in"],
                meta_item["json_out"],
                meta_item["consensus"],
                meta_item["len_filter"],
                meta_item["div_filter"],
                meta_item["exclude_contigs"],
            ]
        )

    # download in parallel
    print("TELR parsing...")
    try:
        pool = Pool(processes=args.threads)
        pool.map(parse_telr_out, pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("TELR parsing failed, exiting...")
        sys.exit(1)

    # merge parsed TELR output
    parsed_data = []
    for meta_item in meta_data:
        json_out = meta_item["json_out"]
        with open(json_out) as f:
            parsed_data_item = json.load(f)
        parsed_data.append(parsed_data_item)
        os.remove(json_out)

    # write seq quality report
    parsed_data_out = args.out + "/data_parsed.json"
    with open(parsed_data_out, "w") as c:
        json.dump(parsed_data, c, indent=4)

    for family in families:
        get_telr_seqs_phylogeny(
            data_json=parsed_data_out,
            family=family,
            consensus=consensus_dict,
            out_dir=args.out,
            method=args.method,
            bootstrap=args.bootstrap,
            thread=args.threads,
            add_consensus=args.add_consensus,
        )


main()
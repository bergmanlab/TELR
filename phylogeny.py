#!/usr/bin/env python3
import argparse
import sys
import os
# import time
import logging
import subprocess
from Bio import SeqIO
import shutil
# from TELR_utility import mkdir

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
    required.add_argument(
        "-s",
        "--seq",
        type=str,
        help="TE sequences in fasta format",
        required=True,
    )

    required.add_argument(
        "-c",
        "--consensus",
        type=str,
        help="TE consensus sequence",
        required=True,
    )

    # optional
    optional.add_argument(
        "-o",
        "--out",
        type=str,
        help="directory to output data (default = '.')",
        required=False,
    )
    optional.add_argument(
        "-t",
        "--thread",
        type=int,
        help="max cpu threads to use (default = '1')",
        required=False,
    )
    optional.add_argument(
        "-f",
        "--family",
        type=str,
        help="TE families (separated by comma, default = 'all)",
        required=False,
    )
    optional.add_argument(
        "-b",
        "--bootstrap",
        type=int,
        help="bootstrap number (only apply when raxml is used for creating phylogeny)",
        required=False,
    )
    optional.add_argument(
        "-m",
        "--method",
        type=str,
        help="method to create phylogeny, raxml/iqtree/both (default: raxml)",
        required=False,
    )
    optional.add_argument(
        "--length_filter",
        type=float,
        help="percentage of TE sequence longer or shorter than consensus sequence (default: 0.1)",
        required=False,
    )
    parser._action_groups.append(optional)
    args = parser.parse_args()
    
    # checks if in files exist
    try:
        test = open(args.seq, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.seq)
        sys.exit(1)

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
        args.method = "raxml"
    elif args.method != "raxml" and args.method != "iqtree" and args.method != "both":
        print("method not recognized, please check help page")
        sys.exit(1)

    if args.family is None:
        args.family = "all"

    if args.length_filter is None:
        args.length_filter = 0.1

    if args.bootstrap is None:
        args.bootstrap = 5

    return args


def get_te_seq(te_fa, family, consensus_length, length_filter, out):
    sample_name = os.path.basename(te_fa).replace(".final.fa", "")
    with open(te_fa, "r") as input, open(out, "w") as output:
        for record in SeqIO.parse(input, "fasta"):
            if record.id.split("#")[1] == family:
                if (
                    ((1 - length_filter) * float(consensus_length))
                    <= len(record.seq)
                    <= ((1 + length_filter) * float(consensus_length))
                ):
                    new_id = "_".join(
                        [sample_name, record.id.split("_")[0], record.id.split("_")[1]]
                    )
                    output.write(">" + new_id + "\n" + str(record.seq) + "\n")


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


def phylogeny_from_telr(
    family, consensus, length_filter, raw_fa, out_dir, method, bootstrap, thread
):
    # create dir for intermediate files
    tree_tmp_dir = os.path.join(out_dir, "tree_tmp_files", family)
    if os.path.isdir(tree_tmp_dir):
        shutil.rmtree(tree_tmp_dir)
    os.makedirs(tree_tmp_dir)

    # get consensus length/merge
    consensus_length = len(consensus[family].seq)

    # build consensus sequence from TELR TE output
    seqs = os.path.join(tree_tmp_dir, family + ".telr.fa")
    get_te_seq(raw_fa, family, consensus_length, length_filter, seqs)

    # align TE sequences
    align_fa = os.path.join(tree_tmp_dir, family + ".align.fa")
    with open(align_fa, "w") as output:
        subprocess.call(
            ["mafft", "--auto", "--quiet", "--thread", str(thread), seqs],
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


def main():
    args = get_args()
    consensus_dict = SeqIO.index(args.consensus, "fasta")
    families = args.family.replace(" ", "").split(",")
    for family in families:
        phylogeny_from_telr(
            family=family,
            consensus=consensus_dict,
            length_filter=args.length_filter,
            raw_fa=args.seq,
            out_dir=args.out,
            method=args.method,
            bootstrap=args.bootstrap,
            thread=args.thread,
        )


main()
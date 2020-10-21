#!/usr/bin/env python3
import argparse
import sys
import os
import glob
import editdistance

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
        "--length_filter",
        type=float,
        help="percentage of TE sequence longer or shorter than consensus sequence (default: 10%%)",
        required=False,
    )
    optional.add_argument(
        "--divergence_filter",
        type=float,
        help="percentage of TE sequence divergent from consensus sequence (default: 10%%)",
        required=False,
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
        for record in te_seqs:
            if record.id.split("_")[3] == family:
                output.write(">" + record.id + "\n" + str(record.seq) + "\n")


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
    family,
    consensus,
    te_seqs,
    out_dir,
    method,
    bootstrap,
    thread,
    add_consensus=False,
):
    # create dir for intermediate files
    tree_tmp_dir = os.path.join(out_dir, "tree_tmp_files", family)
    if os.path.isdir(tree_tmp_dir):
        shutil.rmtree(tree_tmp_dir)
    os.makedirs(tree_tmp_dir)

    # build consensus sequence from TELR TE output
    seqs = os.path.join(tree_tmp_dir, family + ".telr.fa")
    get_te_seq(te_seqs, family, consensus, add_consensus, seqs)

    # align TE sequences
    align_fa = os.path.join(tree_tmp_dir, family + ".align.fa")
    with open(align_fa, "w") as output:
        if thread > 8:
            mafft_thread = 8
        else:
            mafft_thread = thread
        subprocess.call(
            ["mafft", "--quiet", "--auto", "--thread", str(mafft_thread), seqs],
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


def per_base_len_dist(consensus_len, seq_len):
    return abs(seq_len - consensus_len) / consensus_len


def main():
    args = get_args()

    # TODO process consensus sequence, remove stuff after #
    consensus_dict = SeqIO.index(args.consensus, "fasta")

    # TODO add all family option

    meta_out = os.path.join(args.out, "meta.txt")

    families = args.family.replace(" ", "").replace("_", "-").split(",")
    # parse TELR output directories and load consensus
    all_te_seqs = []
    with open(meta_out, "w") as output:
        for telr_out_dir in args.telr_dirs:
            if telr_out_dir[-1] == "/":
                telr_out_dir = telr_out_dir[:-1]
            for telr_out_fa in glob.glob(telr_out_dir + "/**/*telr.fa", recursive=True):
                sample = (
                    os.path.basename(telr_out_fa)
                    .replace(".telr.fa", "")
                    .replace("_", "-")
                )
                # the following section can be parallized
                for record in SeqIO.parse(telr_out_fa, "fasta"):
                    family = record.id.split("#")[1].replace("_", "-")
                    if family in families:
                        consensus_seq = str(consensus_dict[family].seq)
                        consensus_length = len(consensus_seq)
                        sample_seq = str(record.seq)
                        seq_len = len(sample_seq)
                        len_diff = per_base_len_dist(consensus_length, seq_len)
                        div_diff = (
                            editdistance.eval(consensus_seq, sample_seq)
                            / consensus_length
                        )
                        contig = record.id.split("_")[0].replace("_", "-")
                        start = record.id.split("_")[1]
                        meta_out_line = "\t".join(
                            [
                                sample,
                                contig,
                                start,
                                family,
                                str(len_diff),
                                str(div_diff),
                            ]
                        )
                        output.write(meta_out_line + "\n")
                        if div_diff < args.divergence_filter:
                            record.id = "_".join([sample, contig, start, family])
                            all_te_seqs.append(record)

    for family in families:
        phylogeny_from_telr(
            te_seqs=all_te_seqs,
            family=family,
            consensus=consensus_dict,
            out_dir=args.out,
            method=args.method,
            bootstrap=args.bootstrap,
            thread=args.thread,
            add_consensus=args.add_consensus,
        )


main()
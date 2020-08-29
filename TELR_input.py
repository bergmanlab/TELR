import argparse
import sys
import os
import subprocess
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to detect TEs in long read data")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    # required
    required.add_argument("-i", "--reads", type=str,
                          help="reads in fasta/fastq format or read alignments in bam format", required=True)
    required.add_argument("-r", "--reference", type=str,
                          help="reference genome in fasta format", required=True)
    required.add_argument("-l", "--library", type=str,
                          help="TE consensus sequences in fasta format", required=True)

    # optional
    optional.add_argument("-x", "--presets", type=str,
                          help="parameter presets for different sequencing technologies (default = 'pacbio')", required=False)
    optional.add_argument("-p", "--polish", type=int,
                          help="rounds of contig polishing (default = 1)", required=False)
    optional.add_argument("-o", "--out", type=str,
                          help="directory to output data (default = '.')", required=False)
    optional.add_argument("-t", "--thread", type=int,
                          help="max cpu threads to use (default = '1')", required=False)
    optional.add_argument("-g", "--gap", type=int,
                          help="max gap size for flanking sequence alignment (default = '20')", required=False)
    optional.add_argument("-v", "--overlap", type=int,
                          help="max overlap size for flanking sequence alignment (default = '20')", required=False)
    parser._action_groups.append(optional)
    args = parser.parse_args()
    # TODO: remove intermediate files

    # checks if in files exist
    try:
        test = open(args.reads, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.reference, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.library, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    if args.presets is None:
        args.presets = "pacbio"

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.thread is None:
        args.thread = 1

    if args.polish is None:
        args.polish = 1

    if args.gap is None:
        args.gap = 20

    if args.overlap is None:
        args.overlap = 20

    return args


def parse_input(input_reads, input_reference, sample_name, out_dir):
    # create symbolic link for the input file
    input_reads_copy = out_dir + '/' + os.path.basename(input_reads)
    if not os.path.isabs(input_reads):
        input_reads = os.getcwd() + '/' + input_reads
    if os.path.islink(input_reads_copy):
        os.remove(input_reads_copy)
    try:
        os.symlink(input_reads, input_reads_copy)
    except Exception as e:
        print(e)
        sys.exit(1)

    input_reference_copy = out_dir + '/' + os.path.basename(input_reference)
    if not os.path.isabs(input_reference):
        input_reference = os.getcwd() + '/' + input_reference
    if os.path.islink(input_reference_copy):
        os.remove(input_reference_copy)
    try:
        os.symlink(input_reference, input_reference_copy)
    except Exception as e:
        print(e)
        sys.exit(1)

    sample_suffix = os.path.basename(input_reads_copy).split('.')[-1]
    if sample_suffix == 'bam' or sample_suffix == 'cram':
        print("Read alignment file is provided")
        skip_alignment = True
        tmp_fasta = out_dir + '/' + sample_name + '.tmp.fasta'
        bam2fasta(input_reads_copy, tmp_fasta)
        fasta = out_dir + '/' + sample_name + '.telr.fasta'
        rm_fasta_redundancy(tmp_fasta, fasta)
        os.remove(tmp_fasta)
    elif sample_suffix == 'fasta' or sample_suffix == 'fastq' or sample_suffix == 'fa' or sample_suffix == 'fq':
        print("Raw read file is provided")
        skip_alignment = False
        fasta = input_reads_copy
    else:
        print("input reads/alignments format not recognized, exiting...")
        sys.exit(1)

    return input_reads_copy, input_reference_copy, fasta, skip_alignment


def bam2fasta(bam, fasta):
    print("Converting bam to fasta...")
    try:
        with open(fasta, "w") as output:
            subprocess.call(["samtools", "fasta", bam], stdout=output)
    except Exception as e:
        print(e)
        print("Conversion failed, check input bam file, exiting...")
        sys.exit(1)


def rm_fasta_redundancy(fasta, new_fasta):
    records = set()
    with open(new_fasta, "w") as output_handle:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id not in records:
                records.add(record.id)
                SeqIO.write(record, output_handle, 'fasta')

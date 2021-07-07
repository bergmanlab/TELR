import argparse
import sys
import os
import subprocess
import logging
from Bio import SeqIO
from TELR_utility import mkdir


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to detect TEs in long read data"
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    # required
    required.add_argument(
        "-i",
        "--reads",
        type=str,
        help="reads in fasta/fastq format or read alignments in bam format",
        required=True,
    )
    required.add_argument(
        "-r",
        "--reference",
        type=str,
        help="reference genome in fasta format",
        required=True,
    )
    required.add_argument(
        "-l",
        "--library",
        type=str,
        help="TE consensus sequences in fasta format",
        required=True,
    )

    # optional
    optional.add_argument(
        "-m",
        "--method",
        type=str,
        help="method for read alignment, please provide 'nglmr' or 'minimap2' (default = 'nglmr')",
        required=False,
    )
    optional.add_argument(
        "-x",
        "--presets",
        type=str,
        help="parameter presets for different sequencing technologies, please provide 'ont' or 'pacbio' (default = 'pacbio')",
        required=False,
    )
    optional.add_argument(
        "-p",
        "--polish",
        type=int,
        help="rounds of contig polishing (default = 1)",
        required=False,
    )
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
        "-g",
        "--gap",
        type=int,
        help="max gap size for flanking sequence alignment (default = '20')",
        required=False,
    )
    optional.add_argument(
        "-v",
        "--overlap",
        type=int,
        help="max overlap size for flanking sequence alignment (default = '20')",
        required=False,
    )
    optional.add_argument(
        "--repeatmasker_family",
        action='store_true',
        help="If provided then repeatmasker will be used to annotate TE families in the assembled contigs (default: use minimap2 for contig TE annotation)",
        required=False,
    )
    optional.add_argument(
        "-k",
        "--keep_files",
        action='store_true',
        help="If provided then all intermediate files will be kept (default: remove intermediate files)",
        required=False,
    )
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # checks if in files exist
    try:
        test = open(args.reads, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.reads)
        sys.exit(1)

    try:
        test = open(args.reference, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.reference)
        sys.exit(1)

    try:
        test = open(args.library, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.library)
        sys.exit(1)

    if args.method is None:
        args.method = "nglmr"

    if args.presets is None:
        args.presets = "pacbio"

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    mkdir(args.out)

    if args.thread is None:
        args.thread = 1

    if args.polish is None:
        args.polish = 1

    if args.gap is None:
        args.gap = 20

    if args.overlap is None:
        args.overlap = 20

    return args


def parse_input(input_reads, input_reference, input_library, sample_name, out_dir):
    """
    Parse input files. If bam file is provided, convert to fasta format.
    """
    logging.info("Parsing input files...")
    # create symbolic link for the input file
    input_reads_copy = os.path.join(out_dir, os.path.basename(input_reads))
    if not os.path.isabs(input_reads):
        input_reads = os.path.abspath(input_reads)
    if os.path.islink(input_reads_copy):
        os.remove(input_reads_copy)
    try:
        os.symlink(input_reads, input_reads_copy)
    except Exception as e:
        print(e)
        logging.exception("Create symbolic link for " + input_reads + " failed")
        sys.exit(1)

    input_reference_copy = os.path.join(out_dir, os.path.basename(input_reference))
    if not os.path.isabs(input_reference):
        input_reference = os.path.abspath(input_reference)
    if os.path.islink(input_reference_copy):
        os.remove(input_reference_copy)
    try:
        os.symlink(input_reference, input_reference_copy)
    except Exception:
        logging.exception("Create symbolic link for " + input_reference + " failed")
        sys.exit(1)
    
    input_library_copy = os.path.join(out_dir, os.path.basename(input_library))
    if not os.path.isabs(input_library):
        input_library = os.path.abspath(input_library)
    if os.path.islink(input_library_copy):
        os.remove(input_library_copy)
    try:
        os.symlink(input_library, input_library_copy)
    except Exception:
        logging.exception("Create symbolic link for " + input_library + " failed")
        sys.exit(1)

    reads_filename, reads_extension = os.path.splitext(input_reads_copy)
    if reads_extension == ".bam":
        logging.info("BAM file is provided, skip alignment step")
        skip_alignment = True
        fasta = os.path.join(out_dir, sample_name + ".telr.fasta")
        logging.info("Converting input BAM file to fasta...")
        bam2fasta(input_reads_copy, fasta)
    elif (
        reads_extension == ".fasta"
        or reads_extension == ".fastq"
        or reads_extension == ".fa"
        or reads_extension == ".fq"
    ):
        logging.info("Raw reads are provided")
        skip_alignment = False
        fasta = input_reads_copy
    else:
        print("Input reads/alignments format not recognized, exiting...")
        logging.error("Input format not recognized")
        sys.exit(1)

    return input_reads_copy, input_reference_copy, input_library_copy, fasta, skip_alignment


def bam2fasta(bam, fasta):
    """
    Convert bam to fasta.
    """
    fasta_tmp = fasta + ".tmp"
    try:
        with open(fasta_tmp, "w") as output:
            subprocess.call(["samtools", "fasta", bam], stdout=output)
    except Exception as e:
        print(e)
        print("BAM to Fasta conversion failed, check input bam file, exiting...")
        sys.exit(1)

    try:
        rm_fasta_redundancy(fasta_tmp, fasta)
    except Exception as e:
        print(e)
        logging.exception("Remove redundancy in fasta file failed")
        sys.exit(1)
    os.remove(fasta_tmp)


def rm_fasta_redundancy(fasta, new_fasta):
    """
    Remove redundancy in fasta file.
    If there are multiple IDs, keep the first one.
    """
    records = set()
    with open(new_fasta, "w") as output_handle:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id not in records:
                records.add(record.id)
                SeqIO.write(record, output_handle, "fasta")

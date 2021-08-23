import argparse
import sys
import os
import subprocess
import logging
from Bio import SeqIO
from telr.TELR_utility import mkdir


def get_args():
    parser = argparse.ArgumentParser(
        description="Program for detecting non-reference TEs in long read data"
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
        "--aligner",
        type=str,
        help="choose method for read alignment, please provide 'nglmr' or 'minimap2' (default = 'nglmr')",
        required=False,
    )
    optional.add_argument(
        "--assembler",
        type=str,
        help="Choose the method to be used for local contig assembly step, please provide 'wtdbg2' or 'flye' (default = 'wtdbg2')",
        required=False,
    )
    optional.add_argument(
        "--polisher",
        type=str,
        help="Choose the method to be used for local contig polishing step, please provide 'wtdbg2' or 'flye' (default = 'wtdbg2')",
        required=False,
    )
    optional.add_argument(
        "-x",
        "--presets",
        type=str,
        help="parameter presets for different sequencing technologies, please provide 'pacbio' or 'ont' (default = 'pacbio')",
        required=False,
    )
    optional.add_argument(
        "-p",
        "--polish_iterations",
        type=int,
        help="iterations of contig polishing (default = 1)",
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
        "--flank_len",
        type=int,
        help="flanking sequence length (default = '500')",
        required=False,
    )
    optional.add_argument(
        "--af_flank_interval",
        type=int,
        help="5' and 3'flanking sequence interval size used for allele frequency estimation (default = '100')",
        required=False,
    )
    optional.add_argument(
        "--af_flank_offset",
        type=int,
        help="5' and 3' flanking sequence offset size used for allele frequency estimation (default = '200')",
        required=False,
    )
    optional.add_argument(
        "--af_te_interval",
        type=int,
        help="5' and 3' te sequence interval size used for allele frequency estimation (default: '50')",
        required=False,
    )
    optional.add_argument(
        "--af_te_offset",
        type=int,
        help="5' and 3' te sequence offset size used for allele frequency estimation (default: '50')",
        required=False,
    )
    optional.add_argument(
        "--different_contig_name",
        action="store_true",
        help="If provided then TELR does not require the contig name to match before and after annotation liftover (default: require contig name to be the same before and after liftover)",
        required=False,
    )
    optional.add_argument(
        "--minimap2_family",
        action="store_true",
        help="If provided then minimap2 will be used to annotate TE families in the assembled contigs (default: use repeatmasker for contig TE annotation)",
        required=False,
    )
    optional.add_argument(
        "-k",
        "--keep_files",
        action="store_true",
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

    # check if optional arguments are valid
    if args.aligner is None:
        args.aligner = "nglmr"
    elif args.aligner not in ["nglmr", "minimap2"]:
        print("Please provide a valid alignment method (nglmr/minimap2), exiting...")
        sys.exit(1)

    if args.assembler is None:
        args.assembler = "wtdbg2"
    elif args.assembler not in ["wtdbg2", "flye"]:
        print("Please provide a valid assembly method (wtdbg2/flye), exiting...")
        sys.exit(1)

    if args.polisher is None:
        args.polisher = "wtdbg2"
    elif args.polisher not in ["wtdbg2", "flye"]:
        print("Please provide a valid polish method (wtdbg2/flye), exiting...")
        sys.exit(1)

    if args.presets is None:
        args.presets = "pacbio"
    elif args.presets not in ["pacbio", "ont"]:
        print("Please provide a valid preset option (pacbio/ont), exiting...")
        sys.exit(1)

    if args.polish_iterations is None:
        args.polish_iterations = 1
    elif args.polish_iterations < 1:
        print("Please provide a valid number of iterations for polishing, exiting...")

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    mkdir(args.out)

    if args.thread is None:
        args.thread = 1

    if args.flank_len is None:
        args.flank_len = 500

    if args.af_flank_interval is None:
        args.af_flank_interval = 100
    else:
        if args.af_flank_interval <= 0:
            print(
                "Please provide a valid flanking sequence interval size (positive integer) for allele frequency estimation, exiting..."
            )
            sys.exit(1)

    if args.af_flank_offset is None:
        args.af_flank_offset = 200
    else:
        if args.af_flank_offset < 0:
            print(
                "Please provide a valid flanking sequence offset size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.af_te_interval is None:
        args.af_te_interval = 50
    else:
        if args.af_te_interval <= 0:
            print(
                "Please provide a valid TE interval size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.af_te_offset is None:
        args.af_te_offset = 50
    else:
        if args.af_te_offset < 0:
            print(
                "Please provide a valid TE offset size (positive integer) for allele frequency estimation, exiting..."
            )

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

    return (
        input_reads_copy,
        input_reference_copy,
        input_library_copy,
        fasta,
        skip_alignment,
    )


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

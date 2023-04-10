import os
import sys
import logging



def input_reads(out_dir, in_reads):
    reads_extension = in_reads[in_reads.rindex("."):]
    if(reads_extension in [".fasta",".fastq",".fa",".fq"]):
        logging.info("Raw reads are provided")
        in_reads_copy = os.path.join(out_dir, "reads.fasta")
    elif(reads_extension == ".bam"):
        logging.info("BAM file is provided, skip alignment step")
        in_reads_copy = os.path.join(out_dir, "reads.bam")
    else:
        print("Input reads/alignments format not recognized, exiting...")
        logging.error("Input format not recognized")
        sys.exit(1)
    symlink(in_reads, in_reads_copy)

def input_library(out_dir, in_lib):
    if not in_lib[in_lib.rindex("."):] in [".fasta",".fastq",".fa",".fq"]:
        print("Input library format not recognized, exiting...")
        logging.error("Input format not recognized")
        sys.exit(1)
    symlink(in_lib, os.path.join(out_dir,"library.fasta"))

def input_reference(out_dir, in_ref):
    if not in_ref[in_ref.rindex("."):] in [".fasta",".fastq",".fa",".fq"]:
        print("Input library format not recognized, exiting...")
        logging.error("Input format not recognized")
        sys.exit(1)
    symlink(in_ref, os.path.join(out_dir,"reference.fasta"))

def symlink(input, output):
    if os.path.islink(output):
        os.remove(output)
    try:
        os.symlink(input, output)
    except Exception:
        logging.exception(f"Create symbolic link for {input} failed")
        sys.exit(1)

if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])
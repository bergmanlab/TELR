import os
import sys
import logging
import subprocess
#from Bio import SeqIO

def input(file_type, sample_name, in_file):
    print(f"---- {file_type} ----")
    extension = in_file[in_file.rindex("."):]
    extension_matrix = { #valid file extensions for each type of input
        "reads":[".fasta",".fastq",".fa",".fq",".bam"],
        "library":[".fasta",".fastq",".fa",".fq"],
        "reference":[".fasta",".fastq",".fa",".fq"]
        }[file_type]
    if extension not in extension_matrix: 
        print(f"Input {file_type} format not recognized, exiting...")
        logging.error("Input format not recognized")
        sys.exit(1)
    print(f"intermediate_files/input/{file_type}-{sample_name}{extension}")
    symlink(in_file, f"intermediate_files/input/{file_type}-{sample_name}{extension}")
    if extension == ".bam":
        bam2fasta(f"intermediate_files/input/reads-{sample_name}{extension}",f"intermediate_files/input/reads-{sample_name}.fasta")

def symlink(input, output): #create a symbolic link at the output location referencing the input path.
    if os.path.islink(output):
        os.remove(output)
    try:
        os.symlink(input, output)
    except Exception:
        logging.exception(f"Create symbolic link for {input} failed")
        sys.exit(1)

def bam2fasta(bam, fasta):
    """
    Convert bam to fasta.
    """
    fasta_tmp = f"{fasta}.tmp"
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

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
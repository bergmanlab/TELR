import sys
import os
import subprocess
import logging
import time
from Bio import SeqIO
from STELR_utility import format_time

def alignment(read, reference, outfile, sample_name, thread, method, presets):
    """
    This function takes raw reads and performs read alignment using ngmlr or minimap2.
    """
    logging.info("Start alignment...")
    start_time = time.perf_counter()
    method_array = {
        "ngmlr":{"ont":{"presets":"ont","label":"ont"},"pacbio":{"presets":"pacbio","label":"pb"}},
        "minimap2":{"ont":{"presets":"map-ont","label":""},"pacbio":{"presets":"map-pb","label":""}} #still needs label to prevent error in line 27
        }
    if(method in method_array): method_array = method_array[method]
    else: 
        print(f"Read method not recognized, please provide {string_list(method_array)}, exiting...")
        sys.exit(1)
    if(presets in method_array):
        method_array = method_array[presets]
    else: 
        print(f"Read presets not recognized, please provide {string_list(method_array)}, exiting...")
        sys.exit(1)
    try:
        command_array = {
            "ngmlr":["ngmlr","-r",reference,"-q",read,"-x",method_array["presets"],"-t",str(thread),"--rg-id",sample_name,"--rg-sm",sample_name,"--rg-lb",method_array["label"],"--no-progress"],
            "minimap2":["minimap2","--cs","--MD","-Y","-L","-ax",method_array["presets"],reference,read]
            }[method]
        with open(outfile, "w") as output:
            subprocess.call(command_array,stdout=output)
    except Exception as e:
        print(e)
        print("Read alignment failed, check input reads, exiting...")
        sys.exit(1)
    proc_time = time.perf_counter() - start_time
    logging.info(f"First alignment finished in {format_time(proc_time)}")
    
def string_list(list1):
    if list1 is dict:
        list1 = list1.keys()
    if len(list1) < 2:
        return str(list1[0])
    elif len(list1) == 2:
        return " or ".join(list1)
    else:
        return ", ".join(list1[:-1]) + ", or " + list1[-1] 


def sort_index_bam(unsorted_bam, sorted_bam, thread):
    """
    Sort and index bam file
    """
    logging.info("Sort and index BAM...")
    start_time = time.perf_counter()
    try:
        print(["samtools", "sort", "-@", str(thread), "-o", sorted_bam, unsorted_bam])
        subprocess.call(["samtools", "sort", "-@", str(thread), "-o", sorted_bam, unsorted_bam])
        subprocess.call(["samtools", "index", "-@", str(thread), sorted_bam])
    except Exception as e:
        print(e)
        print("Sort and index BAM file failed, exiting...")
        sys.exit(1)
    if os.path.isfile(sorted_bam) is False:
        sys.stderr.write("Sorted and indexed BAM file does not exist, exiting...\n")
        sys.exit(1)
    os.remove(unsorted_bam)

    proc_time = time.perf_counter() - start_time
    logging.info(f"Bam sorted and indexed in {format_time(proc_time)}")

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
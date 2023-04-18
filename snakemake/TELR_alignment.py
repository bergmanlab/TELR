import sys
import os
import subprocess
import logging
import time
from telr.TELR_utility import format_time

def alignment(read, reference, out, sample_name, thread, method, presets):
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
    if(presets in method_array["presets"]):
        method_array = method_array[presets]
    else: 
        print(f"Read presets not recognized, please provide {string_list(method_array)}, exiting...")
        sys.exit(1)
    try:
        align_sam = f"{out}/{sample_name}.sam"
        command_array = {
            "ngmlr":["ngmlr","-r",reference,"-q",read,"-x",method_array["presets"],"-t",str(thread),"--rg-id",sample_name,"--rg-sm",sample_name,"--rg-lb",method_array["label"],"--no-progress"],
            "minimap2":["--cs","--MD","-Y","-L","-ax",method_array["presets"],reference,read]
            }[method]
        with open(align_sam, "w") as output:
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


def alignmenta(bam, read, reference, out, sample_name, thread, method, presets):
    """
    This function takes raw reads and performs read alignment using ngmlr or minimap2.
    """

    sort_index_bam(align_sam, bam, thread)
    if os.path.isfile(bam) is False:
        sys.stderr.write("Sorted and indexed BAM file does not exist, exiting...\n")
        sys.exit(1)
    os.remove(align_sam)

    proc_time = time.time() - start_time
    logging.info("First alignment finished in " + format_time(proc_time))


def sort_index_bam(unsorted_bam, sorted_bam, thread):
    """
    Sort and index bam file
    """
    logging.info("Sort and index BAM...")
    start_time = time.perf_counter()
    try:
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

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
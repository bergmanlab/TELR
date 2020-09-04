import sys
import os
import subprocess
import logging
import time
from TELR_utility import format_time


def alignment(bam, read, reference, out, sample_name, thread, presets):
    """
    This function takes raw reads and performs read alignment using ngmlr.
    """
    if presets == "ont":
        presets = "ont"
        label = "ont"
    else:
        presets = "pacbio"
        label = "pb"
    bam_tmp = out + "/" + sample_name + ".tmp.bam"

    logging.info("Start alignment...")
    start_time = time.time()

    try:
        with open(bam_tmp, "w") as output:
            subprocess.call(
                [
                    "ngmlr",
                    "-r",
                    reference,
                    "-q",
                    read,
                    "-x",
                    presets,
                    "-t",
                    str(thread),
                    "--rg-id",
                    sample_name,
                    "--rg-sm",
                    sample_name,
                    "--rg-lb",
                    label,
                    "--no-progress",
                ],
                stdout=output,
            )
    except Exception as e:
        print(e)
        print("Read alignment (ngmlr) failed, check input reads, exiting...")
        sys.exit(1)
    proc_time = time.time() - start_time
    logging.info("Alignment finished in " + format_time(proc_time))

    sort_index_bam(bam_tmp, bam, thread)
    if os.path.isfile(bam) is False:
        sys.stderr.write("Sorted and indexed BAM file does not exist, exiting...\n")
        sys.exit(1)
    os.remove(bam_tmp)


def sort_index_bam(bam, sorted_bam, thread):
    """
    Sort and index bam file
    """
    logging.info("Sort and index BAM...")
    try:
        subprocess.call(["samtools", "sort", "-@", str(thread), "-o", sorted_bam, bam])
        subprocess.call(["samtools", "index", "-@", str(thread), sorted_bam])
    except Exception as e:
        print(e)
        print("Sort and index BAM file failed, exiting...")
        sys.exit(1)

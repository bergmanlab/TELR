import sys
import os
import subprocess
import logging
import time
from TELR_utility import format_time


def alignment(bam, read, reference, out, sample_name, thread, method, presets):
    """
    This function takes raw reads and performs read alignment using ngmlr or minimap2.
    """
    logging.info("Start alignment...")
    start_time = time.time()
    if method == "nglmr":
        if presets == "ont":
            presets = "ont"
            label = "ont"
        elif presets == "pacbio":
            presets = "pacbio"
            label = "pb"
        else:
            print(
                "Read presets not recognized, please provide ont or pacbio, exiting..."
            )
            sys.exit(1)

        try:
            bam_tmp = out + "/" + sample_name + ".tmp.bam"
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
            print("Read alignment failed, check input reads, exiting...")
            sys.exit(1)
    elif method == "minimap2":
        if presets == "ont":
            presets = "map-ont"
        elif presets == "pacbio":
            presets = "map-pb"
        else:
            print(
                "Read presets not recognized, please provide ont or pacbio, exiting..."
            )
            sys.exit(1)
        try:
            sam = out + "/" + sample_name + ".sam"
            with open(sam, "w") as output:
                subprocess.call(
                    [
                        "minimap2",
                        "--cs",
                        "--MD",
                        "-Y",
                        "-L",
                        "-ax",
                        presets,
                        reference,
                        read,
                    ],
                    stdout=output,
                )
            bam_tmp = out + "/" + sample_name + ".tmp.bam"
            with open(bam_tmp, "w") as output:
                subprocess.call(["samtools", "view", "-bS", sam], stdout=output)
                os.remove(sam)
        except Exception as e:
            print(e)
            print("Read alignment failed, check input reads, exiting...")
            sys.exit(1)
    else:
        print(
            "Alignment method not recognized, please provide ont or pacbio, exiting..."
        )
        sys.exit(1)

    sort_index_bam(bam_tmp, bam, thread)
    if os.path.isfile(bam) is False:
        sys.stderr.write("Sorted and indexed BAM file does not exist, exiting...\n")
        sys.exit(1)
    os.remove(bam_tmp)

    proc_time = time.time() - start_time
    logging.info("First alignment finished in " + format_time(proc_time))


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

import sys
import os
import subprocess
import logging
import time
from telr.TELR_utility import format_time


def alignment(files, out, thread, method, presets):
    """
    This function takes raw reads and performs read alignment using ngmlr or minimap2.
    """
    files.frame = "alignment()"
    logging.info("Start alignment...")
    start_time = time.perf_counter()
    label_presets = {"pacbio":{"ngmlr":"pb","minimap2":"map-pb"},"ont":{"ngmlr":"ont","minimap2":"map-ont"}}
    if(not presets in label_presets):
        print(
            "Read presets not recognized, please provide ont or pacbio, exiting..."
        )
        sys.exit(1)
    method_array = {
            #ngmlr version 0.2.7 | current 0.2.7
            "ngmlr":{
                "file_extension":".tmp.sam",
                "run_array":[
                    "ngmlr",
                    "-r", files.reference.path,
                    "-q", files.reads.path,
                    "-x", presets,
                    "-t", str(thread),
                    "--rg-id", files.sample_name,
                    "--rg-sm", files.sample_name,
                    "--rg-lb", label_presets[presets]["ngmlr"],
                    "--no-progress"
                    ]
                },
            #minimap2 verion 2.22 | current 2.24
            "minimap2":{
                "file_extension":".sam",
                "run_array":[
                    "minimap2",
                    "--cs",
                    "--MD",
                    "-Y",
                    "-L",
                    "-ax",
                    label_presets[presets]["minimap2"],
                    files.reference.path,
                    files.reads.path,
                    ]
                }
            }
    if(not method in method_array):
        print(
            "Alignment method not recognized, please provide ont or pacbio, exiting..."
        )
        sys.exit(1)
    else: method_array = method_array[method]
    try:
        files.add("align_sam", out, method_array["file_extension"])
        with files.align_sam.open("w") as output:
            subprocess.call(
                method_array["run_array"],
                stdout=output,
            )        
    except Exception as e:
        print(e)
        print("Read alignment failed, check input reads, exiting...")
        sys.exit(1)

    sort_index_bam(files.align_sam.path, files.bam.path, thread)
    if files.bam.exists() is False:
        sys.stderr.write("Sorted and indexed BAM file does not exist, exiting...\n")
        sys.exit(1)
    files.align_sam.remove()

    proc_time = time.perf_counter() - start_time
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

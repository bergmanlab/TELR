import sys
import os
import subprocess


def alignment(read, reference, out, sample_name, thread, presets):
    if presets == "ont":
        presets_nglmr = "ont"
    else:
        presets_nglmr = "pacbio"
    bam = out + "/" + sample_name + ".tmp.bam"
    print("Generating raw read alignment...")
    try:
        with open(bam, "w") as output:
            subprocess.call(["ngmlr",
                             "-r", reference,
                             "-q", read,
                             "-x", presets_nglmr,
                             "-t", str(thread),
                             "--rg-id", sample_name,
                             "--rg-sm", sample_name,
                             "--rg-lb", "pb", "--no-progress"], stdout=output)
            # TODO: change rg-lb
    except Exception as e:
        print(e)
        print("Read alignment (ngmlr) failed, check input reads, exiting...")
        sys.exit(1)

    sorted_bam = out + "/" + sample_name + "_sort.bam"
    sort_index_bam(bam, sorted_bam, thread)

    if os.path.isfile(sorted_bam) is False:
        sys.stderr.write(
            "BAM file from read alignment step failed, exiting....\n")
        sys.exit(1)
    else:
        return(sorted_bam)


def sort_index_bam(bam, sorted_bam, thread):
    print("Sort and index read alignment...")
    try:
        subprocess.call(["samtools", "sort", "-@",
                         str(thread), "-o", sorted_bam, bam])
        subprocess.call(["samtools", "index", "-@", str(thread), sorted_bam])
    except Exception as e:
        print(e)
        print("Sort and index bam file failed, exiting....")
        sys.exit(1)
    print("Done\n")

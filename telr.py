#!/usr/bin/env python3

import sys
import os
import time
import logging
from TELR_input import get_args, parse_input
from TELR_alignment import alignment, sort_index_bam
from TELR_sv import detect_sv, vcf_parse_filter
from TELR_assembly import prep_assembly, local_assembly
from TELR_te import annotate_contig, find_te
from TELR_utility import format_time

# python3 telr.py -o $output_dir -i $read_path -r $reference_path -l $te_library_path -t 16 -x pacbio


def main():
    # logging config
    formatstr = "%(asctime)s: %(levelname)s: %(message)s"
    datestr = "%m/%d/%Y %H:%M:%S"
    logging.basicConfig(level=logging.DEBUG,
                        filename="TELR.log",
                        filemode="w",
                        format=formatstr,
                        datefmt=datestr)
    logging.info("CMD: " + ' '.join(sys.argv))
    start_time = time.time()

    # Parse input
    args = get_args()
    sample_name = os.path.splitext(os.path.basename(args.reads))[0]
    reads, reference, fasta, skip_alignment = parse_input(
        args.reads, args.reference, sample_name, args.out)

    # Alignment
    bam = os.path.join(args.out, sample_name + "_sort.bam")
    if not skip_alignment:
        alignment(bam, fasta, reference, args.out,
                  sample_name, args.thread, args.presets)
    else:
        sort_index_bam(reads, bam, args.thread)

    # Detect and parse SV
    vcf_parsed = os.path.join(args.out, sample_name+".vcf_filtered.tsv")
    detect_sv(vcf_parsed, bam, reference, args.library,
              args.out, sample_name, args.thread)

    # Local assembly
    contig_dir = os.path.join(args.out, "contig_assembly")
    local_assembly(contig_dir, vcf_parsed, args.out, sample_name, fasta,
                   args.thread, args.presets, args.polish)

    # find TEs
    # TODO: different approaches
    find_te(contig_dir, vcf_parsed, reference, args.library, args.out,
            sample_name, args.thread, args.gap, args.overlap, args.presets)

    proc_time = time.time() - start_time

    print("TELR finished!")
    logging.info("TELR finished in " + format_time(proc_time))


main()

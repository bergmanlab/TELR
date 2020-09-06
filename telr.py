#!/usr/bin/env python3

import sys
import os
import time
import logging
from TELR_input import get_args, parse_input
from TELR_alignment import alignment, sort_index_bam
from TELR_sv import detect_sv, vcf_parse_filter
from TELR_assembly import local_assembly
from TELR_te import annotate_contig, find_te
from TELR_utility import format_time

# python3 telr.py -o $output_dir -i $reads -r $reference -l $te_library


def main():
    args = get_args()
    # logging config
    formatstr = "%(asctime)s: %(levelname)s: %(message)s"
    datestr = "%m/%d/%Y %H:%M:%S"
    logging.basicConfig(
        level=logging.DEBUG,
        filename=os.path.join(args.out, "TELR.log"),
        filemode="w",
        format=formatstr,
        datefmt=datestr,
    )
    logging.info("CMD: " + " ".join(sys.argv))
    start_time = time.time()

    # Parse input
    sample_name = os.path.splitext(os.path.basename(args.reads))[0]
    reads, reference, fasta, skip_alignment = parse_input(
        args.reads, args.reference, sample_name, args.out
    )

    # # Alignment
    bam = os.path.join(args.out, sample_name + "_sort.bam")
    if not skip_alignment:
        alignment(
            bam, fasta, reference, args.out, sample_name, args.thread, args.presets
        )
    else:
        sort_index_bam(reads, bam, args.thread)

    # initialize loci eveluation file
    loci_eval = os.path.join(args.out, sample_name + ".loci_eval.tsv")
    if os.path.isfile(loci_eval):
        os.remove(loci_eval)

    # Detect and parse SV
    vcf = os.path.join(args.out, sample_name + ".vcf")
    detect_sv(vcf, bam, reference, args.library, args.out, sample_name, args.thread)

    # Parse SV and filter for TE candidate locus
    vcf_parsed = os.path.join(args.out, sample_name + ".vcf_filtered.tsv")
    vcf_parse_filter(
        vcf,
        vcf_parsed,
        bam,
        args.library,
        args.out,
        sample_name,
        args.thread,
        loci_eval,
    )

    # Local assembly
    contig_dir = os.path.join(args.out, "contig_assembly")
    local_assembly(
        contig_dir,
        vcf_parsed,
        args.out,
        sample_name,
        fasta,
        args.thread,
        args.presets,
        args.polish,
    )

    # Annotate contig for TE region
    (
        contig_te_annotation,
        contig_rm_annotation,
        te_freq,
        te_fa,
        merge_contigs,
    ) = annotate_contig(
        contig_dir,
        args.library,
        vcf_parsed,
        args.out,
        sample_name,
        args.thread,
        args.presets,
        loci_eval,
    )

    # find TEs
    find_te(
        contig_te_annotation,
        contig_rm_annotation,
        te_freq,
        te_fa,
        merge_contigs,
        vcf_parsed,
        reference,
        args.out,
        sample_name,
        args.gap,
        args.overlap,
        args.presets,
        loci_eval
    )

    proc_time = time.time() - start_time

    print("TELR finished!")
    logging.info("TELR finished in " + format_time(proc_time))


main()

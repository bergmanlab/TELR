#!/usr/bin/env python3

import sys
import os
import time
import logging
import shutil
from telr.TELR_input import get_args, parse_input
from telr.TELR_alignment import alignment, sort_index_bam
from telr.TELR_sv import detect_sv, vcf_parse_filter
from telr.TELR_assembly import get_local_contigs
from telr.TELR_te import annotate_contig, find_te, get_af, repeatmask, gff3tobed
from telr.TELR_output import generate_output
from telr.TELR_utility import format_time, mkdir, export_env

"""
Author: Shunhua Han <hanshunhua0829@gmail.com>
Contributors: Casey Bergman <cbergman.uga.edu>, Guilherme Dias <guilhermebdias@live.com>
"""


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

    # create directory for intermediate files
    tmp_dir = os.path.join(args.out, "intermediate_files")
    mkdir(tmp_dir)

    # Parse input
    sample_name = os.path.splitext(os.path.basename(args.reads))[0]
    reads, reference, library, fasta, skip_alignment = parse_input(
        args.reads, args.reference, args.library, sample_name, tmp_dir
    )

    # # Alignment
    bam = os.path.join(tmp_dir, f"{sample_name}_sort.bam")
    if not skip_alignment:
        alignment(
            bam,
            fasta,
            reference,
            tmp_dir,
            sample_name,
            args.thread,
            args.aligner,
            args.presets,
        )
    else:
        sort_index_bam(reads, bam, args.thread)

    # initialize loci eveluation file
    loci_eval = os.path.join(args.out, f"{sample_name}.loci_eval.tsv")
    if os.path.isfile(loci_eval):
        os.remove(loci_eval)

    # Detect and parse SV
    vcf = os.path.join(tmp_dir, sample_name + ".vcf")
    detect_sv(vcf, bam, reference, tmp_dir, sample_name, args.thread)

    # Parse SV and filter for TE candidate locus
    vcf_parsed = os.path.join(tmp_dir, f"{sample_name}.vcf_filtered.tsv")
    vcf_parse_filter(
        vcf,
        vcf_parsed,
        bam,
        library,
        tmp_dir,
        sample_name,
        args.thread,
        loci_eval,
    )

    # Local assembly
    contig_dir = os.path.join(tmp_dir, "contig_assembly")
    merged_contigs, assembly_passed_loci = get_local_contigs(
        assembler=args.assembler,
        polisher=args.polisher,
        contig_dir=contig_dir,
        vcf_parsed=vcf_parsed,
        out=tmp_dir,
        sample_name=sample_name,
        bam=bam,
        raw_reads=fasta,
        thread=args.thread,
        presets=args.presets,
        polish_iterations=args.polish_iterations,
    )

    # Annotate contig for TE region
    contig_te_annotation, te_fa = annotate_contig(
        merged_contigs,
        assembly_passed_loci,
        library,
        vcf_parsed,
        tmp_dir,
        sample_name,
        args.thread,
        args.presets,
        args.minimap2_family,
        loci_eval,
    )

    # calculate AF
    te_freq = get_af(
        tmp_dir,
        sample_name,
        bam,
        fasta,
        contig_te_annotation,
        contig_dir,
        vcf_parsed,
        args.af_flank_interval,
        args.af_flank_offset,
        args.af_te_interval,
        args.af_te_offset,
        args.presets,
        args.thread,
    )

    # repeatmask reference genome using custom TE library
    repeatmask_ref_dir = os.path.join(tmp_dir, "ref_repeatmask")
    ref_masked, te_gff = repeatmask(
        ref=reference,
        library=library,
        outdir=repeatmask_ref_dir,
        thread=args.thread,
    )
    ref_te_bed = os.path.join(tmp_dir, f"{os.path.basename(reference)}.te.bed")
    if te_gff is not None:
        gff3tobed(te_gff, ref_te_bed)
    else:
        ref_te_bed = None

    # find TEs
    liftover_json = find_te(
        reference=reference,
        contigs_fa=merged_contigs,
        contig_te_bed=contig_te_annotation,
        ref_te_bed=ref_te_bed,
        out=tmp_dir,
        gap=args.gap,
        overlap=args.overlap,
        flank_len=args.flank_len,
        different_contig_name=args.different_contig_name,
        keep_files=args.keep_files,
        thread=args.thread,
    )

    # generate output files
    if liftover_json:
        generate_output(
            liftover_report_path=liftover_json,
            te_freq_dict=te_freq,
            te_fa=te_fa,
            vcf_parsed=vcf_parsed,
            contig_te_annotation=contig_te_annotation,
            contig_fa=merged_contigs,
            out=args.out,
            sample_name=sample_name,
            ref=reference,
        )
    else:
        print("No non-reference TE insertion found")
        logging.info("TELR found no non-reference TE insertions")

    # clean tmp files
    if not args.keep_files:
        shutil.rmtree(tmp_dir)
    os.remove(loci_eval)

    # export conda environment
    env_file = os.path.join(args.out, "conda_env.yml")
    export_env(env_file)

    proc_time = time.time() - start_time
    print("TELR finished!")
    logging.info("TELR finished in " + format_time(proc_time))


if __name__ == "__main__":
    main()

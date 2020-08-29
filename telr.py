#!/usr/bin/env python3

import sys
import os
import time

from TELR_input import get_args, parse_input
from TELR_alignment import alignment, sort_index_bam
from TELR_sv import detect_sv, sniffle_parse
from TELR_assembly import local_assembly
from TELR_te import annotate_contig, find_te
from TELR_utility import rm_file, mkdir

# python3 telr.py -o $output_dir -i $read_path -r $reference_path -l $te_library_path -t 16 -x pacbio

def main():
    args = get_args()
    print("CMD: ", ' '.join(sys.argv), "\n")
    sample_name = os.path.splitext(os.path.basename(args.reads))[0]
    reads, reference, fasta, skip_alignment = parse_input(args.reads, args.reference, sample_name, args.out)

    # Alignment
    if not skip_alignment:
        start_time = time.time()
        bam = alignment(fasta, reference, args.out, sample_name, args.thread, args.presets)
        proc_time = time.time() - start_time
        print("Alignment time:", proc_time, "\n")
    else:
        print("Alignment provided, sort and index BAM...")
        bam = args.out + "/" + sample_name + "_sort.bam"
        sort_index_bam(reads, bam, args.thread)

    # SV detection
    start_time = time.time()
    vcf = detect_sv(bam, reference, args.out, sample_name, args.thread)
    proc_time = time.time() - start_time
    print("SV detection time:", proc_time, "\n")

    # Sniffle output parsing
    # vcf = args.out+"/"+sample_name+".vcf"
    start_time = time.time()
    vcf_parsed, contig_reads_dir = sniffle_parse(vcf, args.out, sample_name, fasta, args.library, args.thread)
    proc_time = time.time() - start_time
    print("SV parsing time:", proc_time, "\n")

    # local assembly on every insertion cluster
    # contig_reads_dir=args.out+"/"+"contig_reads"
    # vcf_parsed=args.out+"/"+sample_name+".vcf.parsed.filtered"
    start_time = time.time()
    contig_assembly_dir = local_assembly(contig_reads_dir, vcf_parsed, args.out, fasta, args.thread, args.presets, args.polish)
    proc_time = time.time() - start_time
    print("Local assembly time:", proc_time, "\n")
    
    # annotate TEs using different method, extract flanking sequence and annotate TE family
    # contig_assembly_dir=args.out+"/"+"contig_assembly"
    # vcf_parsed=args.out+"/"+sample_name+".vcf.parsed.filtered"
    start_time = time.time()
    te_contigs_annotation, family_annotation, te_freq, te_fa, merge_contigs = annotate_contig(contig_assembly_dir, args.library, vcf_parsed, args.out, sample_name, args.thread, args.presets)
    proc_time = time.time() - start_time
    print("Contig annotation time:", proc_time, "\n")

    # map contig back to the reference
    # te_contigs_annotation = args.out+"/"+sample_name+".te2contig_filter.bed"
    # family_annotation = args.out+"/"+sample_name+".te2contig_rm.merge.bed"
    # te_fa=args.out+"/"+sample_name+".te.fa"
    # merge_contigs=args.out+"/"+sample_name+".contigs.fa"
    start_time = time.time()
    find_te(merge_contigs, reference, te_contigs_annotation, family_annotation, te_freq, te_fa, args.out, sample_name, args.gap, args.overlap, args.presets)
    proc_time = time.time() - start_time
    print("Contig flanks liftover time:", proc_time, "\n")

    print("TELR finished!")

main()
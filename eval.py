#!/usr/bin/env python3
import sys
import argparse
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
import pandas as pd
import re
import liftover
import glob

'''
A script to evaluate TELR predictions by using two different strategies

1. lift over TE annotations from query genome to subject genome, then comprare lifted annotations with ones from TELR prediction.

2. map whole/flanking sequences from contig back to query genome and confirm all the TE predictions overlap the query genome TE annotation.

how to run: python3 eval.py -1 $genome1 -2 $genome2 -a $annotation -p $pred -o $out_dir
'''

def get_args():
    parser = argparse.ArgumentParser(description="Script to evaluate TELR predictions")

    ## required ##
    parser.add_argument("-1", "--fasta1", type=str, help="genome 1", required=True)
    parser.add_argument("-2", "--fasta2", type=str, help="genome 2", required=True)
    parser.add_argument("-a", "--annotation1", type=str, help="annotation file for genome 1 in bed/gff format (default = 'bed')", required=True)
    parser.add_argument("-b", "--annotation2", type=str, help="annotation file for genome 2 in bed/gff format (default = 'bed')", required=True)
    parser.add_argument("-t", "--telr_dir", type=str, help="TELR prediction in bed format", required=True)
    ## optional ##
    parser.add_argument("-o", "--out", type=str, help="directory to output data (default = '.')", required=False)
    parser.add_argument("--region1", type=str, help="focus region in genome1 in the analysis (default = 'none')", required=False)
    parser.add_argument("--region2", type=str, help="focus region in genome2 in the analysis (default = 'none')", required=False)
    parser.add_argument("-d", "--discard_family", type=str, help="discard families in the analysis, separated by comma (default = 'none')", required=False)
    parser.add_argument("-x", "--presets", type=str, help="parameter presets for different sequencing technologies: pacbio or ont (default = 'pacbio')", required=False)

    args = parser.parse_args()

     # checks if in files exist
    if args.presets is None:
        args.presets = "pacbio"

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    
    return args
    

def main():
    args = get_args()

    sample_name = os.path.splitext(os.path.basename(args.fasta1))[0]

    if args.presets == "ont":
        preset = "map-ont"
    elif args.presets == "pacbio":
        preset = "map-pb"
    else:
        print("unrecognized preset provided, used pacbio mode")
        preset = "map-pb"
    overlap = 20
    gap = 20

    stat_summary = args.out + "/" + "stat.summary.txt"
    rm(stat_summary)

    # telr output files
    telr_bed = []
    for file in glob.glob(args.telr_dir + "/*.final.bed"):
        telr_bed.append(file)
    prefix = os.path.basename(telr_bed[0]).replace('.final.bed', '')
    telr_bed = args.telr_dir + "/" + prefix + '.final.bed'
    telr_contigs = args.telr_dir + "/" + prefix + '.contigs.fa'
    telr_bed = args.telr_dir + '/' + prefix + ".final.bed"
    telr_meta = args.telr_dir + '/' + prefix + ".final.meta.bed"
    telr_contig_te = args.telr_dir +"/" + prefix + ".te2contig_filter.bed"
    telr_family = args.telr_dir + "/" + prefix + ".te2contig_rm.merge.bed"

    # filter genome 1 TE annotation by focus regions and families
    ## TODO: check annotation file format?
    annotation_filter = args.out + "/" + sample_name + ".annotation.filter.gff"
    filter_annotation(args.annotation1, annotation_filter, args.out, sample_name, region = args.region1, discard_family = args.discard_family)
    annotation_filter_bed = args.out + "/" + sample_name + ".annotation.filter.bed"
    liftover.gff_to_bed(annotation_filter, annotation_filter_bed)
    annotation_set = get_set(annotation_filter_bed)
    out_line = "annotated TEs in genome 1: " + str(len(annotation_set))
    write_report(stat_summary, out_line)


    # filter TELR predictions
    telr_filter = args.out + "/" + prefix + ".telr.filter.bed"
    filter_annotation(telr_meta, telr_filter, args.out, sample_name, region = args.region1, discard_family = args.discard_family)
    telr_contig_set = set()
    with open(telr_filter, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            telr_contig_set.add(entry[5])
    telr_contig_set = get_set(telr_filter, type = 2)
    out_line = "non-ref TEs predicted by TELR: " + str(len(telr_contig_set))
    write_report(stat_summary, out_line)
    
    ################### liftover and overlap with TELR predictions ###################
    out_line = "\n### stats by liftover TEs from genome1 to genome2 and compare with TELR predictions ###"
    write_report(stat_summary, out_line)

    # liftover from one genome to another
    te_lift, te_lift_meta = liftover.lift_annotation(fasta1 = args.fasta1, fasta2 = args.fasta2, bed = annotation_filter_bed, sample_name = sample_name, out_dir = args.out, preset = preset, overlap = overlap, gap = gap, flank_len = 500)

    # remove reference lifted annotations that overlap with genome2 TE annotations
    te_lift_filter = args.out + "/" + sample_name + ".lift.nonref.bed"
    with open(te_lift_filter, "w") as output:
        subprocess.call(["bedtools", "intersect", "-a", te_lift, "-b", args.annotation2, "-v"], stdout = output)
    lift_set = get_set(te_lift_filter)
    out_line = "non-ref TEs lifted from genome1 to genome2: " + str(len(lift_set))
    write_report(stat_summary, out_line)

    # overlap
    overlap = args.out + "/" + sample_name + ".overlap.bed"
    window = 3
    with open(overlap, "w") as output:
        subprocess.call(["bedtools", "window", "-w", str(window), "-a", te_lift_filter, "-b", telr_filter], stdout = output)
    share_lift_set = set()
    share_telr_contig_set = set()
    with open(overlap, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if entry[3].lower() in entry[9].lower():
                share_lift_set.add('_'.join(entry[0:4]))
                share_telr_contig_set.add(entry[11])
    out_line = "lifted non-ref TEs supported by TELR predictions: " + str(len(share_lift_set))
    write_report(stat_summary, out_line)
    lift_only_set = lift_set.difference(share_lift_set)
    out_line = "non-ref TEs predicted by TELR supported by liftover set: " + str(len(share_telr_contig_set))
    write_report(stat_summary, out_line)
    lift_only_set = lift_set.difference(share_lift_set)
    out_line = "non-ref TEs only in liftover set: " + str(len(lift_only_set))
    write_report(stat_summary, out_line)
    telr_only_contig_set = telr_contig_set.difference(share_telr_contig_set)
    out_line = "non-ref TEs only in TELR set: " + str(len(telr_only_contig_set))
    write_report(stat_summary, out_line)

    rm(overlap)

    ################## TELR flanks mapped to genome 1 ###################
    print("Lift flanking sequences from TELR contigs to genome 1...")
    out_line = "\n### stats by liftover TEs from TELR contigs to genome1 and compare with genome1 annotations ###"
    write_report(stat_summary, out_line)

    # generate bed file
    contig_te_strand_dict = dict()
    with open(telr_family, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_te_strand_dict[entry[0]] = entry[4]

    contig_te_bed = args.out + "/" + prefix + ".contig.te.filter.bed"
    with open(contig_te_bed, "w") as output, open(telr_meta, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if entry[5] in telr_contig_set:
                family = entry[3]
                contig_name = re.sub(":.*", "", entry[5])
                coord = re.sub(".*:", "", entry[5])
                start = coord.split('-')[0]
                end = coord.split('-')[1]
                strand = contig_te_strand_dict[contig_name]
                out_line = '\t'.join([contig_name, start, end, family, strand])
                output.write(out_line + '\n')
    
    sample_name = prefix
    overlap = 20
    gap = 30000
    flank_lift, flank_lift_meta = liftover.lift_annotation(fasta1 = telr_contigs, fasta2 = args.fasta1, bed = contig_te_bed, sample_name = sample_name, out_dir = args.out, preset = preset, overlap = overlap, gap = gap, flank_len = 500)

    telr_lift_contig_set = set()
    with open(flank_lift_meta, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            telr_lift_contig_set.add(entry[5])

    telr_unlift_set = telr_contig_set.difference(telr_lift_contig_set)
    out_line = "TEs can not be lifted from TELR contigs to genome1: " + str(len(telr_unlift_set))
    write_report(stat_summary, out_line)

    out_line = "TEs lifted from TELR contigs to genome1: " + str(len(telr_lift_contig_set))
    write_report(stat_summary, out_line)

    # overlap
    overlap = args.out + "/" + sample_name + ".genome1.overlap.bed"
    window = 3
    with open(overlap, "w") as output:
        subprocess.call(["bedtools", "window", "-w", str(window), "-a", flank_lift_meta, "-b", annotation_filter_bed], stdout = output)
    share_annotation_set = set()
    share_telr_lift_contig_set = set()
    with open(overlap, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            family_telr_lift_genome1 = entry[3].lower()
            family_annotation_genome1 = entry[9].lower()
            if family_annotation_genome1 in family_telr_lift_genome1:
                share_telr_lift_contig_set.add(entry[5])
                share_annotation_set.add('_'.join(entry[7:11]))
    
    out_line = "TEs lifted from TELR contigs supported by genome 1 annotation: " + str(len(share_telr_lift_contig_set))
    write_report(stat_summary, out_line)
    telr_lift_only_set = telr_lift_contig_set.difference(share_telr_lift_contig_set)
    out_line = "TEs lifted from TELR contigs not supported by genome 1 annotation: " + str(len(telr_lift_only_set))
    write_report(stat_summary, out_line)

    rm(overlap)

    #TODO compare sequence quality for all the overlapped?

    # for TELR only and unlifted flanks, map the contig to dm6 and do dnadiff comparison
    print("For TELR unlifted flanks and lifted flanks without genome1 annotation support, map contig to genome1")
    telr_check_set = telr_unlift_set.union(telr_lift_only_set)
    # out_line = "TELR TEs not lifted or not supported by genome1 annotation: " + str(len(telr_check_set))
    # write_report(stat_summary, out_line)

    # then, extract contigs and map to dm6
    telr_check_contig_list = args.out + "/" + sample_name + ".check.contig.txt"
    telr_check_contig_set = set()
    with open(telr_check_contig_list, "w") as output:
        for item in telr_check_set:
            contig = item.split(':')[0]
            telr_check_contig_set.add(contig)
            output.write(contig + '\n')
    
    telr_check_contigs = args.out + "/" + sample_name + ".check.contig.fa"
    with open(telr_check_contigs, "w") as output:
        subprocess.call(["seqtk", "subseq", telr_contigs, telr_check_contig_list], stdout = output)
    rm(telr_check_contig_list)

    telr_check_align = args.out + "/" + sample_name + ".check.contig.paf"
    with open(telr_check_align, "w") as output:
        subprocess.call(["minimap2", "-cx", preset, "-v", "0", "--secondary=no", args.fasta1, telr_check_contigs], stdout = output)
    rm(telr_check_contigs)

    # select full length mapped contig, compare using dnadiff and make mummer plot
    telr_check_align_contig_set = set()
    telr_check_dir = args.out + '/' + 'compare_telr_contig_genome1'
    mkdir(telr_check_dir)
    with open(telr_check_align, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_telr = entry[0].split('_')[0]
            if contig_telr == entry[5] and (int(entry[1]) * 0.8) < (int(entry[3]) - int(entry[2])):
                telr_check_align_contig_set.add(entry[0])
                genome1_extract_fa = telr_check_dir + "/" + '_'.join([entry[5], entry[7], entry[8]]) + '.genome1.fa'
                extract_seq(args.fasta1, entry[5], entry[7], entry[8], genome1_extract_fa)
                contig_fa = telr_check_dir + '/' + entry[0] + '.contig.fa'
                with open(contig_fa, "w") as output:
                    command = "samtools faidx " + telr_contigs + " " + entry[0]
                    subprocess.call(command, shell = True, stdout = output)
                # dnadiff
                prefix_mummer = telr_check_dir + '/' + entry[0]
                subprocess.call(["dnadiff", contig_fa, genome1_extract_fa, "-p", prefix_mummer], stderr = subprocess.DEVNULL)
                mdelta = prefix_mummer + '.mdelta'
                subprocess.call(["mummerplot", "-p", prefix_mummer, "-s", "medium", "-f", "--png", mdelta, "--color"], stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
                clean_mummer(prefix_mummer)
                rm(genome1_extract_fa)
                rm(contig_fa)
                
    telr_no_align_set = telr_check_contig_set.difference(telr_check_align_contig_set)
    rm(telr_check_align)

    # generate contig summary file
    contig_summary = args.out + "/" + "contig.summary.tsv"
    with open(contig_summary, "w") as output:
        out_line = '\t'.join(["contig_name", "genome2_support", "genome1_support"])
        output.write(out_line + '\n')
        for ins in telr_contig_set:
            if ins in share_telr_contig_set:
                genome2_support = "yes"
            else:
                genome2_support = "no"
            if ins in telr_lift_contig_set:
                if ins in share_telr_lift_contig_set:
                    genome1_support = "yes"
                else:
                    genome1_support = "no"
            else:
                genome1_support = "TELR_flank_unlifted"
            contig_name = ins.split(':')[0]
            out_line = '\t'.join([contig_name, genome2_support, genome1_support])
            output.write(out_line + '\n')
    
    # generate stats for contig summary file and write in final stats report
    out_line = "\n### stats by comparing two evaluation methods ###"
    write_report(stat_summary, out_line)
    out_line = "within " + str(len(telr_only_contig_set)) + " non-ref TEs in TELR set not in liftover TE set in genome2"
    write_report(stat_summary, out_line)
    num = len(telr_only_contig_set.intersection(share_telr_lift_contig_set))
    out_line = str(num) + " TEs can be lifted to genome1 and supported by genome1 annotations (real support insertions)"
    write_report(stat_summary, out_line)
    num = len(telr_only_contig_set.intersection(telr_lift_only_set))
    out_line = str(num) + " TEs can be lifted to genome1 but not supported by genome1 annotations (potential new insertions in dataset relative to genome1 or misannotation in genome1)"
    write_report(stat_summary, out_line)
    num = len(telr_only_contig_set.intersection(telr_unlift_set))
    out_line = str(num) + " TEs can be not be lifted to genome1 (potential false positive or new insertions in dataset relative to genome1)"
    write_report(stat_summary, out_line)
    
    print("evaluation workflow finished!")

def write_report(file, line):
    with open(file, "a") as output:
        output.write(line + '\n')

def clean_mummer(prefix):
    rm(prefix + '.mdelta')
    rm(prefix + '.mcoords')
    rm(prefix + '.gp')
    rm(prefix + '.fplot')
    rm(prefix + '.filter')
    rm(prefix + '.delta')
    rm(prefix + '.1delta')
    rm(prefix + '.1coords')
    rm(prefix + '.rplot')
    rm(prefix + '.rdiff')
    rm(prefix + '.qdiff')
    rm(prefix + '.snps')

def extract_seq(fa, chr, start, end, out):
    with open(out, "w") as output:
        command = "samtools faidx " + fa + " " + chr + ":" + start + "-" + end
        subprocess.call(command, shell = True, stdout = output)

def get_set(bed, type = 1):
    ins_set = set()
    if type == 1:
        with open(bed, "r") as input:
            for line in input:
                entry = line.replace('\n', '').split("\t")
                ins_set.add('_'.join(entry[0:4]))
    else:
        with open(bed, "r") as input:
            for line in input:
                entry = line.replace('\n', '').split("\t")
                ins_set.add(entry[5])
    return ins_set

def filter_annotation(file_in, file_out, out_dir, sample_name, region = None, discard_family = None):
    annotation = file_in
    # filter for regions
    if region is not None:
        region_filter = out_dir + "/" + sample_name + ".te.region_filter.gff"
        with open(region_filter, "w") as output:
            subprocess.call(["bedtools", "intersect", "-a", annotation, "-b", region, "-u"], stdout = output)
        annotation = region_filter
    # filter for families
    if discard_family is not None:
        family_filter = out_dir + "/" + sample_name + ".te.family_filter.gff"
        family_list = discard_family.replace(" ", "").split(',')
        with open(family_filter, "w") as output, open(annotation, "r") as input:
            for line in input:
                write = True
                for family in family_list:
                    if family in line:
                        write = False
                if write:
                    output.write(line)
        annotation = family_filter
    shutil.copyfile(annotation, file_out)
    rm(region_filter)
    rm(family_filter)
    

def filter_family(file_in, file_out, families):
    family_list = families.replace(" ", "").split(',')
    with open(file_out, "w") as output, open(file_in, "r") as input:
        for line in input:
            write = True
            for family in family_list:
                if family in line:
                    write = False
            if write:
                output.write(line)

def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        print ("Creation of the directory %s failed" % dir)

def rm(file):
    if os.path.exists(file):
        os.remove(file)

def check_lines(file):
    fname = os.path.basename(file)
    if os.path.isfile(file) and os.stat(file).st_size != 0:
        with open(file) as f:
            for i, l in enumerate(f):
                pass
            n = i + 1
            return n
    else:
        print("file doesn't exist or empty")
        return 0

main()
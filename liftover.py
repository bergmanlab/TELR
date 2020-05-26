#!/usr/bin/env python3
import sys
import argparse
import os
from Bio import SeqIO
import pandas as pd
import subprocess
import numpy as np
import re

def get_args():
    parser = argparse.ArgumentParser(description="Script to lift annotation from one genome to another")

    ## required ##
    parser.add_argument("-a", "--fasta1", type=str, help="genome 1", required=True)
    parser.add_argument("-b", "--fasta2", type=str, help="genome 2", required=True)
    parser.add_argument("-i", "--bed", type=str, help="annotation file for genome 1 in bed/gff format (default = 'bed')", required=True)
    parser.add_argument("-p", "--overlap", type=int, help="maximum overlap (default = '30000')", required=False)
    parser.add_argument("-g", "--gap", type=int, help="maximum gap (default = '25')", required=False)
    parser.add_argument("-x", "--preset", type=str, help="preset (default = 'map-pb')", required=False)
    parser.add_argument("-l", "--length", type=int, help="flanking sequence length (default = '500')", required=False)
    parser.add_argument("-o", "--out", type=str, help="output dir", required=False)

    args = parser.parse_args()

     # checks if in files exist
    try:
        test = open(args.fasta1, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.fasta2, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    if args.overlap is None:
        args.overlap = 25

    if args.gap is None:
        args.gap = 30000
    
    if args.length is None:
        args.length = 500
    
    if args.preset is None:
        args.preset = 'map-pb'
    
    if args.out is None:
        args.out = "."
    
    return args

def main():
    args = get_args()

    # if the annotation is in gff format (like from flybase), then convert to bed format first
    if "gff" in args.bed:
        ins_bed = args.out + "/ins.bed"
        gff_source = determine_gff_source(args.bed)
        gff_to_bed(args.bed, ins_bed, gff_source)
    else:
        ins_bed = args.bed
    
    sample_name = os.path.splitext(os.path.basename(args.fasta1))[0]
    
    # call liftover function
    te_report, te_report_meta = lift_annotation(args.fasta1, args.fasta2, ins_bed, sample_name, args.out, args.preset, args.overlap, args.gap, args.length)

def lift_annotation(fasta1, fasta2, bed, sample_name, out_dir = ".", preset = "map-pb", overlap = 25, gap = 30000, flank_len = 500, family_rm = None):
    # for each annotation, extract flanking sequence, map to another genome, check the distance
    flank_bed = out_dir + "/" + sample_name + ".flank.bed"
    get_flank_bed(bed, fasta1, flank_bed, flank_len)

    ## extract flanking sequences
    flank_fa = out_dir + "/" + sample_name + ".flank.fa"
    print ("Generating flanking sequences...")
    with open(flank_fa, "w") as output:
        subprocess.call(["bedtools", "getfasta", "-fi", fasta1, "-bed", flank_bed], stdout=output)
    print ("Done\n")
    # os.remove(flank_bed)

    ## map to another genome
    # minimap2 way
    mm2_out = out_dir + "/" + sample_name + ".flank.paf"
    print ("Align flanking sequence to reference...")
    with open(mm2_out, "w") as output:
        subprocess.call(["minimap2", "-cx", preset, "-v", "0", "--secondary=no", fasta2, flank_fa], stdout=output)
    print ("Done\n")
    # os.remove(flank_fa)

    # read family, strand and TE length info into dict
    te_len_dict = dict()
    contig_ins_dict = dict()
    with open(bed, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            ins_name = entry[0] + ":" + entry[1] + "-" + entry[2]
            te_len = int(entry[2]) - int(entry[1])
            te_len_dict[ins_name] = te_len
            contig_ins_dict[entry[0]] = ins_name

    if family_rm is not None:
        family_dict = dict()
        strand_dict = dict()
        with open(family_rm, "r") as input:
            for line in input:
                entry = line.replace('\n', '').split("\t")
                ins_name = contig_ins_dict[entry[0]]
                family_dict[ins_name] = entry[3]
                if entry[4] != "+" and entry[4] != "-":
                    strand_dict[ins_name] = "."
                else:
                    strand_dict[ins_name] = entry[4]
    else:
        family_dict = dict()
        strand_dict = dict()
        with open(bed, "r") as input:
            for line in input:
                entry = line.replace('\n', '').split("\t")
                ins_name = entry[0] + ":" + entry[1] + "-" + entry[2]
                family_dict[ins_name] = entry[3]
                if entry[4] != "+" and entry[4] != "-":
                    strand_dict[ins_name] = "."
                else:
                    strand_dict[ins_name] = entry[4]
    
    # read flanking info into dict
    flank_side_dict = dict() # used later for flanking sequence extraction
    flank_ins_dict = dict()
    with open(flank_bed, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            flank_name = entry[0]+":"+entry[1]+"-"+entry[2]
            flank_end = entry[3]
            flank_side_dict[flank_name] = flank_end
            flank_ins_dict[flank_name] = entry[4]
    
    # process the mapped paf
    te_remove_tmp = out_dir + "/" + sample_name + ".lift.remove.tmp.tsv"
    te_report_tmp = out_dir + "/" + sample_name + ".lift.pass.tmp.tsv"
    print ("Generating final TE annotation...")
    header = ["flank_name", "flank_len", "flank_start", "flank_end", "flank_strand", "chr", "start", "end"]
    df = pd.read_csv(mm2_out, delimiter="\t", usecols=[0,1,2,3,4,5,7,8], names=header)
    df['ins_name'] = df['flank_name'].map(flank_ins_dict)
    df['flank_side'] = df['flank_name'].map(flank_side_dict).replace('_.*', '', regex=True)
    # filter out alignments onto wrong chromosome
    if "chr" in df.iloc[0]['flank_name']:
        df['chr_flank'] = df['flank_name'].replace(':.*', '', regex=True).replace('_.*', '', regex=True)
        rm_df = df.query('chr_flank != chr')
        if len(rm_df) > 0:
            rm_df['type'] = "chr_mismatch"
            rm_df.to_csv(te_remove_tmp, sep = '\t', index=False, header=False, columns=['ins_name', 'type'])
        df.query('chr_flank == chr', inplace=True)
        df.drop(['chr_flank'], inplace=True, axis=1)
    # remove multi hit
    rm_df = df[df.duplicated(subset='flank_name', keep='first')]
    if len(rm_df) > 0:
        rm_df['type'] = "multi_hit"
        rm_df.to_csv(te_remove_tmp, sep = '\t', index=False, header=False, columns=['ins_name', 'type'], mode = 'a')
    df.drop_duplicates(subset='flank_name', keep=False, inplace = True)
    # remove if two flank map to different chr or strand
    rm_df = df.groupby('ins_name').filter(lambda x: x['chr'].nunique() > 1 or x['flank_strand'].nunique() > 1)
    if len(rm_df) > 0:
        rm_df['type'] = "multi_strand"
        rm_df.to_csv(te_remove_tmp, sep = '\t', index=False, header=False, columns=['ins_name', 'type'], mode = 'a')
    df = df.groupby('ins_name').filter(lambda x: x['chr'].nunique() == 1 and x['flank_strand'].nunique() == 1)
    # add TE family info
    df['te_family'] = df['ins_name'].map(family_dict)
    df.dropna(subset=["te_family"], inplace = True)  # remove annotation without family annotation
    df['te_strand'] = df['ins_name'].map(strand_dict)
    df['te_len'] = df['ins_name'].map(te_len_dict)
    # group and summarize by contig
    new_df = df.groupby('ins_name').apply(get_coordinate).reset_index()
    # remove if two flank have gap/overlap bigger than threshold
    new_df = new_df[((new_df['score'] >= 0) & (new_df['end'] - new_df['start'] <= gap)) | ((new_df['score'] <= 0) & (new_df['end'] - new_df['start'] <= overlap)) | (new_df['score'] == 0.5)]
    new_df['feature'] = "transposon"
    new_df.to_csv(te_report_tmp, sep = '\t', index=False, header=False, columns=['chr', 'start', 'end', 'family', 'score', 'strand', 'ins_name', 'te_strand', 'te_len'])
    print("Done\n")

    # sort gff
    te_report_tmp_sort = out_dir + "/" + sample_name + ".lift.tmp.sort.bed"
    with open(te_report_tmp_sort, "w") as output:
        command = "bedtools sort -i " + te_report_tmp
        subprocess.call(command, shell=True, stdout=output)

    # merge overlap/identical entries
    te_report_tmp_merge = out_dir + "/" + sample_name + ".lift.tmp.merge.bed"
    with open(te_report_tmp_merge, "w") as output:
        command = "bedtools merge -d 0 -o collapse -c 2,3,4,5,6,7,8,9 -delim \",\" -i " + te_report_tmp_sort
        subprocess.call(command, shell=True, stdout=output)
    
    # output overlapped/identical entries
    with open(te_report_tmp_merge, "r") as input, open(te_remove_tmp, "a") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if "," in entry[3]:
                output.write(entry[8] + '\t' + 'ins_overlap' + '\n')
    
    # find and remove duplicate entries in the remove list
    te_remove_final = out_dir + "/" + sample_name + ".lift.remove.bed"
    te_remove_set = set()
    with open(te_remove_tmp, "r") as input, open(te_remove_final, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if entry[0] not in te_remove_set:
                output.write(line)
                te_remove_set.add(entry[0])

    # final report
    te_report = out_dir + "/" + sample_name + ".final.bed"
    te_report_meta = out_dir + "/" + sample_name + ".final.meta.bed"
    with open(te_report_tmp_merge, "r") as input, open(te_report, "w") as output, open(te_report_meta, "w") as more:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            chr = entry[0]
            if "," in entry[3]:
                len_list = entry[10].split(",")
                idx = len_list.index(max(len_list))
                start = entry[3].split(",")[idx]
                end = entry[4].split(",")[idx]
                family = entry[5].split(",")[idx]
                support_type = entry[6].split(",")[idx]
                strand = entry[7].split(",")[idx]
                ins_name = entry[8].split(",")[idx]
                te_strand = entry[9].split(",")[idx]
            else:
                start = entry[3]
                end = entry[4]
                family = entry[5]
                support_type = entry[6]
                strand = entry[7]
                ins_name = entry[8]
                te_strand = entry[9]
            out_line = '\t'.join([chr, start, end, family, support_type, strand])
            output.write(out_line+"\n")
            out_line = '\t'.join([chr, start, end, family, te_strand, ins_name])
            more.write(out_line+"\n")

    # clean files
    os.remove(te_remove_tmp)
    os.remove(te_report_tmp)
    os.remove(te_report_tmp_sort)
    os.remove(te_report_tmp_merge)

    return te_report, te_report_meta

def determine_gff_source(gff):
    with open(gff, "r") as input:
        for line in input:
            if "#" not in line:
                if "RepeatMasker" in line:
                    source = "RepeatMasker"
                elif "FlyBase" in line:
                    source = "FlyBase"
                else:
                    source = "other"

def get_coordinate(df_group):
    chr = df_group['chr'].tolist()[0]
    family = df_group['te_family'].tolist()[0]
    te_strand = df_group['te_strand'].tolist()[0]
    flank_strand = df_group['flank_strand'].tolist()[0]
    te_len = df_group['te_len'].tolist()[0]
    # determine final strand
    if te_strand == ".":
        strand = "."
    elif flank_strand == te_strand:
        strand = "+"
    else:
        strand = "-"
    # determine coordinate
    if df_group.flank_name.nunique() == 1:
        score = 0.5
        if flank_strand == '+':
            if df_group['flank_side'].tolist()[0] == 'LEFT':
                start = end = df_group['end'].tolist()[0]
            else:
                start = end = df_group['start'].tolist()[0]
        else:  
            if df_group['flank_side'].tolist()[0] == 'LEFT':
                start = end = df_group['start'].tolist()[0]
            else:
                start = end = df_group['end'].tolist()[0]
    else:
        if flank_strand == '+':
            start = df_group.loc[df_group['flank_side'] == 'LEFT', 'end'].iloc[0]
            end = df_group.loc[df_group['flank_side'] == 'RIGHT', 'start'].iloc[0]
        else:
            start = df_group.loc[df_group['flank_side'] == 'LEFT', 'start'].iloc[0]
            end = df_group.loc[df_group['flank_side'] == 'RIGHT', 'end'].iloc[0]
        if start > end: # when there are overlaps
            score = end - start
            start, end = end, start
        else:
            score = end - start
    return pd.Series([chr, start, end, family, score, strand, te_strand, te_len], index=['chr', 'start', 'end', 'family', 'score', 'strand', 'te_strand', 'te_len'])


def gff_to_bed(gff, bed, source):
    with open(gff, "r") as input, open(bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if source == "FlyBase":
                family = re.sub(".*Name=|{}.*", "", entry[8])
            elif source == "RepeatMasker":
                family = re.sub(".*Motif:|\".*", "", entry[8])
            else:
                family = entry[2]
            out_line = '\t'.join([entry[0], entry[3], entry[4], family, entry[6]])
            output.write(out_line + '\n')

def get_flank_bed(gff, genome, flank_bed, flank_len=500):
    # get genome size
    genome_index = genome + ".fai"
    subprocess.call(["samtools", "faidx", genome])
    genome_len_dict = dict()
    with open(genome_index, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            genome_len_dict[entry[0]] = int(entry[1])

    # get flanking coordinates
    write_flank1=False
    write_flank2=False
    with open(gff, "r") as input, open(flank_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = entry[0]
            contig_size = genome_len_dict[contig_name]
            aln_start = int(entry[1])
            aln_end = int(entry[2])
            ins_name = contig_name + ":" + str(aln_start) + "-" + str(aln_end)
            flank1_end = aln_start - 1
            if flank1_end - flank_len > 0:
                flank1_start = flank1_end - flank_len
                write_flank1=True
            else:
                flank1_start=0
                write_flank1=False
                print("flank1 out of boundary for: " + ins_name)
            flank2_start=aln_end + 1
            if flank2_start + flank_len < contig_size - 1:
                flank2_end = flank2_start + flank_len
                write_flank2=True
            else:
                flank2_end = contig_size-1
                write_flank2=False
                print("flank2 out of boundary for: " + ins_name)
                
            if write_flank1 and write_flank2:
                output.write('\t'.join([contig_name, str(flank1_start), str(flank1_end), "LEFT_P", ins_name]) + '\n')
                output.write('\t'.join([contig_name, str(flank2_start), str(flank2_end), "RIGHT_P", ins_name]) + '\n')
            elif write_flank1:
                output.write('\t'.join([contig_name, str(flank1_start), str(flank1_end), "LEFT_S", ins_name]) + '\n')
            elif write_flank2:
                output.write('\t'.join([contig_name, str(flank2_start), str(flank2_end), "RIGHT_S", ins_name]) + '\n')

if __name__ == "__main__":                
    main()
#!/usr/bin/env python3

import sys
import argparse
import os
import os.path
import subprocess
import time
from Bio import SeqIO
import pandas as pd
import numpy as np
import re

# import utility as utl

# python3 main.py -o $output_dir -i $read_path -r $reference_path -l $te_library_path -t 16 -x pacbio

def main():
    args = get_args()
    sample_name = os.path.splitext(os.path.basename(args.read))[0]

    # Alignment
    start_time = time.time()
    bam = alignment(args.read, args.reference, args.out, sample_name, args.thread, args.presets)
    proc_time = time.time() - start_time
    print("Alignment time:", proc_time)

    # SV detection
    start_time = time.time()
    vcf = run_sniffle(bam, args.out, sample_name, args.thread)
    proc_time = time.time() - start_time
    print("SV detection time:", proc_time)

    # Sniffle output parsing
    # vcf = args.out+"/"+sample_name+".vcf"
    start_time = time.time()
    vcf_parsed, contig_reads_dir = sniffle_parse(vcf, args.out, sample_name, args.read, args.library, args.thread)
    proc_time = time.time() - start_time
    print("VCF parsing time:", proc_time)

    # local assembly on every insertion cluster
    # contig_reads_dir=args.out+"/"+"contig_reads"
    # vcf_parsed=args.out+"/"+sample_name+".vcf.parsed.filtered"
    start_time = time.time()
    contig_assembly_dir = local_assembly(contig_reads_dir, vcf_parsed, args.out, args.read, args.thread, args.presets)
    proc_time = time.time() - start_time
    print("Local assembly time:", proc_time)
    
    # annotate TEs using different method, extract flanking sequence and annotate TE family
    # contig_assembly_dir=args.out+"/"+"contig_assembly"
    # vcf_parsed=args.out+"/"+sample_name+".vcf.parsed.filtered"
    start_time = time.time()
    flank_fa, flank_bed, family_annotation, te_fa = annotate_contig(contig_assembly_dir, args.library, vcf_parsed, args.out, sample_name, args.thread, args.presets)
    proc_time = time.time() - start_time
    print("Contig annotation time:", proc_time)

    # map contig back to the reference
    # family_annotation = args.out+"/"+sample_name+".te2contig_rm.merge.bed"
    # flank_bed=args.out+"/"+sample_name+".flank.bed"
    # flank_fa=args.out+"/"+sample_name+".flank.fa"
    # te_fa=args.out+"/"+sample_name+".te.fa"
    start_time = time.time()
    find_te(flank_bed, flank_fa, args.reference, family_annotation, te_fa, args.out, sample_name, args.gap, args.overlap, args.presets)
    proc_time = time.time() - start_time
    print("Flanking alignment time:", proc_time)

def get_args():
    parser = argparse.ArgumentParser(description="Script to detect TEs in long read data")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    ## required ##
    required.add_argument("-i", "--read", type=str, help="reads in fasta/fastq format", required=True)
    required.add_argument("-r", "--reference", type=str, help="reference genome in fasta format", required=True)
    required.add_argument("-l", "--library", type=str, help="TE consensus sequences in fasta format", required=True)

    ## optional ##
    optional.add_argument("-x", "--presets", type=str, help="parameter presets for different sequencing technologies", required=False)
    optional.add_argument("-o", "--out", type=str, help="directory to output data (default = '.')", required=False)
    optional.add_argument("-t", "--thread", type=int, help="max cpu threads to use (default = '1')", required=False)
    optional.add_argument("-g", "--gap", type=int, help="max gap size for flanking sequence alignment (default = '20')", required=False)
    optional.add_argument("-p", "--overlap", type=int, help="max overlap size for flanking sequence alignment (default = '20')", required=False)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # checks if in files exist
    try:
        test = open(args.read, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.reference, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.library, 'r')
    except Exception as e:
        print(e)
        sys.exit(1)

    if args.presets is None:
        args.presets = "pacbio"

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.thread is None:
        args.thread = 1
    
    if args.gap is None:
        args.gap = 20

    if args.overlap is None:
        args.overlap = 20

    return args

def alignment(read, reference, out, sample_name, thread, presets):
    if presets == "ont":
        presets_nglmr = "ont"
    else:
        presets_nglmr = "pacbio"
    tmp_bam=out+"/"+sample_name+"_tmp.bam"
    if os.path.isfile(tmp_bam):
        print ("Raw read alignment exist")
    else:
        print ("Generating raw read alignment...")
        with open(tmp_bam,"w") as output:
            subprocess.call(["ngmlr", \
                                    "-r", reference, \
                                    "-q", read, \
                                    "-x", presets_nglmr, \
                                    "-t", str(thread), \
                                    "--rg-id", sample_name, \
                                    "--rg-sm", sample_name, \
                                    "--rg-lb", "pb", "--no-progress"], stdout=output)
        print ("Done\n")
        
    sorted_bam=out+"/"+sample_name+"_sort.bam"
    sorted_bam_idx=sorted_bam+".bai"
    if os.path.isfile(sorted_bam) and os.path.isfile(sorted_bam_idx):
        print ("Sorted indexed read alignment exist")
    else:
        print ("Generating sorted indexed read alignment...")
        command="samtools sort -@ "+str(thread)+" -o "+sorted_bam+" "+tmp_bam
        subprocess.call(command, shell=True)
        command="samtools index -@ "+str(thread)+" "+sorted_bam
        subprocess.call(command, shell=True)
        os.remove(tmp_bam)
        print ("Done\n")

    out_bam = sorted_bam
    if os.path.isfile(out_bam) == False:
            sys.stderr.write("NGMLR failed...exiting....\n")
            sys.exit(1)
    else:
        return(out_bam)

def run_sniffle(bam, out, sample_name, thread):
    vcf=out+"/"+sample_name+".vcf"
    if os.path.isfile(vcf):
        print ("Sniffle output exist")
    else:
        print ("Generating sniffle output...")
        command="sniffles -n -1 --genotype --report_seq --report_BND --threads "+str(thread)+" -m "+bam+" -v "+vcf
        subprocess.call(command, shell=True)
        print ("Done\n")
    if os.path.isfile(vcf) == False:
            sys.stderr.write("Sniffles failed...exiting....\n")
            sys.exit(1)
    else:
        return(vcf)

def sniffle_parse(vcf, out, sample_name, raw_reads, TE_library, thread):
    contig_reads_dir=out+"/"+"contig_reads"
    mkdir(contig_reads_dir)

    vcf_parsed=out+"/"+sample_name+".vcf.parsed"
    print ("Generating parsed sniffle VCF...")
    query_str="\"%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\n\""
    command="bcftools query -i \'SVTYPE=\"INS\" & ALT!=\"<INS>\"\' -f "+query_str+" "+vcf
    with open(vcf_parsed, "w") as output:
        subprocess.call(command, stdout=output, shell=True)
    # remove redundancy in parsed vcf
    remove_vcf_redundancy(vcf_parsed)
    print ("Done\n")
    
    # constrct fasta from parsed vcf file
    ins_seqs=out+"/"+sample_name+".ins.fasta"
    print ("Generating VCF sequences...")
    vcf2fasta(vcf_parsed, ins_seqs)
    print ("Done\n")

    # run RM on the inserted seqeunce
    ins_rm_out=ins_seqs+".out.gff"
    print ("Generating VCF sequences repeatmasker output ...")
    subprocess.call(["RepeatMasker", "-dir", out, "-gff", "-s", "-nolow", "-no_is", "-xsmall", "-e", "ncbi", "-lib", TE_library, "-pa", str(thread), ins_seqs])
    print ("Done\n")

    # extract VCF sequences that contain TEs 
    with open(ins_rm_out, "r") as input:
        te_seqs_dict = dict()
        for line in input:
            if "RepeatMasker" in line:
                entry = line.replace('\n', '').split("\t")
                seq_len = int(entry[4]) - int(entry[3])
                if entry[0] in te_seqs_dict:
                    te_seqs_dict[entry[0]] = seq_len + te_seqs_dict[entry[0]]
                else:
                    te_seqs_dict[entry[0]] = seq_len

    # filter parsed VCF using TE-sequence list and 
    vcf_parsed_filtered=out+"/"+sample_name+".vcf.parsed.filtered"
    filter_vcf(vcf_parsed, te_seqs_dict, vcf_parsed_filtered)

    # constrct fasta from parsed filtered vcf file
    ins_te_seqs=out+"/"+sample_name+".ins.te.fasta"
    vcf2fasta(vcf_parsed_filtered, ins_te_seqs)

    # extract read IDs
    ids = out+"/"+sample_name+".ids"
    with open(vcf_parsed_filtered, "r") as input, open(ids, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            read_list = entry[8].split(",")
            for read in read_list:
                output.write(read+"\n")
    
    # generate unique ID list
    ids_unique = ids+".unique"
    command = "cat "+ids+" |sort|uniq"
    with open(ids_unique, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # filter raw reads using read list
    print ("Filter raw reads around TE insertions...")
    sub_fa=out+"/"+sample_name+".subreads.fa"
    start_time = time.time()
    command = "seqtk subseq "+raw_reads+" "+ids_unique+" | seqtk seq -a"
    with open(sub_fa, "w") as output:
        subprocess.call(command, stdout=output, shell=True)
    proc_time = time.time() - start_time
    print("filter raw read time:", proc_time)
    print ("Done\n")

    # reorder
    print ("Reorder reads...")
    sub_fa_reorder=out+"/"+sample_name+".subreads.reorder.fa"
    start_time = time.time()
    extract_reads(sub_fa, ids, sub_fa_reorder)
    proc_time = time.time() - start_time
    print("reorder filtered reads time:", proc_time)
    print ("Done\n")

    # separate reads into multiple files, using csplit
    print ("Split reads by cluster...")
    os.chdir(contig_reads_dir)
    print(os.getcwd())
    m = []
    k = 1
    with open(vcf_parsed_filtered,"r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            read_list = entry[8].split(",")
            k = k + 2*(len(read_list))
            m.append(k)
    m = m[:-1]
    index=" ".join(str(i) for i in m)
    command="csplit -s -f contig -n 1 "+sub_fa_reorder+" "+index
    subprocess.call(command, shell=True)
    print ("Done\n")

    os.remove(sub_fa)
    os.remove(sub_fa_reorder)

    return vcf_parsed_filtered, contig_reads_dir

def local_assembly(contig_reads_dir, vcf, out, raw_reads, thread, presets):
    contig_assembly_dir=out+"/"+"contig_assembly"
    mkdir(contig_assembly_dir)

    if presets == "ont":
        presets_wtdbg2 = "ont"
    else:
        presets_wtdbg2 = "rs"

    print ("Assemble contigs...")
    k = 0
    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            wtdbg2_out = contig_assembly_dir + "/" + contig_name + ".fa"
            contig_reads=contig_reads_dir + "/contig" + str(k)
            k=k+1
            run_wtdbg2(contig_reads, wtdbg2_out, thread, presets_wtdbg2)
    print ("Done\n")
    return(contig_assembly_dir)

def run_wtdbg2(fastq, wtdbg2_out, thread, presets):
    prefix = os.path.splitext(fastq)[0]
    command = "wtdbg2 -x "+presets+" -q -g 30k -t "+str(thread)+" -i "+fastq+" -fo "+prefix
    subprocess.call(command, shell = True)
    # subprocess.call(["wtdbg2", "-x", "rs", "-q", "-g", "30k", "-t", str(thread), "-i", fastq, "-fo", prefix])
    contig_layout=prefix+".ctg.lay.gz"
    if os.path.isfile(contig_layout) and os.stat(contig_layout).st_size != 0:
        command = "wtpoa-cns -t "+"4"+" -i "+contig_layout+" -fo "+wtdbg2_out
        subprocess.call(command, shell = True)
    else:
        print("assembly failed for "+fastq+"\n")
    # subprocess.call(["wtpoa-cns", "-q", "-t", str(thread), "-i", contig_layout, "-fo", wtdbg2_out])

    # if os.path.isfile(wtdbg2_out) == False:
    #     sys.stderr.write("wtdbg2 failed for "+prefix+" , exiting..."+"\n")
    #     sys.exit(1)
    
def annotate_contig(asm_dir, TE_library, vcf, out, sample_name, thread, presets):
    if presets == "ont":
        presets_minimap2 = "map-ont"
    else:
        presets_minimap2 = "map-pb"

    # merge all contigs into a single file
    merge_contigs=out+"/"+sample_name+".contigs.fa"
    contig_list=out+"/"+sample_name+".contigs.list"
    # print ("Generate merged contig file...")
    with open(merge_contigs, "w") as output_contigs, open(contig_list, "w") as output_list:
        for file in os.listdir(asm_dir):
            if ".fa" in file and os.stat(asm_dir+"/"+file).st_size > 0:
                contig_name=os.path.splitext(file)[0]
                with open(asm_dir+"/"+file, "rU") as handle:
                    records = SeqIO.parse(handle, "fasta")
                    for record in records:
                        if record.id == 'ctg1':
                            record.id = contig_name
                            record.description = "len="+str(len(record.seq))
                            SeqIO.write(record, output_contigs, 'fasta')
                            output_list.write(contig_name+"\n")

    # map sequence to contigs
    print ("Map VCF sequences to contigs...")
    seq2contig_out=out+"/"+"seq2contig.paf"
    if os.path.isfile(seq2contig_out):
        os.remove(seq2contig_out)

    # constrct fasta from parsed vcf file
    ins_seqs=out+"/"+sample_name+".ins.filter.fasta"
    print ("Generating filtered VCF sequences...")
    vcf2fasta(vcf, ins_seqs)
    print ("Done\n")

    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            vcf_seq=entry[7]
            sv_len=entry[3]
            query=out+"/"+contig_name+".seq.fa"
            create_fa(contig_name, vcf_seq, query)
            subject=asm_dir+"/"+contig_name+".fa"
            if os.path.isfile(subject):
                with open(seq2contig_out, "a") as output:
                    subprocess.call(["minimap2", "-cx", presets_minimap2, "--secondary=no", subject, query], stdout=output)
            # if sv_len=="999999999":
            #     seq2contig(query, subject, seq2contig_out)
            # else:
            #     seq2contig(query, subject, seq2contig_out)
            os.remove(query)
    seq2contig_bed=out+"/"+"seq2contig.bed"
    ## covert to bed format
    with open(seq2contig_out, "r") as input, open(seq2contig_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            seq_id = entry[0]
            contig_id = entry[5]
            bed_line = "\t".join([entry[0],entry[7],entry[8],entry[5],entry[11],entry[4]])
            output.write(bed_line+"\n")
    print ("Done\n")

    # map TE library to contigs using minimap2
    ## TE-contig alignment
    te2contig_out=out+"/"+sample_name+".te2contig.paf"
    print ("Generating TE-contig alignment...")
    if os.path.isfile(te2contig_out):
        os.remove(te2contig_out)
    with open(contig_list, "r") as input:
        for line in input:
            contig_name = line.replace('\n', '')
            contig = out+"/"+contig_name+".fa"
            with open(contig, "w") as output:
                    subprocess.call(["samtools", "faidx", merge_contigs, contig_name], stdout=output)
            # map TE library to contig using minimap2 map-pb -p 0.8 -c
            with open(te2contig_out, "a") as output:
                    subprocess.call(["minimap2", "-cx", presets_minimap2, contig, TE_library, "-t", str(thread)], stdout=output, stderr=subprocess.DEVNULL)
            # remove contig file
            os.remove(contig)
    # convert to bed format
    te2contig_bed=out+"/"+sample_name+".te2contig.bed"
    with open(te2contig_out, "r") as input, open(te2contig_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            bed_line = "\t".join([entry[5],entry[7],entry[8],entry[0],entry[11],entry[4]])
            output.write(bed_line+"\n")
    print ("Done\n")

    # Use VCF sequence alignment to filter minimap2 TE-contig alignment
    te2contig_filter_raw=out+"/"+sample_name+".te2contig_filter.tsv"
    command="bedtools intersect -a "+te2contig_bed+" -b "+seq2contig_bed+" -wao"
    print ("Filter minTE-contigimap2 alignment...")
    with open(te2contig_filter_raw, "w") as output:
        subprocess.call(command, shell = True, stdout=output)
    print ("Done")

    # filter and merge
    ## get rid of -1 and make it into bed format
    te2contig_filter_tmp_bed=out+"/"+sample_name+".te2contig_filter.tmp.bed"
    with open(te2contig_filter_raw, "r") as input, open(te2contig_filter_tmp_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if int(entry[12]) > 10: # the overlap between VCF sequence alignment and TE-contig alingment has to be over 10bp
                out_line = "\t".join([entry[0], entry[1], entry[2], entry[3], entry[4], entry[5]])
                output.write(out_line+"\n")
    ## sort
    te2contig_filter_tmp_sort_bed=out+"/"+sample_name+".te2contig_filter.tmp.sort.bed"
    command="bedtools sort -i "+te2contig_filter_tmp_bed
    with open(te2contig_filter_tmp_sort_bed, "w") as output:
        subprocess.call(command, shell = True, stdout=output)

    ## merge
    te2contig_filter_bed=out+"/"+sample_name+".te2contig_filter.bed"
    command="bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim \"|\" -i "+te2contig_filter_tmp_sort_bed
    with open(te2contig_filter_bed, "w") as output:
        subprocess.call(command, shell = True, stdout=output)

    ## remove tmp files
    os.remove(te2contig_filter_raw)
    os.remove(te2contig_filter_tmp_bed)
    os.remove(te2contig_filter_tmp_sort_bed)

    # get contig length
    contig_len_dict = dict() # used later for flanking sequence extraction
    with open(merge_contigs, "r") as output_contigs:
        for record in SeqIO.parse(merge_contigs, "fasta"):
            contig_len_dict[record.id] = len(record.seq)

    ## extract flanking, need bed and contig length dict
    flank_bed=out+"/"+sample_name+".flank.bed"
    if os.path.isfile(flank_bed):
        os.remove(flank_bed)
    get_flank_seqs(contig_len_dict, te2contig_filter_bed, merge_contigs, flank_bed, flank_len=500)

    ## extract flanking sequences
    flank_fa=out+"/"+sample_name+".flank.fa"
    print ("Generating flanking sequences...")
    with open(flank_fa, "w") as output:
        subprocess.call(["bedtools", "getfasta", "-fi", merge_contigs, "-bed", flank_bed], stdout=output)
    print ("Done\n")

    ## extract sequence and RM
    te_fa=out+"/"+sample_name+".te.fa"
    with open(te_fa, "w") as output:
        subprocess.call(["bedtools", "getfasta", "-fi", merge_contigs, "-bed", te2contig_filter_bed], stdout=output)
    te_rm_out=te_fa+".out.gff"
    print ("Generating TE sequences repeatmasking output...")
    subprocess.call(["RepeatMasker", "-dir", out, "-gff", "-s", "-nolow", "-no_is", "-xsmall", "-e", "ncbi", "-lib", TE_library, "-pa", str(thread), te_fa])
    print ("Done\n")

    ## parse and merge
    te2contig_rm=out+"/"+sample_name+".te2contig_rm.bed"
    print ("Repeatmask TE sequences repeatmasker...")
    with open(te_rm_out, "r") as input, open(te2contig_rm, "w") as output:
        for line in input:
            if "##" not in line:
                entry = line.replace('\n', '').split("\t")
                contigs = entry[0].replace(':', '-').split("-")
                family = re.sub('Target \"Motif:|\".*', '', entry[8])
                strand = entry[6]
                score = entry[5]
                out_line="\t".join([contigs[0], contigs[1], contigs[2], family, score, strand])
                output.write(out_line+"\n")
    print("Done\n")

    te2contig_rm_merge=out+"/"+sample_name+".te2contig_rm.merge.bed"
    command="bedtools merge -c 4,6 -o distinct -delim \"|\" -i "+te2contig_rm
    with open(te2contig_rm_merge, "w") as output:
        subprocess.call(command, shell = True, stdout=output)
    
    return flank_fa, flank_bed, te2contig_rm_merge, te_fa
    
def find_te(flank_bed, flank_seq, ref, family_annotation, te_fa, out, sample_name, gap, overlap, presets):
    if presets == "ont":
        presets_minimap2 = "map-ont"
    else:
        presets_minimap2 = "map-pb"

    # minimap2 way
    mm2_out=out+"/"+sample_name+".mm2.paf"
    print ("Align flanking sequence to reference...")
    with open(mm2_out, "w") as output:
        subprocess.call(["minimap2", "-cx", presets_minimap2, ref, flank_seq], stdout=output)
    print ("Done\n")

    # read family and strand annotation into dict
    family_dict = dict()
    strand_dict = dict()
    with open(family_annotation, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            family_dict[entry[0]] = entry[3]
            if entry[4] != "+" and entry[4] != "-":
                strand_dict[entry[0]] = "."
            else:
                strand_dict[entry[0]] = entry[4]

    # read TE sequence length into dict
    te_len_dict = dict()
    with open(te_fa, "r") as input:
        for record in SeqIO.parse(input, "fasta"):
            contig_name = re.sub('>|:.*', '', record.id)
            te_len_dict[contig_name] = len(record.seq)
    
    # read flanking info into dict
    flank_dict = dict() # used later for flanking sequence extraction
    with open(flank_bed, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            flank_name = entry[0]+":"+entry[1]+"-"+entry[2]
            flank_end = entry[3]
            flank_dict[flank_name] = flank_end

    # prase minimap2 output and report
    te_report_tmp = out+"/"+sample_name+".tmp.gff"
    print ("Generating final TE annotation...")
    header = ["flank_name", "flank_len", "flank_start", "flank_end", "flank_strand", "chr", "start", "end"]
    df = pd.read_csv(mm2_out, delimiter="\t", usecols=[0,1,2,3,4,5,7,8], names=header)
    df['contig_name'] = df.flank_name.str.replace(':.*', '', regex=True)
    df['flank_end'] = df['flank_name'].map(flank_dict).replace('_.*', '', regex=True)
    # remove multi hit
    df.drop_duplicates(subset='flank_name', keep=False, inplace = True)
    # remove if two flank map to different chr or strand
    df = df.groupby('contig_name').filter(lambda x: x['chr'].nunique() == 1 and x['flank_strand'].nunique() == 1)
    # add TE family info
    df['te_family'] = df['contig_name'].map(family_dict)
    # add TE strand info
    df['te_strand'] = df['contig_name'].map(strand_dict)
    # add TE sequence length info
    df['te_len'] = df['contig_name'].map(te_len_dict)
    # group and summarize by contig
    new_df = df.groupby('contig_name').apply(get_coordinate).reset_index()
    # remove if two flank have gap/overlap bigger than threshold
    new_df = new_df[((new_df['score'] == 2) & (new_df['end'] - new_df['start'] <= gap)) | ((new_df['score'] == 3) & (new_df['end'] - new_df['start'] <= overlap)) | (new_df['score'] == 1)]
    new_df['feature'] = "transposon"
    # new_df['phase'] = "."
    # merge entries with duplicate coordinates
    # new_df = new_df.groupby(['chr', 'start', 'end']).apply(merge_prediction).reset_index()
    # output
    new_df.to_csv(te_report_tmp, sep = '\t', index=False, header=False, columns=['chr', 'contig_name', 'family', 'start', 'end', 'score', 'strand', 'te_len', 'te_strand'])
    print("Done\n")

    # sort gff
    te_report_tmp_sort = out+"/"+sample_name+".tmp.sort.gff"
    with open(te_report_tmp_sort, "w") as output:
        command = "bedtools sort -i " + te_report_tmp
        subprocess.call(command, shell=True, stdout=output)

    # merge overlap/identical entries
    te_report_tmp_merge = out+"/"+sample_name+".tmp.merge.gff"
    with open(te_report_tmp_merge, "w") as output:
        command = "bedtools merge -d 0 -o collapse -c 2,3,4,5,6,7,8,9 -delim \",\" -i " + te_report_tmp_sort
        subprocess.call(command, shell=True, stdout=output)

    # output overlapped/identical entries
    te_overlap_final = out+"/"+sample_name+".final.overlap"
    with open(te_report_tmp_merge, "r") as input, open(te_overlap_final, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            if "," in entry[3]:
                output.write(line+"\n")
    
    te_report_final = out+"/"+sample_name+".final.bed"
    contig_dict = dict()
    with open(te_report_tmp_merge, "r") as input, open(te_report_final, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            chr = entry[0]
            if "," in entry[3]:
                start = min(entry[5].split(','))
                end = max(entry[6].split(','))

                len_list = entry[9].split(",")
                idx = len_list.index(min(len_list))
                contig_name = entry[3].split(",")[idx]
                family = entry[4].split(",")[idx]
                support_type = entry[7].split(",")[idx]
                strand = entry[8].split(",")[idx]
                te_strand = entry[10].split(",")[idx]
            else:
                start = entry[5]
                end = entry[6]
                contig_name = entry[3]
                family = entry[4]
                support_type = entry[7]
                strand = entry[8]
                te_strand = entry[10]
            contig_dict[contig_name] = [chr, start, end, family, te_strand]
            out_line = '\t'.join([chr, start, end, family, support_type, strand])
            output.write(out_line+"\n")

    # generate TE sequence fasta
    final_te_seqs = out+"/"+sample_name+".final.fa"
    if os.path.isfile(final_te_seqs):
        os.remove(final_te_seqs)
    with open(te_fa, "r") as input, open(final_te_seqs, "a") as output:
        for record in SeqIO.parse(input, "fasta"):
            contig_name = (record.id).split(":")[0]
            if contig_name in contig_dict:
                chr = contig_dict[contig_name][0]
                start = contig_dict[contig_name][1]
                end = contig_dict[contig_name][2]
                family = contig_dict[contig_name][3]
                te_strand = contig_dict[contig_name][4]
                record.id = chr+"_"+str(start)+"_"+str(end)+"#"+family
                if te_strand == "+" or te_strand == ".":
                    output.write(">"+record.id+"\n"+str(record.seq)+"\n")
                else:
                    output.write(">"+record.id+"\n"+str(record.seq.reverse_complement())+"\n")

def get_coordinate(df_group):
    chr = df_group['chr'][0]
    family = df_group['te_family'][0]
    te_strand = df_group['te_strand'][0]
    flank_strand = df_group['flank_strand'][0]
    te_len = df_group['te_len'][0]
    # determine final strand
    if te_strand == ".":
        strand = "."
    elif flank_strand == te_strand:
        strand = "+"
    else:
        strand = "-"
    # determine coordinate
    if df_group.flank_name.nunique() == 1:
        score = 1
        if flank_strand == '+':
            if df_group['flank_end'][0] == 'LEFT':
                start = end = df_group['end'][0]
            else:
                start = end = df_group['start'][0]
        else:
            if df_group['flank_end'][0] == 'LEFT':
                start = end = df_group['start'][0]
            else:
                start = end = df_group['end'][0]
    else:
        if flank_strand == '+':
            start = df_group.loc[df_group['flank_end'] == 'LEFT', 'end'].iloc[0]
            end = df_group.loc[df_group['flank_end'] == 'RIGHT', 'start'].iloc[0]
        else:
            start = df_group.loc[df_group['flank_end'] == 'LEFT', 'start'].iloc[0]
            end = df_group.loc[df_group['flank_end'] == 'RIGHT', 'end'].iloc[0]
        if start > end: # when there are gaps
            start, end = end, start
            score = 2
        else:
            score = 3
    return pd.Series([chr, start, end, family, score, strand, te_strand, te_len], index=['chr', 'start', 'end', 'family', 'score', 'strand', 'te_strand', 'te_len'])

def create_fa(header, seq, out):
    with open(out, "w") as output:
        output.write(">"+header+"\n")
        output.write(seq)

def seq2contig(seq, contig, out):
        with open(out, "a") as output:
            subprocess.call(["minimap2", "-cx", "map-pb", "--secondary=no", contig, seq], stdout=output) # only retain primary alignment

def extract_reads(reads, list, out):
    record_dict = SeqIO.index(reads, "fasta")
    with open(out, "wb") as output_handle, open(list, "r") as ID:
        for entry in ID:
            entry = entry.replace('\n', '')
            # sequences=record_dict[entry]
            output_handle.write(record_dict.get_raw(entry))
            # SeqIO.write(sequences, output_handle, "fasta")

def vcf2fasta(vcf, out):
    with open(vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            coord='_'.join([entry[0],entry[1],entry[2]])
            output.write(">"+coord+"\n")
            output.write(entry[7]+"\n")

def unique(list1): 
    # insert the list to the set 
    list_set = set(list1) 
    # convert the set to the list 
    unique_list = (list(list_set)) 
    return(unique_list)

def af_sum(nums):
    sum_nums=sum(nums)
    if sum_nums>1:
        sum_nums=int(sum_nums)
    return(sum_nums)

def id_merge(strings):
    string_merged=','.join(strings)
    ids=string_merged.split(",")
    id_string=','.join(unique(ids))
    return(id_string)

def remove_vcf_redundancy(vcf):
    header = ["chr", "start", "end", "length", "coverage", "AF", "ID", "seq", "ids", "filter"]
    df = pd.read_csv(vcf, delimiter="\t", names=header)
    df2 = df.groupby(['chr', 'start', 'end']).agg({'length': 'first', 'coverage': 'sum', 'AF': af_sum, 'ID': 'first', 'seq': 'first', 'ids': id_merge, 'filter': 'first'}).reset_index()
    df2.to_csv(vcf, sep='\t', header=False, index=False)

def filter_vcf(vcf, dict, out, cov=0):
    with open(vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")  
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            length = entry[3]
            if length=="999999999":
                length = 999
            else:
                length = int(length)
            if contig_name in dict:
                if dict[contig_name] >= cov * length: # no VCF sequence TE coverage filtering for now
                    output.write(line)

def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        print ("Creation of the directory %s failed" % dir)
    else:
        print ("Successfully created the directory %s " % dir)

def get_flank_seqs(contig_len, bed, contigs, flank_bed, flank_len=500):
    # get flanking coordinates
    write_flank1=False
    write_flank2=False
    with open(bed, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = entry[0]
            contig_size = contig_len[contig_name]
            aln_start = int(entry[1])
            aln_end = int(entry[2])
            flank1_end = aln_start - 1
            if flank1_end - flank_len > 0:
                flank1_start = flank1_end - flank_len
                write_flank1=True
            else:
                flank1_start=0
                write_flank1=False
            flank2_start=aln_end + 1
            if flank2_start + flank_len < contig_size - 1:
                flank2_end = flank2_start + flank_len
                write_flank2=True
            else:
                flank2_end = contig_size-1
                write_flank2=False
                
            if write_flank1 and write_flank2:
                with open(flank_bed, "a") as output:
                    output.write(contig_name+"\t"+str(flank1_start)+"\t"+str(flank1_end)+"\tLEFT_P\n")
                with open(flank_bed, "a") as output:
                    output.write(contig_name+"\t"+str(flank2_start)+"\t"+str(flank2_end)+"\tRIGHT_P\n")
            elif write_flank1:
                with open(flank_bed, "a") as output:
                    output.write(contig_name+"\t"+str(flank1_start)+"\t"+str(flank1_end)+"\tLEFT_S\n")
            elif write_flank2:
                with open(flank_bed, "a") as output:
                    output.write(contig_name+"\t"+str(flank2_start)+"\t"+str(flank2_end)+"\tRIGHT_S\n")
                

main()
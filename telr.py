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
from multiprocessing import Process, Pool
import liftover

# import utility as utl

# python3 telr.py -o $output_dir -i $read_path -r $reference_path -l $te_library_path -t 16 -x pacbio

def main():
    args = get_args()
    print("CMD: ", ' '.join(sys.argv), "\n")
    sample_name = os.path.splitext(os.path.basename(args.read))[0]

    # Alignment
    start_time = time.time()
    bam = alignment(args.read, args.reference, args.out, sample_name, args.thread, args.presets)
    proc_time = time.time() - start_time
    print("Alignment time:", proc_time, "\n")

    # SV detection
    start_time = time.time()
    vcf = run_sniffle(bam, args.reference, args.out, sample_name, args.thread)
    proc_time = time.time() - start_time
    print("SV detection time:", proc_time, "\n")

    # Sniffle output parsing
    # vcf = args.out+"/"+sample_name+".vcf"
    start_time = time.time()
    vcf_parsed, contig_reads_dir = sniffle_parse(vcf, args.out, sample_name, args.read, args.library, args.thread)
    proc_time = time.time() - start_time
    print("SV parsing time:", proc_time, "\n")

    # local assembly on every insertion cluster
    # contig_reads_dir=args.out+"/"+"contig_reads"
    # vcf_parsed=args.out+"/"+sample_name+".vcf.parsed.filtered"
    start_time = time.time()
    contig_assembly_dir = local_assembly(contig_reads_dir, vcf_parsed, args.out, args.read, args.thread, args.presets, args.polish)
    proc_time = time.time() - start_time
    print("Local assembly time:", proc_time, "\n")
    
    # annotate TEs using different method, extract flanking sequence and annotate TE family
    # contig_assembly_dir=args.out+"/"+"contig_assembly"
    # vcf_parsed=args.out+"/"+sample_name+".vcf.parsed.filtered"
    start_time = time.time()
    te_contigs_annotation, family_annotation, te_fa, merge_contigs = annotate_contig(contig_assembly_dir, args.library, vcf_parsed, args.out, sample_name, args.thread, args.presets)
    proc_time = time.time() - start_time
    print("Contig annotation time:", proc_time, "\n")

    # map contig back to the reference
    # te_contigs_annotation = args.out+"/"+sample_name+".te2contig_filter.bed"
    # family_annotation = args.out+"/"+sample_name+".te2contig_rm.merge.bed"
    # te_fa=args.out+"/"+sample_name+".te.fa"
    # merge_contigs=args.out+"/"+sample_name+".contigs.fa"
    start_time = time.time()
    find_te(merge_contigs, args.reference, te_contigs_annotation, family_annotation, te_fa, args.out, sample_name, args.gap, args.overlap, args.presets)
    proc_time = time.time() - start_time
    print("Flanking alignment time:", proc_time, "\n")

    print("TELR finished!")

def get_args():
    parser = argparse.ArgumentParser(description="Script to detect TEs in long read data")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    ## required ##
    required.add_argument("-i", "--read", type=str, help="reads in fasta/fastq format", required=True)
    required.add_argument("-r", "--reference", type=str, help="reference genome in fasta format", required=True)
    required.add_argument("-l", "--library", type=str, help="TE consensus sequences in fasta format", required=True)

    ## optional ##
    optional.add_argument("-x", "--presets", type=str, help="parameter presets for different sequencing technologies (default = 'pacbio')", required=False)
    optional.add_argument("-p", "--polish", type=int, help="rounds of contig polishing (default = 1)", required=False)
    optional.add_argument("-o", "--out", type=str, help="directory to output data (default = '.')", required=False)
    optional.add_argument("-t", "--thread", type=int, help="max cpu threads to use (default = '1')", required=False)
    optional.add_argument("-g", "--gap", type=int, help="max gap size for flanking sequence alignment (default = '20')", required=False)
    optional.add_argument("-v", "--overlap", type=int, help="max overlap size for flanking sequence alignment (default = '20')", required=False)
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

    if args.polish is None:
        args.polish = 1
    
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
        # print ("Done\n")
        
    sorted_bam=out+"/"+sample_name+"_sort.bam"
    sorted_bam_idx=sorted_bam+".bai"
    print ("Sort and index read alignment...")
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

def run_sniffle(bam, reference, out, sample_name, thread, svim = False):
    print ("Generating SV output...")
    if svim:
        subprocess.call(["svim", "alignment", "--insertion_sequences", "--read_names", "--sample", sample_name, "--duplications_as_insertions", out, bam, reference])
        vcf = out+"/"+"variants.vcf"
    else:
        vcf = out+"/"+sample_name+".vcf"
        command="sniffles -n -1 --genotype --report_seq --report_BND --threads "+str(thread)+" -m "+bam+" -v "+vcf
        subprocess.call(command, shell=True)
    print ("Done\n")
    if os.path.isfile(vcf) == False:
            sys.stderr.write("SV detection failed...exiting....\n")
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
    #TODO check whether vcf file contains insertions, quit if 0
    # remove redundancy in parsed vcf
    vcf_unique=out+"/"+sample_name+".vcf.reduced"
    remove_vcf_redundancy(vcf_parsed, vcf_unique)
    os.remove(vcf_parsed)
    vcf_parsed=vcf_unique
    print ("Done\n")
    
    # constrct fasta from parsed vcf file
    ins_seqs=out+"/"+sample_name+".ins.fasta"
    vcf2fasta(vcf_parsed, ins_seqs)

    # run RM on the inserted seqeunce
    print ("Repeatmask VCF sequences...")
    subprocess.call(["RepeatMasker", "-dir", out, "-gff", "-s", "-nolow", "-no_is", "-xsmall", "-e", "ncbi", "-lib", TE_library, "-pa", str(thread), ins_seqs])
    print ("Done\n")
    ins_rm_out=ins_seqs+".out.gff"

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
    #TODO check how many TEs in filtered set, quit if 0
    os.remove(vcf_parsed)

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
    sub_fa=out+"/"+sample_name+".subreads.fa"
    command = "seqtk subseq "+raw_reads+" "+ids_unique+" | seqtk seq -a"
    with open(sub_fa, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # reorder reads
    sub_fa_reorder=out+"/"+sample_name+".subreads.reorder.fa"
    extract_reads(sub_fa, ids, sub_fa_reorder)

    # separate reads into multiple files, using csplit
    work_dir=os.getcwd()
    # print ("current working dir:"+work_dir)
    os.chdir(contig_reads_dir)
    # print ("New working dir:"+os.getcwd())
    m = []
    k = 1
    with open(vcf_parsed_filtered,"r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            read_list = entry[8].split(",")
            k = k + 2*(len(read_list))
            m.append(k)
    if len(m) == 1:
        subprocess.call(["cp", sub_fa_reorder, "contig0"])
    elif len(m) == 0:
        print ("No insertion detected, exiting...")
    else:
        m = m[:-1]
        index=" ".join(str(i) for i in m)
        command="csplit -s -f contig -n 1 "+sub_fa_reorder+" "+index
        subprocess.call(command, shell=True)
        # print ("Done\n")

    os.remove(ids)
    os.remove(ids_unique)
    os.remove(sub_fa)
    os.remove(sub_fa_reorder)
    os.chdir(work_dir)
    # print ("New working dir:"+os.getcwd())

    return vcf_parsed_filtered, contig_reads_dir

def local_assembly(contig_reads_dir, vcf, out, raw_reads, thread, presets, polish):
    contig_assembly_dir=out+"/"+"contig_assembly"
    mkdir(contig_assembly_dir)

    if presets == "ont":
        presets_wtdbg2 = "ont"
        presets_minimap2 = "map-ont"
    else:
        presets_wtdbg2 = "rs"
        presets_minimap2 = "map-pb"

    print ("Assemble contigs...")
    k = 0
    asm_pa_list=[]
    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            contig_reads=contig_reads_dir + "/contig" + str(k)
            # rename contig reads
            contig_reads_rename=contig_reads_dir + "/" + contig_name + ".reads.fa"
            os.rename(contig_reads, contig_reads_rename)
            thread_asm=1
            asm_pa = [contig_reads_rename, contig_assembly_dir, contig_name, thread_asm, presets_wtdbg2, presets_minimap2, polish]
            asm_pa_list.append(asm_pa)
            k = k + 1
    # run assembly in parallel
    pool = Pool(processes=thread)
    pool.map(run_wtdbg2, asm_pa_list)
    pool.close()
    pool.join()
    print ("Done\n")
    return(contig_assembly_dir)

def run_wtdbg2(args):
    reads = args[0]
    asm_dir = args[1]
    contig_name = args[2]
    thread = args[3]
    presets_wtdbg2 = args[4]
    presets_minimap2 = args[5]
    polish = args[6]

    prefix = reads.replace('.reads.fa', '')
    command = "wtdbg2 -x "+presets_wtdbg2+" -q -AS 1 -g 30k -t "+str(thread)+" -i "+reads+" -fo "+prefix
    try:
        subprocess.run(command, shell = True, timeout = 300)
    except subprocess.TimeoutExpired:
        print("fail to build contig layout for contig: " + contig_name)
        return

    contig_layout=prefix+".ctg.lay.gz"
    if check_exist(contig_layout):
        # derive consensus
        cns_thread = str(min(thread, 4))
        contig_raw=prefix+".ctg.lay.gz"
        consensus = prefix + ".raw.fa"
        command = "wtpoa-cns -q -t "+cns_thread+" -i "+contig_layout+" -fo "+consensus
        try:
            subprocess.run(command, shell = True, timeout = 300)
        except subprocess.TimeoutExpired:
            print("fail to assemble contig: " + contig_name)
            return
        if check_exist(consensus):
            consensus_final = asm_dir + "/" + contig_name + ".cns.fa"
            if polish > 0:
                polish_contig(prefix, reads, consensus, consensus_final, thread, presets_minimap2, polish, contig_name)
            else:
                # move raw consensus to asm dir
                os.rename(consensus, consensus_final)
        else:
            print("Initial assembly fail for "+prefix+"\n")
    else:
        print("Build contig layout fail for "+prefix+"\n")

def check_exist(file):
    if os.path.isfile(file) and os.stat(file).st_size != 0:
        return True
    else:
        return False

def polish_contig(prefix, reads, raw_contig, polished_contig, thread, preset, round, contig_name):
    # polish consensus
    polish_thread = str(min(thread, 4))
    bam = raw_contig + ".bam"
    k = 0
    while True:
        tmp_contig = raw_contig + ".tmp"
        command = "minimap2 -t " + polish_thread + " -ax " + preset + " -r2k " + raw_contig + " " + reads + " | samtools sort -@" + polish_thread + " > " + bam
        try:
            subprocess.run(command, shell = True, timeout = 300, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
        except subprocess.TimeoutExpired:
            print("fail to map reads to contig: " + contig_name)
            return
        command = "samtools view -F0x900 " + bam + " | wtpoa-cns -t " + polish_thread + " -d " + raw_contig + " -i - -fo " + tmp_contig
        try:
            subprocess.run(command, shell = True, timeout = 300, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
        except subprocess.TimeoutExpired:
            print("fail to polish contig: " + contig_name)
            return
        raw_contig = tmp_contig
        k = k + 1
        if k >= round:
            break
    if os.path.isfile(tmp_contig) and os.stat(tmp_contig).st_size != 0:
        os.rename(tmp_contig, polished_contig)
    else:
        print("polishing failed for "+contig_name+"\n")
    return
    
def annotate_contig(asm_dir, TE_library, vcf, out, sample_name, thread, presets):
    print("Generating contig annotation...")
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
            if ".cns.fa" in file and os.stat(asm_dir+"/"+file).st_size > 0:
                contig_name = file.replace('.cns.fa', '')
                with open(asm_dir+"/"+file, "r") as handle:
                    records = SeqIO.parse(handle, "fasta")
                    for record in records:
                        if record.id == 'ctg1':
                            record.id = contig_name
                            record.description = "len="+str(len(record.seq))
                            SeqIO.write(record, output_contigs, 'fasta')
                            output_list.write(contig_name+"\n")

    # map sequence to contigs
    # print ("Map VCF sequences to contigs...")
    seq2contig_out=out+"/"+"seq2contig.paf"
    if os.path.isfile(seq2contig_out):
        os.remove(seq2contig_out)

    # constrct fasta from parsed vcf file
    ins_seqs=out+"/"+sample_name+".ins.filter.fasta"
    # print ("Generating filtered VCF sequences...")
    vcf2fasta(vcf, ins_seqs)
    # print ("Done\n")

    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            vcf_seq=entry[7]
            sv_len=entry[3]
            query=out+"/"+contig_name+".seq.fa"
            create_fa(contig_name, vcf_seq, query)
            subject = out+"/"+contig_name+".fa"
            with open(subject, "w") as output:
                    subprocess.call(["samtools", "faidx", merge_contigs, contig_name], stdout=output)
            # subject=asm_dir+"/"+contig_name+".cns.fa"
            if os.path.isfile(subject):
                with open(seq2contig_out, "a") as output:
                    subprocess.call(["minimap2", "-cx", presets_minimap2, "--secondary=no", "-v", "0", subject, query], stdout=output)
            os.remove(query)
            os.remove(subject)
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
                    subprocess.call(["minimap2", "-cx", presets_minimap2, contig, TE_library, "-v", "0", "-t", str(thread)], stdout=output)
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
    # print ("Filter TE-contig alignment...")
    with open(te2contig_filter_raw, "w") as output:
        subprocess.call(command, shell = True, stdout=output)
    # print ("Done")

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
    print ("Repeatmask TE sequences...")
    with open(te_rm_out, "r") as input, open(te2contig_rm, "w") as output:
        for line in input:
            if "##" not in line:
                entry = line.replace('\n', '').split("\t")
                contig_name = entry[0].rsplit(':', 1)[0]
                start = entry[0].rsplit(':', 1)[1].split("-")[0]
                end = entry[0].rsplit(':', 1)[1].split("-")[1]
                # contigs = entry[0].replace(':', '-').split("-")
                family = re.sub('Target \"Motif:|\".*', '', entry[8])
                strand = entry[6]
                score = entry[5]
                out_line="\t".join([contig_name, start, end, family, score, strand])
                output.write(out_line+"\n")
    print("Done\n")

    te2contig_rm_merge=out+"/"+sample_name+".te2contig_rm.merge.bed"
    command="bedtools merge -c 4,6 -o distinct -delim \"|\" -i "+te2contig_rm
    with open(te2contig_rm_merge, "w") as output:
        subprocess.call(command, shell = True, stdout=output)
    
    
    return te2contig_filter_bed, te2contig_rm_merge, te_fa, merge_contigs
    
def find_te(contigs, ref, te_contigs_annotation, family_annotation, te_fa, out, sample_name, gap, overlap, presets):
    if presets == "ont":
        presets_minimap2 = "map-ont"
    else:
        presets_minimap2 = "map-pb"

    # lift over
    te_report, te_report_meta = liftover.lift_annotation(fasta1 = contigs, fasta2 = ref, bed = te_contigs_annotation, sample_name = sample_name, out_dir = out, preset = presets_minimap2, overlap = overlap, gap = gap, flank_len = 500, family_rm = family_annotation)

    # convert meta to dict
    ins_dict = dict()
    with open(te_report_meta, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            ins_dict[entry[5]] = [entry[0], entry[1], entry[2], entry[3], entry[4]]

    # generate TE sequence fasta
    final_te_seqs = out+"/"+sample_name+".final.fa"
    if os.path.isfile(final_te_seqs):
        os.remove(final_te_seqs)

    with open(te_fa, "r") as input, open(final_te_seqs, "a") as output:
        for record in SeqIO.parse(input, "fasta"):
            ins_name = record.id
            if ins_name in ins_dict:
                chr = ins_dict[ins_name][0]
                start = ins_dict[ins_name][1]
                end = ins_dict[ins_name][2]
                family = ins_dict[ins_name][3]
                te_strand = ins_dict[ins_name][4]
                record.id = chr+"_"+str(start)+"_"+str(end)+"#"+family
                if te_strand == "+" or te_strand == ".":
                    output.write(">"+record.id+"\n"+str(record.seq)+"\n")
                else:
                    output.write(">"+record.id+"\n"+str(record.seq.reverse_complement())+"\n")

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

def remove_vcf_redundancy(vcf_in, vcf_out):
    header = ["chr", "start", "end", "length", "coverage", "AF", "ID", "seq", "ids", "filter"]
    df = pd.read_csv(vcf_in, delimiter="\t", names=header)
    df2 = df.groupby(['chr', 'start', 'end']).agg({'length': 'first', 'coverage': 'sum', 'AF': af_sum, 'ID': 'first', 'seq': 'first', 'ids': id_merge, 'filter': 'first'}).reset_index()
    df2.to_csv(vcf_out, sep='\t', header=False, index=False)

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

main()
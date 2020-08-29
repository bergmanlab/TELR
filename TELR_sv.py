import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
from TELR_utility import mkdir


def detect_sv(bam, reference, out, sample_name, thread, svim=False):
    print("Generating SV output...")
    if svim:
        try:
            subprocess.call(["svim", "alignment", "--insertion_sequences", "--read_names",
                             "--sample", sample_name, "--duplications_as_insertions", out, bam, reference])
        except Exception as e:
            print(e)
            print("svim run failed, please check input bam file, exiting....")
            sys.exit(1)
        vcf = out + "/" + "variants.vcf"
    else:
        vcf = out + "/" + sample_name + ".vcf"
        # command="sniffles -n -1 --genotype --report_seq --report_BND --threads " + str(thread) + " -m " + bam + " -v " + vcf
        command = "sniffles -n -1 --report_BND --threads " + \
            str(thread) + " -m " + bam + " -v " + vcf
        try:
            subprocess.call(command, shell=True)
        except Exception as e:
            print(e)
            print("sniffles run failed, please check input bam file, exiting....")
            sys.exit(1)
    print("Done\n")
    if os.path.isfile(vcf) == False:
        sys.stderr.write("VCF file from SV detection not found, exiting....\n")
        sys.exit(1)
    else:
        return(vcf)


def sniffle_parse(vcf, out, sample_name, raw_reads, TE_library, thread):
    contig_reads_dir = out+"/"+"contig_reads"
    mkdir(contig_reads_dir)

    vcf_parsed = out+"/"+sample_name+".vcf.parsed"
    print("Generating parsed sniffle VCF...")
    query_str = "\"%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\n\""
    command = "bcftools query -i \'SVTYPE=\"INS\"\' -f "+query_str+" "+vcf
    print(command)
    with open(vcf_parsed, "w") as output:
        subprocess.call(command, stdout=output, shell=True)
    # TODO check whether vcf file contains insertions, quit if 0
    # remove redundancy in parsed vcf
    vcf_unique = out+"/"+sample_name+".vcf.reduced"
    rm_vcf_redundancy(vcf_parsed, vcf_unique)
    # os.remove(vcf_parsed)
    vcf_parsed = vcf_unique
    print("Done\n")

    # constrct fasta from parsed vcf file
    ins_seqs = out+"/"+sample_name+".ins.fasta"
    vcf2fasta(vcf_parsed, ins_seqs)

    # run RM on the inserted seqeunce
    print("Repeatmask VCF sequences...")
    subprocess.call(["RepeatMasker", "-dir", out, "-gff", "-s", "-nolow", "-no_is",
                     "-xsmall", "-e", "ncbi", "-lib", TE_library, "-pa", str(thread), ins_seqs])
    print("Done\n")
    ins_rm_out = ins_seqs+".out.gff"

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
    vcf_parsed_filtered = out+"/"+sample_name+".vcf.parsed.filtered"
    filter_vcf(vcf_parsed, te_seqs_dict, vcf_parsed_filtered)
    # TODO check how many TEs in filtered set, quit if 0
    # os.remove(vcf_parsed)

    # constrct fasta from parsed filtered vcf file
    ins_te_seqs = out+"/"+sample_name+".ins.te.fasta"
    vcf2fasta(vcf_parsed_filtered, ins_te_seqs)

    # extract read IDs
    # TODO: convert to set?
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
    sub_fa = out+"/"+sample_name+".subreads.fa"
    command = "seqtk subseq "+raw_reads+" "+ids_unique+" | seqtk seq -a"
    with open(sub_fa, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # reorder reads
    sub_fa_reorder = out+"/"+sample_name+".subreads.reorder.fa"
    extract_reads(sub_fa, ids, sub_fa_reorder)

    # separate reads into multiple files, using csplit
    work_dir = os.getcwd()
    # print ("current working dir:"+work_dir)
    os.chdir(contig_reads_dir)
    # print ("New working dir:"+os.getcwd())
    m = []
    k = 1
    with open(vcf_parsed_filtered, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            read_list = entry[8].split(",")
            k = k + 2*(len(read_list))
            m.append(k)
    if len(m) == 1:
        subprocess.call(["cp", sub_fa_reorder, "contig0"])
    elif len(m) == 0:
        print("No insertion detected, exiting...")
    else:
        m = m[:-1]
        index = " ".join(str(i) for i in m)
        command = "csplit -s -f contig -n 1 "+sub_fa_reorder+" "+index
        subprocess.call(command, shell=True)
        # print ("Done\n")

    os.remove(ids)
    os.remove(ids_unique)
    os.remove(sub_fa)
    os.remove(sub_fa_reorder)
    os.chdir(work_dir)
    # print ("New working dir:"+os.getcwd())

    return vcf_parsed_filtered, contig_reads_dir


def rm_vcf_redundancy(vcf_in, vcf_out):
    header = ["chr", "start", "end", "length",
              "coverage", "AF", "ID", "seq", "ids", "filter"]
    df = pd.read_csv(vcf_in, delimiter="\t", names=header)
    df2 = df.groupby(['chr', 'start', 'end']).agg(
        {'length': 'first', 'coverage': 'sum', 'AF': af_sum, 'ID': 'first', 'seq': 'first', 'ids': id_merge, 'filter': 'first'}).reset_index()
    df2.to_csv(vcf_out, sep='\t', header=False, index=False)


def filter_vcf(vcf, dict, out, cov=0):
    with open(vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            length = entry[3]
            if length == "999999999":
                length = 999
            else:
                length = int(length)
            if contig_name in dict:
                # no VCF sequence TE coverage filtering for now
                if dict[contig_name] >= cov * length:
                    output.write(line)


def extract_reads(reads, list, out):
    record_dict = SeqIO.index(reads, "fasta")
    with open(out, "wb") as output_handle, open(list, "r") as ID:
        for entry in ID:
            entry = entry.replace('\n', '')
            output_handle.write(record_dict.get_raw(entry))


def vcf2fasta(vcf, out):
    with open(vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            coord = '_'.join([entry[0], entry[1], entry[2]])
            output.write(">"+coord+"\n")
            output.write(entry[7]+"\n")


def id_merge(strings):
    string_merged = ','.join(strings)
    ids = string_merged.split(",")
    id_string = ','.join(get_unique_list(ids))
    return(id_string)


def get_unique_list(list1):
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))
    return(unique_list)


def af_sum(nums):
    sum_nums = sum(nums)
    if sum_nums > 1:
        sum_nums = int(sum_nums)
    return(sum_nums)

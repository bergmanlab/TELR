import sys
import os
import subprocess
import pandas as pd
import logging
import time
from Bio import SeqIO
from STELR_utility import mkdir, format_time, create_loci_set, get_contig_name


def detect_sv(
    bam,
    reference,
    vcf_out,
    sample_name,
    thread,
    sv_detector = "Sniffles"
):
    """
    Detect structural variants from BAM file using Sniffles or SVIM
    """
    logging.info("Detecting SVs from BAM file...")
    start_time = time.perf_counter()
    command = {
        #SVIM: Not implemented
        "SVIM":["svim","alignment","--insertion_sequences","--read_names","--sample",sample_name,"--interspersed_duplications_as_insertions",'.',bam,reference],
        #Sniffles: Version 1.0.12 | current 2.0.7
        "Sniffles":["sniffles", "-n", "-1", "--threads", thread, "-m", bam, "-v", vcf_out]
    }[sv_detector]
    try:
        subprocess.call(command)
    except Exception as e:
        print(e)
        logging.exception(f"Detecting SVs using {sv_detector} failed, exiting...")
        sys.exit(1)
    if(sv_detector == "SVIM"):
        vcf_tmp = os.path.join("variants.vcf")
        os.rename(vcf_tmp, vcf)
    proc_time = time.perf_counter() - start_time
    if os.path.isfile(vcf_out) is False:
        sys.stderr.write("SV detection output not found, exiting...\n")
        sys.exit(1)
    else:
        logging.info("SV detection finished in " + format_time(proc_time))

def merge_vcf(merge_in, merge_processed):
    with open(merge_in, "r") as input, open(merge_processed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if ";" in entry[3]:
                chr = entry[0]
                start = average(entry[3])
                end = average(entry[4])
                len_list = entry[5].split(";")
                idx = len_list.index(max(len_list))
                length = len_list[idx]
                coverage = sum(string2int(entry[6].split(";"), integer=False))
                af = af_sum(string2int(entry[7].split(";"), integer=False))
                sv_id = entry[8].split(";")[idx]
                ins_seq = entry[9].split(";")[idx]
                ids = entry[10].replace(";", ",").split(",")
                reads = ",".join(get_unique_list(ids))
                sv_filter = entry[11].split(";")[idx]
                gt = entry[12].split(";")[idx]
                ref_count = entry[13].split(";")[idx]
                alt_count = len(reads.split(","))
                ins_te_prop = entry[15].split(";")[idx]
                out_line = "\t".join(
                    [
                        chr,
                        str(start),
                        str(end),
                        str(length),
                        str(coverage),
                        str(af),
                        sv_id,
                        ins_seq,
                        reads,
                        sv_filter,
                        gt,
                        str(ref_count),
                        str(alt_count),
                        str(ins_te_prop),
                    ]
                )
            else:
                entry = [entry[0]] + entry[3:]
                out_line = "\t".join(entry)
            output.write(out_line + "\n")

def string2int(lst, integer=True):
    if integer:
        for i in range(0, len(lst)):
            lst[i] = int(lst[i])
    else:
        for i in range(0, len(lst)):
            lst[i] = float(lst[i])
    return lst

def average(lst):
    num_list = lst.split(";")
    num_list = string2int(num_list)
    return round(sum(num_list) / len(num_list))

def swap_coordinate(vcf_in, vcf_out, sv_detector="Sniffles1"):
    with open(vcf_in, "r") as input, open(vcf_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if int(entry[2]) < int(entry[1]):
                entry[1], entry[2] = entry[2], entry[1]
            if sv_detector == "Sniffles2":
                entry.insert(4,str(len(entry[7].split(","))))
            out_line = "\t".join(entry)
            output.write(f"{out_line}\n")

#TODO reimplement sniffles1 version

def rm_vcf_redundancy(vcf_in, vcf_out):
    header = [
        "chr",
        "start",
        "end",
        "length",
        "coverage",
        "AF",
        "ID",
        "seq",
        "reads",
        "filter",
        "genotype",
        "ref_count",
        "alt_count",
    ]
    df = pd.read_csv(vcf_in, delimiter="\t", names=header)
    df2 = (
        df.groupby(["chr", "start", "end"])
        .agg(
            {
                "length": "first",
                "coverage": "sum",
                "AF": af_sum,
                "ID": "first",
                "seq": "first",
                "reads": id_merge,
                "filter": "first",
                "genotype": "first",
                "ref_count": "sum",
                "alt_count": "sum",
            }
        )
        .reset_index()
    )
    df2.to_csv(vcf_out, sep="\t", header=False, index=False)

def repeatmask(repeatmasker_dir, ins_seqs, te_library, thread):
    mkdir(repeatmasker_dir)
    try:
        subprocess.call(
            [
                "RepeatMasker",
                "-dir",
                repeatmasker_dir,
                "-gff",
                "-s",
                "-nolow",
                "-no_is",
                "-xsmall",
                "-e",
                "ncbi",
                "-lib",
                te_library,
                "-pa",
                str(thread),
                ins_seqs,
            ]
        )
        ins_repeatmasked = f"{repeatmasker_dir}/{os.path.basename(ins_seqs)}.out.gff"
        open(ins_repeatmasked, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking VCF insertion sequences failed, exiting...")
        sys.exit(1)

def te_extract(parsed_vcf, ins_seqs, ins_rm_merge, ins_filtered, loci_eval):
    # get the length of the insertion sequence TODO: this can be generalized
    contig_len = dict()
    if os.path.isfile(ins_seqs):
        with open(ins_seqs, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                contig_len[record.id] = len(record.seq)
    
    # extract VCF sequences that contain TEs
    ins_te_loci = dict()
    with open(ins_rm_merge, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            length = int(entry[2]) - int(entry[1])
            ins_te_prop = round(length / contig_len[contig_name], 2)
            if contig_name in ins_te_loci:
                ins_te_loci[contig_name] = ins_te_loci[contig_name] + ins_te_prop
            else:
                ins_te_loci[contig_name] = ins_te_prop

    with open(parsed_vcf, "r") as input, open(ins_filtered, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = get_contig_name(entry)
            # TODO: maybe add filter for insertion sequences covered by TE?
            if contig_name in ins_te_loci:
                out_line = line.replace('\n', '') + f"\t{ins_te_loci[contig_name]}"
                output.write(f"{out_line}\n")
                
    with open(loci_eval, "a") as output:
        for locus in create_loci_set(parsed_vcf):
            if locus not in ins_te_loci:
                output.write("\t".join([locus, "VCF sequence not repeatmasked"]) + "\n")

def write_ins_seqs(parsed_vcf, out):
    with open(parsed_vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            output.write(f">{get_contig_name(entry)}\n")
            output.write(f"{entry[7]}\n")

def id_merge(strings):
    string_merged = ",".join(strings)
    ids = string_merged.split(",")
    id_string = ",".join(get_unique_list(ids))
    return id_string


def get_unique_list(list1):
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = list(list_set)
    return unique_list


def af_sum(nums):
    sum_nums = sum(nums)
    if sum_nums > 1:
        sum_nums = 1
    return sum_nums

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
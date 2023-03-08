import sys
import os
import subprocess
import pandas as pd
import logging
import time
from Bio import SeqIO
from telr.TELR_utility import mkdir, format_time, create_loci_set


def detect_sv(
    vcf,
    bam,
    reference,
    out,
    sample_name,
    thread,
    sv_detector = "Sniffles", #used to be svim = False, need to make sure this doesn't break anything
):
    """
    Detect structural variants from BAM file using Sniffles or SVIM

    Notes: this function used to have SVIM=False as an input parameter instead of sv_detector = "Sniffles":
    - Sniffles version in telr.yml is 1.0.12; current version is 2.0.7
    - As far as I can tell, SVIM is not implemented at all; it is not present in the conda environment; it is not possible to call it using telr commands.

    """
    logging.info("Detecting SVs from BAM file...")
    start_time = time.perf_counter()
    process_args = {
        #SVIM: Not implemented
        "SVIM":["svim","alignment","--insertion_sequences","--read_names","--sample",sample_name,"--interspersed_duplications_as_insertions",out,bam,reference],
        #Sniffles: Version 1.0.12 | current 2.0.7
        "Sniffles":["sniffles", "-n", "-1", "--threads", str(thread), "-m", bam, "-v", vcf]
    }[sv_detector]
    try:
        subprocess.call(process_args)
    except Exception as e:
        print(e)
        logging.exception(f"Detecting SVs using {sv_detector} failed, exiting...")
        sys.exit(1)
    if(sv_detector == "SVIM"):
        vcf_tmp = os.path.join(out, "variants.vcf")
        os.rename(vcf_tmp, vcf)
    proc_time = time.perf_counter() - start_time
    if os.path.isfile(vcf) is False:
        sys.stderr.write("SV detection output not found, exiting...\n")
        sys.exit(1)
    else:
        logging.info("SV detection finished in " + format_time(proc_time))


def vcf_parse_filter(
    vcf_in, vcf_out, bam, te_library, out, sample_name, thread, loci_eval
):
    """Parse and filter for insertions from VCF file"""
    logging.info("Parse structural variant VCF...")

    vcf_parsed = f"{vcf_in}.parsed.tmp.tsv"
    parse_vcf(vcf_in, vcf_parsed, bam)

    vcf_filtered = f"{vcf_in}.filtered.tmp.tsv"
    filter_vcf(
        vcf_parsed, vcf_filtered, te_library, out, sample_name, thread, loci_eval
    )

    # merge entries
    vcf_merged = vcf_in + ".merged.tmp.tsv"
    merge_vcf(vcf_filtered, vcf_merged)
    os.rename(vcf_merged, vcf_out)


def merge_vcf(vcf_in, vcf_out, window=20):
    vcf_tmp = vcf_in + ".tmp"
    with open(vcf_tmp, "w") as output:
        command = (
            'bedtools merge -o collapse -c 2,3,4,5,6,7,8,9,10,11,12,13,14 -delim ";"'
            + " -d "
            + str(window)
            + " -i "
            + vcf_in
        )
        subprocess.call(command, shell=True, stdout=output)

    with open(vcf_tmp, "r") as input, open(vcf_out, "w") as output:
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
                reads = ",".join(set(ids)) #changed from get_unique_list to set()
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

    os.remove(vcf_tmp)


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


def parse_vcf(vcf_in, vcf_out, bam):
    '''Parse vcf file using bcftools'''
    query_str = '"%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n"'
    command = (f'bcftools query -i \'SVTYPE="INS" & ALT!="<INS>"\' -f {query_str} {vcf_in}')
    #bcftools version 1.9 | current 1.16
    with open(f"{vcf_in}.tmp", "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # check start and end, swap if needed
    swap_coordinate(f"{vcf_in}.tmp", f"{vcf_in}.swap")

    # sort bed file

    # TODO check whether vcf file contains insertions, quit if 0
    rm_vcf_redundancy(f"{vcf_in}.swap", vcf_out)  # remove redundancy in parsed vcf
    os.remove(f"{vcf_in}.swap")
    os.remove(f"{vcf_in}.tmp")


def swap_coordinate(vcf_in, vcf_out):
    '''checks if column 2 (0 index) of the parsed vcf file > column 1; if so, swaps columns 1 and 2 (start and end pos).'''
    with open(vcf_in, "r") as input, open(vcf_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if int(entry[2]) < int(entry[1]):
                entry[1], entry[2] = entry[2], entry[1]
            out_line = "\t".join(entry)
            output.write(out_line + "\n")


def rm_vcf_redundancy(vcf_in, vcf_out):
    '''Takes the output of parse_vcf, uses pandas dataframe to intelligently 
    remove redundancy without removing additional information (eg counts), and 
    returns a non-redundant vcf file containing the relevant information.'''
    header = [
        "chr",
        "start",
        "end",
        "length",
        "coverage", #from sniffles, RE = read support
        "AF", #allele frequency
        "ID",
        "seq",
        "reads", #ids of all mapped reads
        "filter",
        "genotype",
        "ref_count", # of high quality reference reads
        "alt_count", # of high quality variant reads
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


def filter_vcf(ins, ins_filtered, te_library, out, sample_name, thread, loci_eval):
    """
    Filter insertion sequences from Sniffles VCF by repeatmasking with TE consensus
    """
    # construct fasta from parsed vcf file
    ins_seqs = os.path.join(out, sample_name.replace("+","plus") + ".vcf_ins.fasta")
    write_ins_seqs(ins, ins_seqs)

    # get the length of the insertion sequence TODO: this can be generalized
    contig_len = dict()
    if os.path.isfile(ins_seqs):
        with open(ins_seqs, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                contig_len[record.id] = len(record.seq)

    # run RM on the inserted seqeunce
    repeatmasker_dir = os.path.join(out, "vcf_ins_repeatmask")
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
        ins_repeatmasked = os.path.join(
            repeatmasker_dir, os.path.basename(ins_seqs) + ".out.gff"
        )
        open(ins_repeatmasked, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking VCF insertion sequences failed, exiting...")
        sys.exit(1)

    # sort RM gff
    ins_rm_sort = os.path.join(
        repeatmasker_dir, os.path.basename(ins_seqs) + ".out.sort.gff"
    )
    with open(ins_rm_sort, "w") as output:
        subprocess.call(["bedtools", "sort", "-i", ins_repeatmasked], stdout=output)

    # merge RM gff
    ins_rm_merge = os.path.join(
        repeatmasker_dir, os.path.basename(ins_seqs) + ".out.merge.bed"
    )
    with open(ins_rm_merge, "w") as output:
        subprocess.call(["bedtools", "merge", "-i", ins_rm_sort], stdout=output)

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

    with open(ins, "r") as input, open(ins_filtered, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            # TODO: maybe add filter for insertion sequences covered by TE?
            if contig_name in ins_te_loci:
                out_line = line.replace("\n", "") + "\t" + str(ins_te_loci[contig_name])
                output.write(out_line + "\n")
    # os.remove(ins_seqs)

    # report removed loci
    with open(loci_eval, "a") as output:
        for locus in create_loci_set(ins):
            if locus not in ins_te_loci:
                output.write("\t".join([locus, "VCF sequence not repeatmasked"]) + "\n")


def write_ins_seqs(vcf, out):
    '''Takes file in intermediate '''
    with open(vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            coord = "_".join([entry[0], entry[1], entry[2]])
            output.write(f">{coord}\n{entry[7]}\n")


def id_merge(strings):
    '''returns comma separated list of unique read IDs.'''
    string_merged = ",".join(strings)
    ids = set(string_merged.split(","))
    id_string = ",".join(ids)
    return id_string


def af_sum(nums):
    '''Returns sum of allele frequencies (max 1)'''
    sum_nums = sum(nums)
    if sum_nums > 1:
        sum_nums = 1
    return sum_nums

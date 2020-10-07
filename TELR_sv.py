import sys
import os
import subprocess
import pandas as pd
import logging
import time
from TELR_utility import mkdir, format_time, create_loci_set


def detect_sv(
    vcf,
    bam,
    reference,
    te_library,
    out,
    sample_name,
    thread,
    svim=False,
):
    """
    Detect structural variants from BAM file using Sniffles or SVIM
    """
    logging.info("Detecting SVs from BAM file...")
    start_time = time.time()
    if svim:
        try:
            subprocess.call(
                [
                    "svim",
                    "alignment",
                    "--insertion_sequences",
                    "--read_names",
                    "--sample",
                    sample_name,
                    "--interspersed_duplications_as_insertions",
                    out,
                    bam,
                    reference,
                ]
            )
        except Exception as e:
            print(e)
            logging.exception("Detecting SVs using SVIM failed, exiting...")
            sys.exit(1)
        vcf_tmp = os.path.join(out, "variants.vcf")
        os.rename(vcf_tmp, vcf)
    else:
        try:
            subprocess.call(
                ["sniffles", "-n", "-1", "--threads", str(thread), "-m", bam, "-v", vcf]
            )
        except Exception as e:
            print(e)
            logging.exception("Detecting SVs using Sniffles failed, exiting...")
            sys.exit(1)
    proc_time = time.time() - start_time
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

    vcf_parsed = vcf_in + ".parsed.tmp.tsv"
    parse_vcf(vcf_in, vcf_parsed, bam)

    vcf_filtered = vcf_in + ".filtered.tmp.tsv"
    filter_vcf(
        vcf_parsed, vcf_filtered, te_library, out, sample_name, thread, loci_eval
    )

    # merge entries
    vcf_merged = vcf_in + ".merged.tmp.tsv"
    merge_vcf(vcf_filtered, vcf_merged)

    os.rename(vcf_merged, vcf_out)
    # os.remove(vcf_parsed_tmp)


def merge_vcf(vcf_in, vcf_out, window=20):
    vcf_tmp = vcf_in + ".tmp"
    with open(vcf_tmp, "w") as output:
        command = (
            'bedtools merge -o collapse -c 2,3,4,5,6,7,8,9,10,11,12,13 -delim ";"'
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
                reads = ",".join(get_unique_list(ids))
                sv_filter = entry[11].split(";")[idx]
                gt = entry[12].split(";")[idx]
                ref_count = entry[13].split(";")[idx]
                alt_count = len(reads.split(","))
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
    vcf_tmp = vcf_in + ".tmp"
    query_str = '"%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n"'
    command = (
        'bcftools query -i \'SVTYPE="INS" & ALT!="<INS>"\' -f '
        + query_str
        + " "
        + vcf_in
    )
    with open(vcf_tmp, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # check start and end, swap if needed
    vcf_swap = vcf_in + ".swap"
    swap_coordinate(vcf_tmp, vcf_swap)

    # sort bed file

    # TODO check whether vcf file contains insertions, quit if 0
    rm_vcf_redundancy(vcf_swap, vcf_out)  # remove redundancy in parsed vcf
    os.remove(vcf_swap)
    os.remove(vcf_tmp)


def swap_coordinate(vcf_in, vcf_out):
    with open(vcf_in, "r") as input, open(vcf_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if entry[2] < entry[1]:
                entry[1], entry[2] = entry[2], entry[1]
            out_line = "\t".join(entry)
            output.write(out_line + "\n")


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


def filter_vcf(ins, ins_filtered, te_library, out, sample_name, thread, loci_eval):
    """
    Filter insertion sequences from Sniffles VCF by repeatmasking with TE concensus
    """
    # constrct fasta from parsed vcf file
    ins_seqs = os.path.join(out, sample_name + ".vcf_ins.fasta")
    write_ins_seqs(ins, ins_seqs)

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

    # extract VCF sequences that contain TEs
    with open(ins_repeatmasked, "r") as input:
        ins_te_loci = {
            line.replace("\n", "").split("\t")[0]
            for line in input
            if "RepeatMasker" in line
        }

    with open(ins, "r") as input, open(ins_filtered, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            # TODO: maybe add filter for insertion sequences covered by TE?
            if contig_name in ins_te_loci:
                output.write(line)
    os.remove(ins_seqs)

    # report removed loci
    with open(loci_eval, "a") as output:
        for locus in create_loci_set(ins):
            if locus not in ins_te_loci:
                output.write("\t".join([locus, "VCF sequence not repeatmasked"]) + "\n")


def write_ins_seqs(vcf, out):
    with open(vcf, "r") as input, open(out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            coord = "_".join([entry[0], entry[1], entry[2]])
            output.write(">" + coord + "\n")
            output.write(entry[7] + "\n")


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

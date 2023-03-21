import sys
import os
import subprocess
import pandas as pd
import logging
import time
from Bio import SeqIO
from telr.TELR_utility import mkdir, format_time, create_loci_set

class file_container:
    def __init__(self, sample_name, directory_dict):
        self.sample_name = sample_name
        self.directories = directory_dict

    def add(self, key, directory, extension, **kwargs):
        file_name = self.sample_name
        if "new_dir" in kwargs:
            self.directories.update(kwargs.pop("new_dir"))
        if("file_name" in kwargs): file_name = kwargs["file_name"]
        self.__dict__[key] = file(
            self, directory,
            name=f"{file_name}{extension}",
            **kwargs
        )
    
    def add_dir(self, directory_dict):
        self.directories.update(directory_dict)
    
    def extend(self, file, key, extension, **kwargs):
        file_name = file.name
        if "file_name" in kwargs: file_name = kwargs["file_name"]
        directory = file.directory
        if "new_dir" in kwargs:
            directory = kwargs.pop("new_dir")
            if type(directory) is dict:
                self.directories.update(directory)
                directory = directory[next(iter(directory))]
        self.add(key, directory, extension, file_name, **kwargs)
        return self.__dict__[key]

    

class file:
    def __init__(self, container, directory, name, **kwargs):
        self.container = container
        self.directory = directory
        self.name = name
        self.path = os.path.join(self.container.directories[directory], self.name)
        self.file_format = self.name[self.name.rindex(".")+1:]
        self.__dict__.update(kwargs)
    
    def add(self, **kwargs):
        self.__dict__.update(kwargs)
    
    def extend(self, key, extension, **kwargs):
        return self.container.extend(self, key, extension, **kwargs)
    
    def exists(self):
        return os.path.isfile(self.path)
    
    def open(self, options="r"):
        return open(self.path, options)

def detect_sv(
    tmp_dir,
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
    sv_files = file_container(sample_name,{"tmp":tmp_dir})
    sv_files.add(key="sv_raw", directory="tmp", extension=".vcf", file_format="vcf", note=f"raw output from {sv_detector}")
    start_time = time.perf_counter()
    process_args = {
        #SVIM: Not implemented
        "SVIM":["svim","alignment","--insertion_sequences","--read_names","--sample",sample_name,"--interspersed_duplications_as_insertions",out,bam,reference],
        #Sniffles: Version 1.0.12 | current 2.0.7
        "Sniffles":["sniffles", "-n", "-1", "--threads", str(thread), "-m", bam, "-v", sv_files.sv_raw]
    }[sv_detector]
    try:
        subprocess.call(process_args)
    except Exception as e:
        print(e)
        logging.exception(f"Detecting SVs using {sv_detector} failed, exiting...")
        sys.exit(1)
    if(sv_detector == "SVIM"):
        vcf_tmp = os.path.join(out, "variants.vcf")
        os.rename(vcf_tmp, sv_files.sv_raw)
    proc_time = time.perf_counter() - start_time
    if os.path.isfile(sv_files.sv_raw) is False:
        sys.stderr.write("SV detection output not found, exiting...\n")
        sys.exit(1)
    else:
        logging.info("SV detection finished in " + format_time(proc_time))
        return sv_files


def vcf_parse_filter(
    sv_files, bam, te_library, out, sample_name, thread, loci_eval
):
    """Parse and filter for insertions from VCF file"""
    logging.info("Parse structural variant VCF...")

    sv_files.add(key="parsed",directory="tmp",extension=".parsed.tmp.tsv",file_format="tsv")
    parse_vcf(sv_files.sv_raw, sv_files.parsed)
    
    sv_files.add(key="parsed",directory="tmp",extension=".filtered.tmp.tsv",file_format="tsv")
    sv_files.add_dir()
    filter_vcf(
        sv_files, te_library, sample_name, thread, loci_eval
    )

    # merge entries
    vcf_merged = f"{vcf_in}.merged.tmp.tsv"
    merge_vcf(vcf_filtered, vcf_merged)
    os.rename(vcf_merged, vcf_out)


def merge_vcf(vcf_in, vcf_out, window=20):
    '''merges insertions within 20 bp of each other using bedtools + post-processing in TELR'''
    #as far as okg can tell, this merge loses quite a bit of information about all but the longest read being merged.
    vcf_tmp = f"{vcf_in}.tmp"
    with open(vcf_tmp, "w") as output:
        command = [
            "bedtools", "merge",
            "-o", "collapse",
            "-c", "2,3,4,5,6,7,8,9,10,11,12,13,14",
            "-d", str(window),
            "-delim",'";"',
            "-i", vcf_in
        ]
        subprocess.call(command, stdout=output)

    with open(vcf_tmp, "r") as input, open(vcf_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if ";" in entry[3]:
                chr = entry[0]
                start = average(entry[3])
                end = average(entry[4])
                len_list = entry[5].split(";")
                longest_insertion = len_list.index(max(len_list))
                length = len_list[longest_insertion]
                coverage = sum(string2int(entry[6], integer=False))
                af = af_sum(string2int(entry[7], integer=False))
                sv_id = entry[8].split(";")[longest_insertion]
                ins_seq = entry[9].split(";")[longest_insertion]
                ids = entry[10].replace(";", ",").split(",")
                reads = ",".join(set(ids)) #changed from get_unique_list to set()
                sv_filter = entry[11].split(";")[longest_insertion]
                gt = entry[12].split(";")[longest_insertion]
                ref_count = entry[13].split(";")[longest_insertion]
                alt_count = len(reads.split(","))
                ins_te_prop = entry[15].split(";")[longest_insertion]
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


def string2int(string, integer=True, delimiter=";"):
    '''takes a ; (by default) delimited string and returns as a list of either ints or floats'''
    lst = string.split(delimiter)
    if integer:
        for i in range(len(lst)):
            lst[i] = int(lst[i])
    else:
        for i in range(len(lst)):
            lst[i] = float(lst[i])
    return lst


def average(lst):
    '''returns the average of a ; delimited string with a list of ints'''
    num_list = string2int(lst)
    return round(sum(num_list) / len(num_list))


def parse_vcf(vcf_in, vcf_out):
    '''Parse vcf file using bcftools'''
    query_str = '"%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n"'
    command = (f'bcftools query -i \'SVTYPE="INS" & ALT!="<INS>"\' -f {query_str} {vcf_in.path}')
    bcftools_raw_output = vcf_in.extend("bcftools_raw", ".tmp", file_format = "nonstandard")
    #bcftools version 1.9 | current 1.16
    with open(bcftools_raw_output.path, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # check start and end, swap if needed
    bcftools_swapped_output = vcf_in.extend("bcftools_swap", ".swap", file_format = "nonstandard")
    swap_coordinate(bcftools_raw_output.path, bcftools_swapped_output.path)

    # sort bed file

    # TODO check whether vcf file contains insertions, quit if 0
    rm_vcf_redundancy(bcftools_swapped_output.path, vcf_out.path)  # remove redundancy in parsed vcf
    os.remove(bcftools_swapped_output.path)
    os.remove(bcftools_raw_output.path)


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


def filter_vcf(sv_files, te_library, sample_name, thread, loci_eval):
    """
    Filter insertion sequences from Sniffles VCF by repeatmasking with TE consensus
    """
    # construct fasta from parsed vcf file
    sv_files.add("ins_seqs",directory="tmp",extension=".vcf_ins.fasta", file_name = sample_name.replace("+","plus"), file_format = "fasta")
    write_ins_seqs(sv_files.parsed.path, sv_files.ins_seqs.path)

    # get the length of the insertion sequence TODO: this can be generalized
    contig_len = dict()
    if sv_files.ins_seqs.exists():
        with sv_files.ins_seqs.open() as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                contig_len[record.id] = len(record.seq)

    # run RM on the inserted seqeunce
    sv_files.add_dir({"rpmask":os.path.join(sv_files.directories["tmp"],"vcf_ins_repeatmask")})
    mkdir(sv_files.directories["rpmask"])
    try:
        subprocess.call(
            [
                "RepeatMasker",
                "-dir",
                sv_files.directories["rpmask"],
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
                sv_files.ins_seqs.path
            ]
        )
        sv_files.ins_seqs.extend("ins_repeatmasked", ".out.gff", new_dir = "rpmask")
        sv_files.ins_repeatmasked.open()
    except Exception as e:
        print(e)
        print("Repeatmasking VCF insertion sequences failed, exiting...")
        sys.exit(1)

    # sort RM gff
    sv_files.ins_seqs.extend("ins_rm_sort", ".out.sort.gff", new_dir = "rpmask")
    with sv_files.ins_rm_sort.open("w") as output:
        subprocess.call(["bedtools", "sort", "-i", sv_files.ins_repeatmasked.path], stdout=output)

    # merge RM gff
    sv_files.ins_seqs.extend("ins_rm_merge", ".out.merge.bed", new_dir = "rpmask")
    with sv_files.ins_rm_merge.path.open("w") as output:
        subprocess.call(["bedtools", "merge", "-i", sv_files.ins_rm_sort.path], stdout=output)

    # extract VCF sequences that contain TEs
    ins_te_loci = dict()
    with sv_files.ins_rm_merge.open() as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            length = int(entry[2]) - int(entry[1])
            ins_te_prop = round(length / contig_len[contig_name], 2)
            if contig_name in ins_te_loci:
                ins_te_loci[contig_name] = ins_te_loci[contig_name] + ins_te_prop
            else:
                ins_te_loci[contig_name] = ins_te_prop

    with sv_files.parsed.open() as input, sv_files.ins_filtered.open("w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            # TODO: maybe add filter for insertion sequences covered by TE?
            if contig_name in ins_te_loci:
                output.write(f"{line.strip()}\t{ins_te_loci[contig_name]}\n")
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

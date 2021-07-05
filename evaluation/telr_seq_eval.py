#!/usr/bin/env python3
import sys
import argparse
import os
import json
import glob
import subprocess
import pandas as pd
import time
from multiprocessing import Pool
from eval_utility import filter_annotation, check_exist, format_time, create_soft_link


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to evaluate TELR contig and TE sequence quality"
    )

    ## required ##
    parser.add_argument(
        "-i",
        "--telr_out_dir",
        type=str,
        help="directory that include telr runs",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        help="sample genome assembly in fasta format",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--annotation",
        type=str,
        help="TE annotation for sample genome assembly",
        required=True,
    )
    parser.add_argument(
        "--region1",
        type=str,
        help="Regions in reference genome assembly to include in the analysis in BED format",
        required=False,
    )
    parser.add_argument(
        "--region2",
        type=str,
        help="Regions in sample genome assembly to include in the analysis in BED format",
        required=False,
    )
    parser.add_argument(
        "--include_families",
        type=str,
        help="TE families to use in the analysis (separated by comma)",
        required=False,
    )
    parser.add_argument(
        "--exclude_families",
        type=str,
        help="TE families to exclude in the analysis (separated by comma)",
        required=False,
    )
    parser.add_argument(
        "--flank_len",
        type=int,
        help="Size of flanking sequences to include when aligning contig to sample genome assembly (default = '500bp')",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        type=str,
        help="output directory (default = 'current working directory')",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="number of cores (default = '1')",
        required=False,
    )

    args = parser.parse_args()

    # setup default value for optional arguments
    if args.out_dir is None:
        args.out_dir = "."
    args.out_dir = os.path.abspath(args.out_dir)
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    if args.flank_len is None:
        args.flank_len = 500

    if args.threads is None:
        args.threads = 1

    return args


def extract_contig_seqs(
    contigs, contig_name, contig_len, start, end, out_dir, flank_len
):
    subject = os.path.join(out_dir, contig_name + ".contig.fa")
    start_new = start - flank_len
    if start_new < 0:
        start_new = 0
    end_new = end + flank_len
    if end_new > contig_len:
        end_new = contig_len

    coord = contig_name + ":" + str(start_new + 1) + "-" + str(end_new)
    with open(subject, "w") as output:
        subprocess.call(
            ["samtools", "faidx", contigs, coord],
            stdout=output,
        )
    return subject, start_new, end_new, end_new - start_new


def get_ref_fa(ref, chromosome, start, end, prefix, out_dir):
    subject = os.path.join(out_dir, prefix + ".fa")
    coord = chromosome + ":" + str(start + 1) + "-" + str(end)
    with open(subject, "w") as output:
        subprocess.call(
            ["samtools", "faidx", ref, coord],
            stdout=output,
        )
    return subject


def make_fa(contig_name, seq, out_dir):
    contig_te_seq_fa = out_dir + "/" + contig_name + ".te.seq.fa"
    with open(contig_te_seq_fa, "w") as output:
        output.write(">" + contig_name + "\n")
        output.write(seq)
    return contig_te_seq_fa


def dict_merge(dict1, dict2):
    merged_dict = {**dict1, **dict2}
    return merged_dict


def get_paf_info(paf):
    """
    Extract information from minimap2 PAF file, only output info from longest alignment block
    """
    query_length = None
    num_query_subject_hits = None
    subject_chrom = None
    subject_start = None
    subject_end = None
    query_mapp_qual = None
    query_max_base_mapped_prop = None
    num_residue_matches = None
    alignment_block_length = None
    blast_identity = None

    if check_exist(paf):
        with open(paf, "r") as input:
            k = 0  # number of alignments
            max_alignment_block_length = 0  # proportion of querys that can be aligned
            for line in input:
                entry = line.replace("\n", "").split("\t")
                query_length = int(entry[1])
                alignment_block_length = int(entry[10])
                # choose the longest alignment
                if max_alignment_block_length < alignment_block_length:
                    max_alignment_block_length = alignment_block_length
                    num_residue_matches = int(entry[9])
                    query_max_base_mapped_prop = float(
                        num_residue_matches / query_length
                    )
                    query_mapp_qual = int(entry[11])
                    alignment_block_length = int(entry[10])
                    blast_identity = float(num_residue_matches / alignment_block_length)
                    subject_chrom = entry[5]
                    subject_start = int(entry[7])
                    subject_end = int(entry[8])
                k = k + 1
            num_query_subject_hits = k
    return (
        query_length,
        num_query_subject_hits,
        subject_chrom,
        subject_start,
        subject_end,
        query_mapp_qual,
        query_max_base_mapped_prop,
        num_residue_matches,
        alignment_block_length,
        blast_identity,
    )


def get_paftools_var_summary(paf, contig_name, out_dir):
    """
    Use Paftools.js to call variants, output variant summary
    """
    num_covered = None
    num_snv = None
    num_1bp_ins = None
    num_1bp_del = None
    num_2bp_ins = None
    num_2bp_del = None
    num_50bp_ins = None
    num_50bp_del = None
    num_1kb_ins = None
    num_1kb_del = None
    num_ins = None
    num_del = None
    num_indel = None

    if check_exist(paf):
        variant = out_dir + "/" + contig_name + ".var.txt"
        summary = out_dir + "/" + contig_name + ".var.summary.txt"
        command = "paftools.js call -l 100 -L 100 " + paf
        with open(variant, "w") as VAR, open(summary, "w") as SUM:
            subprocess.call(command, stdout=VAR, shell=True, stderr=SUM)

        with open(summary, "r") as input:
            for line in input:
                num = int(line.split(" ")[0])
                if "reference bases covered" in line:
                    num_covered = num
                if "substitutions" in line:
                    num_snv = num
                if "1bp deletions" in line:
                    num_1bp_del = num
                if "1bp insertions" in line:
                    num_1bp_ins = num
                if "2bp deletions" in line:
                    num_2bp_del = num
                if "2bp insertions" in line:
                    num_2bp_ins = num
                if "[3,50) deletions" in line:
                    num_50bp_del = num
                if "[3,50) insertions" in line:
                    num_50bp_ins = num
                if "[50,1000) deletions" in line:
                    num_1kb_del = num
                if "[50,1000) insertions" in line:
                    num_1kb_ins = num

        num_ins = num_1bp_ins + num_2bp_ins + num_50bp_ins + num_1kb_ins
        num_del = num_1bp_del + num_2bp_del + num_50bp_del + num_1kb_del
        num_indel = num_ins + num_del

    return (
        num_covered,
        num_snv,
        num_1bp_ins,
        num_1bp_del,
        num_2bp_ins,
        num_2bp_del,
        num_50bp_ins,
        num_50bp_del,
        num_1kb_ins,
        num_1kb_del,
        num_ins,
        num_del,
        num_indel,
    )


def compare_contig_ref_te_seqs(fa1, fa2, contig_name, out_dir):
    """
    Compare TE sequences from contigs and sample genome assembly
    """
    paf = out_dir + "/" + contig_name + ".diff.paf"
    with open(paf, "w") as output:
        subprocess.call(
            ["minimap2", "-cx", "asm5", "--cs", fa1, fa2], stdout=output
        )  # TODO: which preset I should use
    sort_paf = out_dir + "/" + contig_name + ".sort.paf"
    with open(sort_paf, "w") as output:
        command = "sort -k6,6 -k8,8n " + paf
        subprocess.call(command, shell=True, stdout=output)

    # get alignment stats from PAF file
    (
        contig_te_length,
        num_contig_ref_te_hits,
        ref_te_aligned_chrom,
        ref_te_aligned_start,
        ref_te_aligned_end,
        contig_te_mapp_qual,
        contig_te_max_base_mapped_prop,
        contig_te_num_residue_matches,
        contig_te_alignment_block_length,
        contig_te_blast_identity,
    ) = get_paf_info(sort_paf)

    # check if file exists
    (
        ref_te_num_bases_covered,
        ref_te_num_snvs,
        ref_te_num_1bp_ins,
        ref_te_num_1bp_del,
        ref_te_num_2bp_ins,
        ref_te_num_2bp_del,
        ref_te_num_50bp_ins,
        ref_te_num_50bp_del,
        ref_te_num_1kb_ins,
        ref_te_num_1kb_del,
        ref_te_num_ins,
        ref_te_num_del,
        ref_te_num_indel,
    ) = get_paftools_var_summary(paf, contig_name, out_dir)

    # seqs_similarity = dict_merge(paf_info, paftools_var_summary)

    seqs_similarity = {
        "contig_te_length": contig_te_length,
        "num_contig_ref_te_hits": num_contig_ref_te_hits,
        "ref_te_aligned_chrom": ref_te_aligned_chrom,
        "ref_te_aligned_start": ref_te_aligned_start,
        "ref_te_aligned_end": ref_te_aligned_end,
        "contig_te_mapp_qual": contig_te_mapp_qual,
        "contig_te_max_base_mapped_prop": contig_te_max_base_mapped_prop,
        "contig_te_num_residue_matches": contig_te_num_residue_matches,
        "contig_te_alignment_block_length": contig_te_alignment_block_length,
        "contig_te_blast_identity": contig_te_blast_identity,
        "ref_te_num_bases_covered": ref_te_num_bases_covered,
        "ref_te_num_snvs": ref_te_num_snvs,
        "ref_te_num_bases_covered": ref_te_num_bases_covered,
        "ref_te_num_1bp_ins": ref_te_num_1bp_ins,
        "ref_te_num_1bp_del": ref_te_num_1bp_del,
        "ref_te_num_2bp_ins": ref_te_num_2bp_ins,
        "ref_te_num_2bp_del": ref_te_num_2bp_del,
        "ref_te_num_50bp_ins": ref_te_num_50bp_ins,
        "ref_te_num_50bp_del": ref_te_num_50bp_del,
        "ref_te_num_1kb_ins": ref_te_num_1kb_ins,
        "ref_te_num_1kb_del": ref_te_num_1kb_del,
        "ref_te_num_ins": ref_te_num_ins,
        "ref_te_num_del": ref_te_num_del,
        "ref_te_num_indel": ref_te_num_indel,
    }

    return seqs_similarity


def check_contig_overlap_annotation(
    contig_name,
    annotation,
    ref_aligned_chrom,
    ref_aligned_start,
    ref_aligned_end,
    contig_te_family,
):
    """Check if TELR contigs cover TE annotations with same family in the sample genome assembly"""
    ref_te_chrom = None
    ref_te_start = None
    ref_te_end = None
    ref_te_family = None
    pass_filter = False
    if (
        ref_aligned_chrom is not None
        and ref_aligned_start is not None
        and ref_aligned_end is not None
    ):
        df = pd.read_csv(annotation, sep="\t", comment="t", header=None)
        header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
        df.columns = header[: len(df.columns)]

        # also check same families?
        df_filtered = df[
            (
                (
                    (df["chrom"] == ref_aligned_chrom)
                    & (df["chromStart"] >= int(ref_aligned_start))
                    & (df["chromEnd"] <= int(ref_aligned_end))
                    & (df["name"] == contig_te_family)
                )
                | (
                    (df["chrom"] == ref_aligned_chrom)
                    & (df["chromStart"] <= int(ref_aligned_end))
                    & (df["chromEnd"] >= int(ref_aligned_end))
                    & (df["name"] == contig_te_family)
                )
                | (
                    (df["chrom"] == ref_aligned_chrom)
                    & (df["chromStart"] <= int(ref_aligned_start))
                    & (df["chromEnd"] >= int(ref_aligned_start))
                    & (df["name"] == contig_te_family)
                )
            )
        ]
        # num_contig_ref_tes_covered = int(df_filtered.shape[0])
        if df_filtered.shape[0] > 1:
            print(contig_name + ":more than one TEs in the aligned region")
        elif df_filtered.shape[0] == 0:
            print(contig_name + ":no ref TEs in the aligned region")
        else:
            ref_te_chrom = df_filtered["chrom"].iloc[0]
            ref_te_start = int(df_filtered["chromStart"].iloc[0])
            ref_te_end = int(df_filtered["chromEnd"].iloc[0])
            ref_te_family = df_filtered["name"].iloc[0]
            pass_filter = True

    return (
        pass_filter,
        ref_te_chrom,
        ref_te_start,
        ref_te_end,
        ref_te_family,
        # num_contig_ref_tes_covered,
    )


def run_eval(args):
    contigs_fa = args[0]
    contig_name = args[1]
    contig_len = args[2]
    contig_te_start = args[3]
    contig_te_end = args[4]
    contig_te_family = args[5]
    seq = args[6]
    reference = args[7]
    annotation = args[8]
    out_dir = args[9]
    flank_len = args[10]

    # initiate final report
    eval_report = {
        "contig_name": contig_name,
        "contig_length": contig_len,
        "align_contig_start": None,
        "align_contig_end": None,
        "align_contig_size": None,
        "num_contig_ref_hits": None,
        "ref_aligned_chrom": None,
        "ref_aligned_start": None,
        "ref_aligned_end": None,
        "contig_max_base_mapped_prop": None,
        "contig_mapp_qual": None,
        "contig_num_residue_matches": None,
        "contig_alignment_block_length": None,
        "contig_blast_identity": None,
        "contig_te_family": contig_te_family,
        "contig_te_start": contig_te_start,
        "contig_te_end": contig_te_end,
        "contig_te_length": None,
        "ref_te_family": None,
        "ref_te_length": None,
        "num_contig_ref_te_hits": None,
        "ref_te_aligned_chrom": None,
        "ref_te_aligned_start": None,
        "ref_te_aligned_end": None,
        "contig_te_mapp_qual": None,
        "contig_te_max_base_mapped_prop": None,
        "contig_te_num_residue_matches": None,
        "contig_te_alignment_block_length": None,
        "contig_te_blast_identity": None,
        "ref_te_num_bases_covered": None,
        "ref_te_num_snvs": None,
        "ref_te_num_1bp_del": None,
        "ref_te_num_1bp_ins": None,
        "ref_te_num_2bp_del": None,
        "ref_te_num_2bp_ins": None,
        "ref_te_num_50bp_del": None,
        "ref_te_num_50bp_ins": None,
        "ref_te_num_1kb_del": None,
        "ref_te_num_1kb_ins": None,
        "ref_te_num_ins": None,
        "ref_te_num_del": None,
        "ref_te_num_indel": None,
    }

    # extract subsequences from contigs including flanks
    (
        contig_fa,
        eval_report["align_contig_start"],
        eval_report["align_contig_end"],
        eval_report["align_contig_size"],
    ) = extract_contig_seqs(
        contigs_fa,
        contig_name,
        contig_len,
        contig_te_start,
        contig_te_end,
        out_dir,
        flank_len,
    )

    # align TELR contig to sample genome assembly
    contig2asm_out = out_dir + "/" + contig_name + ".paf"
    preset = "asm10"  # TODO: think about the best one
    with open(contig2asm_out, "w") as output:
        subprocess.call(
            [
                "minimap2",
                "-cx",
                preset,
                "-v",
                "0",
                "--secondary=no",
                reference,
                contig_fa,
            ],
            stdout=output,
        )
    # get info from contig to ref alignment
    (
        contig_length,
        num_contig_ref_hits,
        ref_aligned_chrom,
        ref_aligned_start,
        ref_aligned_end,
        contig_mapp_qual,
        contig_max_base_mapped_prop,
        contig_num_residue_matches,
        contig_alignment_block_length,
        contig_blast_identity,
    ) = get_paf_info(contig2asm_out)

    contig_ref_align_info = {
        "num_contig_ref_hits": num_contig_ref_hits,
        "ref_aligned_chrom": ref_aligned_chrom,
        "ref_aligned_start": ref_aligned_start,
        "ref_aligned_end": ref_aligned_end,
        "contig_mapp_qual": contig_mapp_qual,
        "contig_max_base_mapped_prop": contig_max_base_mapped_prop,
        "contig_num_residue_matches": contig_num_residue_matches,
        "contig_alignment_block_length": contig_alignment_block_length,
        "contig_blast_identity": contig_blast_identity,
    }

    eval_report.update(contig_ref_align_info)

    # check if contig can be uniquely aligned to reference
    # if contig can be aligned to reference, do something
    (
        pass_filter,
        ref_te_chrom,
        ref_te_start,
        ref_te_end,
        ref_te_family,
        # num_contig_ref_tes_covered,
    ) = check_contig_overlap_annotation(
        contig_name,
        annotation,
        ref_aligned_chrom,
        ref_aligned_start,
        ref_aligned_end,
        contig_te_family,
    )

    # if contig is uniquely aligned to reference, compare TE sequeces betweeen assembly and reference
    if pass_filter:
        prefix = contig_name + "_refTE"
        ref_te_fa = get_ref_fa(
            reference,
            ref_te_chrom,
            ref_te_start,
            ref_te_end,
            prefix,
            out_dir,
        )
        len_ref_te = int(ref_te_end - ref_te_start)

        contig_te_fa = make_fa(contig_name, seq, out_dir)

        # Get alignment info from PAF file and VAR summary from PAFtools call
        seqs_similarity = compare_contig_ref_te_seqs(
            ref_te_fa, contig_te_fa, contig_name, out_dir
        )
        os.remove(ref_te_fa)
        os.remove(contig_te_fa)

        eval_report.update(seqs_similarity)
        eval_report["ref_te_length"] = len_ref_te
        eval_report["ref_te_family"] = ref_te_family

    contig_stats_json = out_dir + "/" + contig_name + ".eval.json"
    with open(contig_stats_json, "w") as output:
        json.dump(eval_report, output, indent=4, sort_keys=False)


def main():
    args = get_args()

    # get prefix
    prefix = "telr_seq_eval"

    # Get contig TE annotation files
    pattern = "/**/*te2contig_rm.merge.bed"
    contig_te_annotation_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    contig_te_annotation = contig_te_annotation_list[0]

    # Get contig fasta file
    pattern = "/**/*contigs.fa"
    contigs_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    contigs_file = contigs_list[0]
    # Create soft link for contig fasta, build index
    contigs = create_soft_link(contigs_file, args.out_dir)
    subprocess.call(["samtools", "faidx", contigs])

    # Get TELR final output
    pattern = "/**/*telr.json"
    json_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    json_file = json_list[0]

    # Get TELR final bed
    pattern = "/**/*telr.bed"
    bed_list = glob.glob(args.telr_out_dir + pattern, recursive=True)
    pred_file = bed_list[0]

    # filter TELR predictions by region and family
    pred_filtered = args.out_dir + "/" + prefix + ".parse.filter.bed"
    filter_annotation(
        bed_in=pred_file,
        bed_out=pred_filtered,
        filter_region=args.region1,
        include_families=args.include_families,
        exclude_families=args.exclude_families,
    )
    # get the coordinate ID set
    coord_ids = set()
    with open(pred_filtered, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[0]
            start = entry[1]
            end = entry[2]
            family = entry[3]
            coord_id = "_".join([chrom, str(start), str(end), family])
            coord_ids.add(coord_id)

    # get contig TE annotation info
    contig_info = dict()
    with open(contig_te_annotation, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            contig_te_start = int(entry[1])
            contig_te_end = int(entry[2])
            contig_info[contig_name] = dict()
            contig_info[contig_name]["start"] = contig_te_start
            contig_info[contig_name]["end"] = contig_te_end

    # get contig length
    with open(contigs, "r") as input:
        for line in input:
            if ">" in line:
                entry = line.replace("\n", "").split(" ")
                contig_name = entry[0].replace(">", "")
                contig_len = int(entry[1].replace("len=", ""))
                if contig_name in contig_info:
                    contig_info[contig_name]["length"] = contig_len

    # load TELR json output
    with open(json_file) as f:
        telr_out = json.load(f)
    # look through each contig in TELR output and prepare for parallel run
    contig_ids = set()
    eval_pa_list = []
    for item in telr_out:
        contig_name = item["ID"]
        family = item["family"]
        chrom = item["chr"]
        start = item["start"]
        end = item["end"]
        coord_id = "_".join([chrom, str(start), str(end), family])
        if coord_id in coord_ids:
            contig_ids.add(contig_name)
            contig_len = contig_info[contig_name]["length"]
            contig_te_start = contig_info[contig_name]["start"]
            contig_te_end = contig_info[contig_name]["end"]
            seq = item["sequence"]
            eval_pa = [
                contigs,
                contig_name,
                contig_len,
                contig_te_start,
                contig_te_end,
                family,
                seq,
                args.ref,
                args.annotation,
                args.out_dir,
                args.flank_len,
            ]
            eval_pa_list.append(eval_pa)
    # run contig evaluation in parallel
    print("Perform TELR sequence evaluation...")
    start_time = time.time()
    try:
        pool = Pool(processes=args.threads)
        pool.map(run_eval, eval_pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("TELR sequence evaluation failed, exiting...")
        sys.exit(1)

    proc_time = time.time() - start_time
    print("TELR sequence evaluation finished in " + format_time(proc_time))

    # report stats
    contig_stats_all = []
    for item in telr_out:
        contig_name = item["ID"]
        if contig_name in contig_ids:
            contig_stats_json = args.out_dir + "/" + contig_name + ".eval.json"
            with open(contig_stats_json) as f:
                contig_stats = json.load(f)
                contig_stats["telr_family"] = item["family"]
                contig_stats["telr_support_type"] = item["support_type"]
                contig_stats["telr_frequency"] = item["frequency"]
                contig_stats["telr_alt_count"] = item["alt_count"]
                contig_stats["telr_ref_count"] = item["gt"].split(":")[1]
                contig_stats_all.append(contig_stats)
            os.remove(contig_stats_json)

    contig_stats_all_json = args.out_dir + "/" + "eval_stat.json"
    with open(contig_stats_all_json, "w") as output:
        json.dump(contig_stats_all, output, indent=4, sort_keys=False)


if __name__ == "__main__":
    main()

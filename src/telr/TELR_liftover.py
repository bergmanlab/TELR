#!/usr/bin/env python3

__version__ = "1.1"
__author__ = "Shunhua Han"

import sys
import argparse
import os
import json
import subprocess
import shutil
from multiprocessing import Pool


def get_args(program_version, arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description="Program for lifting TE annotations from one contig/assembly to another, this program requires SAMtools, bedtools and minimap2."
    )

    ## required ##
    parser.add_argument("--fasta1", type=str, help="genome 1", required=True)
    parser.add_argument("--fasta2", type=str, help="genome 2", required=True)
    parser.add_argument(
        "-1",
        "--bed1",
        type=str,
        help="annotation file for genome 1 in bed format",
        required=True,
    )
    parser.add_argument(
        "-2",
        "--bed2",
        type=str,
        help="annotation file for genome 2 in bed format",
        required=False,
    )

    ## optional ##
    parser.add_argument(
        "-x",
        "--preset",
        type=str,
        help="minimap2 preset used in the flank to reference alignment step (default = 'asm10')",
        required=False,
    )
    parser.add_argument(
        "-l",
        "--flank_len",
        type=int,
        help="flanking sequence length (default = '500')",
        required=False,
    )
    parser.add_argument(
        "-g",
        "--flank_gap_max",
        type=int,
        help="maximum gap size between 5p and 3p flanking sequence alignments (default = '50')",
        required=False,
    )
    parser.add_argument(
        "-p",
        "--flank_overlap_max",
        type=int,
        help="maximum overlap size between 5p and 3p flanking sequence alignments (default = '50')",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        help="output directory (default = 'Current working directory')",
        required=False,
    )
    parser.add_argument(
        "-k",
        "--keep_files",
        action="store_true",
        help="If provided then all intermediate files will be kept (default: remove intermediate files)",
        required=False,
    )
    # parser.add_argument(
    #     "--single_flank",
    #     action="store_true",
    #     help="If provided then allow single flank support (default: requires both flanks to be aligned to reference genome)",
    #     required=False,
    # )
    parser.add_argument(
        "--different_contig_name",
        action="store_true",
        help="If provided then TELR does not require the contig name to match before and after annotation liftover (default: require contig name to be the same before and after liftover)",
        required=False,
    )
    parser.add_argument(
        "--telr_mode",
        action="store_true",
        help="If provided then run as a module in TELR (default: run in standalone mode)",
        required=False,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="number of cores (default = '1')",
        required=False,
    )

    ## info ##
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )

    args = parser.parse_args(arguments)

    # checks if required files exist
    try:
        test = open(args.fasta1, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.fasta2, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        test = open(args.bed1, "r")
    except Exception as e:
        print(e)
        sys.exit(1)

    # set default values for optional arguments
    if args.flank_gap_max is None:
        args.flank_gap_max = 50

    if args.flank_overlap_max is None:
        args.flank_overlap_max = 50

    if args.flank_len is None:
        args.flank_len = 500

    if args.preset is None:
        args.preset = "asm10"

    if args.threads is None:
        args.threads = 1

    return args


def get_ref_seq(ref, chrom, start, end, out_dir):
    """
    Extract subsequence from fasta based on coordinates
    """
    prefix = "_".join([chrom, str(start), str(end), "subseq"])
    bed = os.path.join(out_dir, prefix + ".bed")
    subject = os.path.join(out_dir, prefix + ".fa")
    with open(bed, "w") as output:
        out_line = "\t".join([chrom, str(start), str(end)])
        output.write(out_line + "\n")
    with open(subject, "w") as output:
        subprocess.call(
            ["bedtools", "getfasta", "-fi", ref, "-bed", bed],
            stdout=output,
        )
    subseq = ""
    with open(subject, "r") as input:
        for line in input:
            if ">" not in line:
                entry = line.replace("\n", "")
                subseq = subseq + entry
    # os.remove(bed)
    # os.remove(subject)
    return subseq


def check_exist(file):
    if file:
        if os.path.isfile(file) and os.stat(file).st_size != 0:
            return True
        else:
            return False
    else:
        return False


def extract_genome_seqs(ref, ref_size_dict, chrom, start, end, prefix, out_dir):
    """
    Extract subsequence from fasta based on coordinates
    """
    ref_size = ref_size_dict[chrom]

    if end > ref_size or start < 0:
        return None
    else:
        prefix_bed = "_".join([chrom, str(start), str(end)])
        bed = os.path.join(out_dir, prefix_bed + ".fa")
        subject = os.path.join(out_dir, prefix + ".fa")
        with open(bed, "w") as output:
            out_line = "\t".join([chrom, str(start), str(end)])
            output.write(out_line + "\n")
        with open(subject, "w") as output:
            subprocess.call(
                ["bedtools", "getfasta", "-fi", ref, "-bed", bed],
                stdout=output,
            )
        os.remove(bed)
        return subject


def paf_to_bed(paf, bed, filter=None):
    """
    convert PAF to sorted BED format
    """
    bed_tmp = bed + ".tmp"
    with open(paf, "r") as input, open(bed_tmp, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[5]
            if filter is not None:
                if chrom == filter:
                    start = entry[7]
                    end = entry[8]
                    name = entry[0]
                    score = entry[11]
                    strand = entry[4]
                    out_line = "\t".join([chrom, start, end, name, score, strand])
                    output.write(out_line + "\n")
            else:
                start = entry[7]
                end = entry[8]
                name = entry[0]
                score = entry[11]
                strand = entry[4]
                out_line = "\t".join([chrom, start, end, name, score, strand])
                output.write(out_line + "\n")

    if check_file(bed_tmp):
        with open(bed, "w") as output:
            subprocess.call(["bedtools", "sort", "-i", bed_tmp], stdout=output)
            os.remove(bed_tmp)


def align_flank(flank_fa, ref_fa, paf_out, preset, num_secondary):
    """
    Align flanking sequences to genome 2
    """
    with open(paf_out, "w") as output:
        subprocess.call(
            [
                "minimap2",
                "-cx",
                preset,
                "-v",
                "0",
                "-N",
                str(num_secondary),
                ref_fa,
                flank_fa,
            ],
            stdout=output,
        )


def get_coord(start_3p, end_3p, start_5p, end_5p, strand):
    if strand == "+":
        start = end_3p
        end = start_5p
    else:
        start = start_3p
        end = end_5p
    gap = end - start
    if start > end:
        start, end = end, start
    return start, end, gap


def create_bed(chrom, start, end, family, strand, filename):
    with open(filename, "w") as output:
        out_line = "\t".join([chrom, str(start), str(end), family, ".", strand])
        output.write(out_line)


def check_nearby_ref(
    chrom, start_query, end_query, family, strand, ref_bed, out_dir, threshold=5000
):
    """
    Check if flanking seqeunce alignments have nearby TE annotations in genome 2
    """
    distance = None

    if check_exist(ref_bed):
        prefix_query = "_".join([chrom, str(start_query), str(end_query), "query"])
        query_bed = os.path.join(out_dir, prefix_query + ".bed")
        create_bed(chrom, start_query, end_query, family, strand, query_bed)

        prefix_overlap = "_".join([chrom, str(start_query), str(end_query), "nearby"])
        overlap = os.path.join(out_dir, prefix_overlap + ".bed")
        with open(overlap, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "closest",
                    "-a",
                    query_bed,
                    "-b",
                    ref_bed,
                    "-d",
                    "-D",
                    "ref",
                    "-k",
                    "5"
                    # "-s",
                ],
                stdout=output,
            )

        with open(overlap, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                chrom1 = entry[0]
                chrom2 = entry[6]
                family1 = entry[3]
                family2 = entry[9]
                strand1 = entry[5]
                strand2 = entry[11]
                if chrom1 == chrom2 and family1 == family2 and strand1 == strand2:
                    distance_new = int(entry[12])
                    if distance is None:
                        distance = distance_new
                    else:
                        distance = absmin(distance, distance_new)
    if distance is not None:
        if abs(distance) > threshold:
            distance = None
    return distance


def absmin(num1, num2):
    """
    Return the value with min aboslute
    """
    num1_abs = abs(num1)
    num2_abs = abs(num2)
    num_abs_min = min(num1_abs, num2_abs)
    if num_abs_min == num1_abs:
        return num1
    else:
        return num2


def get_paf_info(paf):
    """
    Extract information from minimap2 PAF file
    """
    paf_info = dict()
    if check_exist(paf):
        with open(paf, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                query_length = int(entry[1])
                num_residue_matches = int(entry[9])
                alignment_block_length = int(entry[10])
                query_mapp_qual = int(entry[11])
                alignment_block_length = int(entry[10])
                sequence_identity = float(num_residue_matches / alignment_block_length)
                paf_entry_id = "_".join([entry[0], entry[5], entry[7], entry[8]])
                paf_info[paf_entry_id] = {
                    "query_length": query_length,
                    "query_mapp_qual": query_mapp_qual,
                    "num_residue_matches": num_residue_matches,
                    "alignment_block_length": alignment_block_length,
                    "sequence_identity": sequence_identity,
                }

    return paf_info


def get_genome_size(fasta):
    genome_index = fasta + ".fai"
    genome_size_dict = dict()
    with open(genome_index, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            genome_size_dict[entry[0]] = int(entry[1])
    return genome_size_dict


def run_liftover_single_annotation(input_json):
    """
    Main liftover function for a single annotation
    """
    # initiate liftover report dictionary
    lift_entries = dict()

    # load input data from json
    with open(input_json) as f:
        input_values = json.load(f)
    chrom = input_values["chrom"]
    start = input_values["start"]
    end = input_values["end"]
    family = input_values["family"]
    strand = input_values["strand"]
    fasta1 = input_values["fasta1"]
    fasta2 = input_values["fasta2"]
    out_dir = input_values["out_dir"]
    flank_len = input_values["flank_len"]
    flank_gap_max = input_values["flank_gap_max"]
    flank_overlap_max = input_values["flank_overlap_max"]
    bed2 = input_values["bed2"]  # could be None
    preset = input_values["preset"]
    # single_flank = input_values["single_flank"]
    different_contig_name = input_values["different_contig_name"]
    telr_mode = input_values["telr_mode"]

    prefix = "_".join([chrom, str(start), str(end), family])
    prefix = prefix.replace("|", "_")
    lift_entries["ID"] = prefix
    genome1_coord = chrom + ":" + str(start) + "-" + str(end)
    lift_entries["genome1_coord"] = genome1_coord

    te_length = int(end) - int(start)
    lift_entries["te_length"] = te_length

    # get genome 1 size
    genome1_size_dict = get_genome_size(fasta1)

    # extract 5 prime flank
    start_5p_flank = int(start) - flank_len + 1
    end_5p_flank = int(start)
    prefix_5p = prefix + "_5p"
    fa_5p = extract_genome_seqs(
        fasta1,
        genome1_size_dict,
        chrom,
        start_5p_flank,
        end_5p_flank,
        prefix_5p,
        out_dir,
    )
    # extract 3 prime flank
    start_3p_flank = int(end) + 1
    end_3p_flank = int(end) + 1 + flank_len
    prefix_3p = prefix + "_3p"
    fa_3p = extract_genome_seqs(
        fasta1,
        genome1_size_dict,
        chrom,
        start_3p_flank,
        end_3p_flank,
        prefix_3p,
        out_dir,
    )

    # intiate minimap2 parameters
    num_secondary = 10  # save 10 secondary alignments
    if not different_contig_name:
        if telr_mode:
            filter_chrom = "_".join(chrom.split("_")[:-2])
        else:
            filter_chrom = chrom
    else:
        filter_chrom = None

    # map 5 prime flanks to genome2
    if fa_5p:
        align_5p_flank_paf = out_dir + "/" + prefix_5p + ".paf"
        align_flank(fa_5p, fasta2, align_5p_flank_paf, preset, num_secondary)
        align_5p_flank_qcs = dict()
        if check_exist(align_5p_flank_paf):
            align_5p_flank_qcs = get_paf_info(align_5p_flank_paf)

        align_5p_flank_bed = out_dir + "/" + prefix_5p + ".bed"

        paf_to_bed(align_5p_flank_paf, align_5p_flank_bed, filter=filter_chrom)
        os.remove(align_5p_flank_paf)
        os.remove(fa_5p)
    else:
        align_5p_flank_bed = None

    # map 3 prime flanks to genome2
    if fa_3p:
        align_3p_flank_paf = out_dir + "/" + prefix_3p + ".paf"
        align_flank(fa_3p, fasta2, align_3p_flank_paf, preset, num_secondary)
        align_3p_flank_qcs = dict()
        if check_exist(align_3p_flank_paf):
            align_3p_flank_qcs = get_paf_info(align_3p_flank_paf)

        align_3p_flank_bed = out_dir + "/" + prefix_3p + ".bed"
        paf_to_bed(align_3p_flank_paf, align_3p_flank_bed, filter=None)
        os.remove(align_3p_flank_paf)
        os.remove(fa_3p)
    else:
        align_3p_flank_bed = None

    # find closest entries between 3 prime and 5 prime flank alignments, require two alignments to be on the same strand, also report distance
    overlap = out_dir + "/" + prefix + ".overlap"
    if check_file(align_5p_flank_bed) and check_file(align_3p_flank_bed):
        with open(overlap, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "closest",
                    "-a",
                    align_5p_flank_bed,
                    "-b",
                    align_3p_flank_bed,
                    "-s",
                    "-d",
                    "-t",
                    "all",
                ],
                stdout=output,
            )

    # identify liftover coordinates and report in json format
    lift_entries["report"] = []
    lift_start = 0
    lift_end = 0
    num_hits = 0  # number of non-reference liftover hits
    reported = False
    if check_file(overlap):
        with open(overlap, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")

                chrom_genome2_5p = entry[0]
                chrom_genome2_3p = entry[6]
                flank_strand = entry[5]
                mapp_quality_5p = int(entry[4])
                mapp_quality_3p = int(entry[10])
                if (
                    chrom_genome2_5p == chrom_genome2_3p and chrom_genome2_3p != "."
                ):  # to make sure the entry exists
                    lift_chrom = chrom_genome2_5p
                    # get flank qc
                    align_5p_flank_id = "_".join(
                        [entry[3], entry[0], entry[1], entry[2]]
                    )
                    align_5p_flank_qc = align_5p_flank_qcs[align_5p_flank_id]
                    align_3p_flank_id = "_".join(
                        [entry[9], entry[6], entry[7], entry[8]]
                    )
                    align_3p_flank_qc = align_3p_flank_qcs[align_3p_flank_id]

                    start_5p = int(entry[1])
                    end_5p = int(entry[2])
                    start_3p = int(entry[7])
                    end_3p = int(entry[8])
                    # get final liftover coordiantes
                    lift_start, lift_end, lift_gap = get_coord(
                        start_5p, end_5p, start_3p, end_3p, flank_strand
                    )
                    # figure out the strand of the lifted annotation
                    if flank_strand == strand:
                        lift_strand = "+"
                    else:
                        lift_strand = "-"
                    # report the flanking sequence alignments
                    align_5p_coord = (
                        chrom_genome2_5p + ":" + str(start_5p) + "-" + str(end_5p)
                    )
                    align_3p_coord = (
                        chrom_genome2_3p + ":" + str(start_3p) + "-" + str(end_3p)
                    )
                    lift_entry = {
                        "type": None,
                        "family": family,
                        "chrom": lift_chrom,
                        "start": int(lift_start),
                        "end": int(lift_end),
                        "strand": lift_strand,
                        "gap": lift_gap,
                        "TSD_length": None,
                        "TSD_sequence": None,
                        "5p_flank_align_coord": align_5p_coord,
                        "5p_flank_mapping_quality": mapp_quality_5p,
                        "5p_flank_num_residue_matches": align_5p_flank_qc[
                            "num_residue_matches"
                        ],
                        "5p_flank_alignment_block_length": align_5p_flank_qc[
                            "alignment_block_length"
                        ],
                        "5p_flank_sequence_identity": align_5p_flank_qc[
                            "sequence_identity"
                        ],
                        "3p_flank_align_coord": align_3p_coord,
                        "3p_flank_mapping_quality": mapp_quality_3p,
                        "3p_flank_num_residue_matches": align_3p_flank_qc[
                            "num_residue_matches"
                        ],
                        "3p_flank_alignment_block_length": align_3p_flank_qc[
                            "alignment_block_length"
                        ],
                        "3p_flank_sequence_identity": align_3p_flank_qc[
                            "sequence_identity"
                        ],
                        "distance_5p_flank_ref_te": None,
                        "distance_3p_flank_ref_te": None,
                        "comment": None,
                    }
                    # check the distance between flank alignment and nearast TE annotation
                    distance_5p = check_nearby_ref(
                        lift_chrom,
                        start_5p,
                        end_5p,
                        family,
                        lift_strand,
                        bed2,
                        out_dir,
                    )
                    distance_3p = check_nearby_ref(
                        lift_chrom,
                        start_3p,
                        end_3p,
                        family,
                        lift_strand,
                        bed2,
                        out_dir,
                    )
                    if distance_5p is not None:
                        lift_entry["distance_5p_flank_ref_te"] = distance_5p
                    if distance_3p is not None:
                        lift_entry["distance_3p_flank_ref_te"] = distance_3p

                    # if the overlap between flank alignments is greater than 50bp, don't report
                    if lift_gap < -flank_overlap_max:
                        if not reported:
                            reported = False  # TODO: redundant
                    elif lift_gap >= -flank_overlap_max and lift_gap <= flank_gap_max:
                        # if the gap between flank alignments is smaller than 50bp
                        # report as reference if 1) there is a reference 2 TE in between two flanks or 2) the gap size and original TE size is similar or 3) the gap size is bigger than the size of the TE.
                        # report as non-reference otherwise.
                        if (
                            (
                                distance_5p is not None
                                and distance_5p >= 0
                                and distance_5p <= lift_gap
                                and distance_3p is not None
                                and distance_3p <= 0
                                and -distance_3p <= lift_gap
                            )
                            or check_nums_similar(lift_gap, te_length)
                            or (lift_gap >= te_length)
                        ):
                            lift_entry["type"] = "reference"
                            lift_entry[
                                "comment"
                            ] = "overlap/gap size between 3p and 5p flanks within threshold, include genome2 TE in between"
                        else:
                            lift_entry["type"] = "non-reference"
                            lift_entry[
                                "comment"
                            ] = "overlap/gap size between 3p and 5p flanks within threshold"
                            # get TSD length and sequence
                            if lift_gap == 0:
                                lift_entry["TSD_length"] = 0
                                lift_entry["TSD_sequence"] = None
                            if lift_gap < 0:
                                lift_entry["TSD_length"] = -lift_gap
                                lift_entry["TSD_sequence"] = get_ref_seq(
                                    fasta2, lift_chrom, lift_start, lift_end, out_dir
                                )

                            num_hits = num_hits + 1
                        lift_entries["report"].append(lift_entry)
                        reported = True
                    else:  # if the gap between flanks is greater than 50bp
                        if lift_gap > flank_gap_max and lift_gap <= 0.5 * te_length:
                            # if the gap size is smaller than half of the original TE size
                            # report as reference if 1) there is a reference 2 TE in between two flanks
                            # report as non-reference otherwise.
                            if (
                                distance_5p is not None
                                and distance_5p >= 0
                                and distance_5p <= lift_gap
                                and distance_3p is not None
                                and distance_3p <= 0
                                and -distance_3p <= lift_gap
                            ):
                                lift_entry["type"] = "reference"
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size less than half of TE annotation, include genome2 TE in between"
                            else:
                                lift_entry["type"] = "non-reference"
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size exceeds threshold but less than half of TE annotation, no genome2 TE in between"
                                num_hits = num_hits + 1
                            lift_entries["report"].append(lift_entry)
                            reported = True
                        elif lift_gap >= 0.5 * te_length and lift_gap <= 20000:
                            # if the gap size is greater than half of the original TE size and smaller than 20kb
                            lift_entry["type"] = "reference"
                            if (
                                distance_5p is not None
                                and distance_5p >= 0
                                and distance_5p <= lift_gap
                                and distance_3p is not None
                                and distance_3p <= 0
                                and -distance_3p <= lift_gap
                            ):
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size greater than half of TE annotation, include genome2 TE in between"  # TODO: check same family?
                            else:
                                lift_entry[
                                    "comment"
                                ] = "flanks gap size greater than half of TE annotation, no genome2 TE in between"
                            lift_entries["report"].append(lift_entry)
                            reported = True
                        else:
                            # if the gap size is greater than 20kb
                            if not reported:
                                reported = False  # redundant

    # more than one report per annotation, need to filter
    # first step is to choose the best reference annotation by comparing gap size with TE length
    if len(lift_entries["report"]) > 1:
        best_ref_entry = dict()
        best_nonref_entry = dict()
        for report in lift_entries["report"]:
            if report["type"] == "reference":
                if not best_ref_entry:
                    best_ref_entry = report
                else:
                    current_gap_size = best_ref_entry["gap"]
                    new_gap_size = report["gap"]
                    if choose_new_size(te_length, current_gap_size, new_gap_size):
                        best_ref_entry = report
            if report["type"] == "non-reference":
                if not best_nonref_entry:
                    best_nonref_entry = report
                else:
                    reported = False
        lift_entries["report"] = None
        if reported:
            if best_ref_entry and best_nonref_entry:
                # if both reference and non-reference liftover can be found, report the non-reference one
                # lift_entries["report"].append(best_ref_entry)
                lift_entries["report"] = best_nonref_entry
            elif best_ref_entry:
                lift_entries["report"] = best_ref_entry
            elif best_nonref_entry:
                lift_entries["report"] = best_nonref_entry
            else:
                reported = False
    elif len(lift_entries["report"]) == 1:
        lift_entries["report"] = lift_entries["report"][0]

    if not reported:
        lift_entry = {
            "type": "unlifted",
            "family": family,
            "chrom": None,
            "start": None,
            "end": None,
            "strand": None,
            "gap": None,
            "TSD_length": None,
            "TSD_sequence": None,
            "5p_flank_align_coord": None,
            "5p_flank_mapping_quality": None,
            "5p_flank_num_residue_matches": None,
            "5p_flank_alignment_block_length": None,
            "5p_flank_sequence_identity": None,
            "3p_flank_align_coord": None,
            "3p_flank_mapping_quality": None,
            "3p_flank_num_residue_matches": None,
            "3p_flank_alignment_block_length": None,
            "3p_flank_sequence_identity": None,
            "distance_5p_flank_ref_te": None,
            "distance_3p_flank_ref_te": None,
            "comment": "flank alignments not nearby each other / only one flank aligned",
        }
        # TODO: if only one flank can be lifted, check to see if there is a reference TE nearby (same family, same strand, similar size?)
        align_5p_coords = []
        align_3p_coords = []
        if check_file(align_5p_flank_bed):
            with open(align_5p_flank_bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    align_5p_coords.append(entry[0] + ":" + entry[1] + "-" + entry[2])
        if check_file(align_3p_flank_bed):
            with open(align_3p_flank_bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    align_3p_coords.append(entry[0] + ":" + entry[1] + "-" + entry[2])

        if len(align_5p_coords) == 1:
            lift_entry["5p_flank_align_coord"] = align_5p_coords[0]
        elif len(align_5p_coords) > 1:
            lift_entry["5p_flank_align_coord"] = align_5p_coords

        if len(align_3p_coords) == 1:
            lift_entry["3p_flank_align_coord"] = align_3p_coords[0]
        elif len(align_3p_coords) > 1:
            lift_entry["3p_flank_align_coord"] = align_3p_coords

        # if single flank mode is turned on, inspect single flanks and report as lifted
        # if single_flank:
        if len(align_5p_coords) == 1 and len(align_3p_coords) == 0:
            with open(align_5p_flank_bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    flank_chrom = entry[0]
                    flank_start = int(entry[1])
                    flank_end = int(entry[2])
                    flank_mapping_quality = int(entry[4])
                    flank_strand = entry[5]
                    align_5p_flank_id = "_".join(
                        [entry[3], entry[0], entry[1], entry[2]]
                    )
            if flank_strand == strand:
                lift_strand = "+"
            else:
                lift_strand = "-"
            if flank_strand == "+":
                lift_start = lift_end = flank_end
            else:
                lift_start = lift_end = flank_start

            align_5p_flank_qc = align_5p_flank_qcs[align_5p_flank_id]

            lift_entry["chrom"] = flank_chrom
            lift_entry["start"] = int(lift_start)
            lift_entry["end"] = int(lift_end)
            lift_entry["mapp_quality_5p"] = flank_mapping_quality
            lift_entry["strand"] = lift_strand
            lift_entry["5p_flank_num_residue_matches"] = align_5p_flank_qc[
                "num_residue_matches"
            ]
            lift_entry["5p_flank_alignment_block_length"] = align_5p_flank_qc[
                "alignment_block_length"
            ]
            lift_entry["5p_flank_sequence_identity"] = align_5p_flank_qc[
                "sequence_identity"
            ]

            distance_5p = check_nearby_ref(
                flank_chrom,
                flank_start,
                flank_end,
                family,
                lift_strand,
                bed2,
                out_dir,
            )
            lift_entry["distance_5p_flank_ref_te"] = distance_5p
            # if distance between flank and ref is small, report as ref, otherwise as non-ref
            if distance_5p is not None and abs(distance_5p) <= 5:
                lift_entry["type"] = "reference"
                lift_entry[
                    "comment"
                ] = "only one flank aligned, flank alignment adjacent to reference TE"
            else:
                lift_entry["type"] = "non-reference"
                lift_entry[
                    "comment"
                ] = "only one flank aligned, flank alignment not adjacent to reference TE"
                num_hits = 1

        elif len(align_5p_coords) == 0 and len(align_3p_coords) == 1:
            with open(align_3p_flank_bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    flank_chrom = entry[0]
                    flank_start = int(entry[1])
                    flank_end = int(entry[2])
                    flank_mapping_quality = int(entry[4])
                    flank_strand = entry[5]
                    align_3p_flank_id = "_".join(
                        [entry[3], entry[0], entry[1], entry[2]]
                    )
            if flank_strand == strand:
                lift_strand = "+"
            else:
                lift_strand = "-"
            if flank_strand == "+":
                lift_start = lift_end = flank_start
            else:
                lift_start = lift_end = flank_end

            align_3p_flank_qc = align_3p_flank_qcs[align_3p_flank_id]

            lift_entry["chrom"] = flank_chrom
            lift_entry["start"] = int(lift_start)
            lift_entry["end"] = int(lift_end)
            lift_entry["mapp_quality_5p"] = flank_mapping_quality
            lift_entry["strand"] = lift_strand
            lift_entry["5p_flank_num_residue_matches"] = align_3p_flank_qc[
                "num_residue_matches"
            ]
            lift_entry["5p_flank_alignment_block_length"] = align_3p_flank_qc[
                "alignment_block_length"
            ]
            lift_entry["5p_flank_sequence_identity"] = align_3p_flank_qc[
                "sequence_identity"
            ]

            distance_3p = check_nearby_ref(
                flank_chrom,
                flank_start,
                flank_end,
                family,
                lift_strand,
                bed2,
                out_dir,
            )
            lift_entry["distance_3p_flank_ref_te"] = distance_3p
            # if distance between flank and ref is small, report as ref, otherwise as non-ref
            if distance_3p is not None and abs(distance_3p) <= 5:
                lift_entry["type"] = "reference"
                lift_entry[
                    "comment"
                ] = "only one flank aligned, flank alignment adjacent to reference TE"
            else:
                lift_entry["type"] = "non-reference"
                lift_entry[
                    "comment"
                ] = "only one flank aligned, flank alignment not adjacent to reference TE"
                num_hits = 1

        lift_entries["report"] = lift_entry
    lift_entries["num_hits"] = num_hits

    # write
    liftover_report = out_dir + "/" + prefix + "_report.json"
    with open(liftover_report, "w") as output:
        json.dump(lift_entries, output, indent=4, sort_keys=False)

    return liftover_report


def choose_new_size(size_ref, size_old, size_new):
    if size_ref - size_old > size_ref - size_new:
        return True
    else:
        return False


def check_nums_similar(num1, num2):
    normalize_diff = abs(num1 - num2) / num2
    if normalize_diff <= 0.1:
        return True
    else:
        return False


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        print("Creation of the directory %s failed" % dir)


def check_file(file):
    if file:
        if os.path.isfile(file) and os.stat(file).st_size != 0:
            return True
        else:
            return False
    else:
        return False


def build_index(fa):
    subprocess.call(["samtools", "faidx", fa])


def liftover(
    fasta1,
    fasta2,
    bed1,
    bed2,
    preset,
    flank_len,
    flank_gap_max,
    flank_overlap_max,
    out,
    threads,
    keep_files,
    # single_flank,
    different_contig_name,
    telr_mode,
):
    """
    Core function to do annotation liftover from one assembly/contig to another
    """

    # generate genome1 index file if not present
    genome1_index = fasta1 + ".fai"
    if not check_file(genome1_index):
        build_index(fasta1)

    # generate genome2 index file if not present
    genome2_index = fasta2 + ".fai"
    if not check_file(genome2_index):
        build_index(fasta2)

    # loop through TE annotations, prepare data for parallel liftover jobs
    input_json_dir = os.path.join(out, "input_json")
    mkdir(input_json_dir)
    tmp_dir = os.path.join(out, "tmp")
    mkdir(tmp_dir)

    liftover_pa_list = []
    with open(bed1, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[0]
            start = int(entry[1])
            end = int(entry[2])
            family = entry[3]
            strand = entry[5]

            prefix = "_".join([chrom, str(start), str(end), family])
            liftover_entry = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "family": family,
                "strand": strand,
                "fasta1": fasta1,
                "fasta2": fasta2,
                "out_dir": tmp_dir,
                "flank_len": flank_len,
                "flank_gap_max": flank_gap_max,
                "flank_overlap_max": flank_overlap_max,
                "bed2": bed2,
                "preset": preset,
                # "single_flank": single_flank,
                "different_contig_name": different_contig_name,
                "telr_mode": telr_mode,
            }
            liftover_input = input_json_dir + "/" + prefix + "_input.json"
            with open(liftover_input, "w") as output:
                json.dump(liftover_entry, output, indent=4, sort_keys=False)
            liftover_pa_list.append(liftover_input)

    # run liftover jobs for all annotations in parallel
    print("Perform liftover...")
    try:
        pool = Pool(processes=threads)
        liftover_report_list = pool.map(
            run_liftover_single_annotation, liftover_pa_list
        )
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("liftover failed, exiting...")
        sys.exit(1)

    # final report
    print("Generate final report...")
    data = []
    for liftover_report in liftover_report_list:
        with open(liftover_report) as f:
            report = json.load(f)
            data.append(report)

    json_report_tmp = out + "/" + "liftover_report.tmp.json"
    with open(json_report_tmp, "w") as output:
        json.dump(data, output, indent=4, sort_keys=False)

    # add a step to merge close annotations if inserions are : 1) supported by only one flank 2) same family 3) same strand

    # merge entries with overlapping or identical coordinates, choose one with longest TE length
    ## step one: write json into bed format
    bed_report_tmp = out + "/" + "liftover_report.tmp.bed"
    with open(bed_report_tmp, "w") as output:
        for entry in data:
            num_hits = entry["num_hits"]
            te_id = entry["ID"]
            te_length = entry["te_length"]
            if num_hits == 1:
                report = entry["report"]
                te_type = report["type"]
                if te_type == "non-reference":
                    chrom = report["chrom"]
                    start = report["start"]
                    end = report["end"]
                    family = report["family"]
                    strand = report["strand"]
                    score = "."
                    out_line = "\t".join(
                        [
                            chrom,
                            str(start),
                            str(end),
                            family,
                            score,
                            strand,
                            str(te_length),
                            te_id,
                        ]
                    )
                    output.write(out_line + "\n")

    ## step two: sort bed file
    bed_report_sort = out + "/" + "liftover_report.sort.bed"
    command = "bedtools sort -i " + bed_report_tmp
    with open(bed_report_sort, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    ## step three: merge entries with same coordinates
    bed_report_merge = out + "/" + "liftover_report.merge.bed"
    with open(bed_report_merge, "w") as output:
        command = (
            'bedtools merge -d 0 -o collapse -c 2,3,4,5,6,7,8 -delim "," -i '
            + bed_report_sort
        )
        subprocess.call(command, shell=True, stdout=output)

    ## step four: find overlapped entries and pick ones to filter out
    remove_ids = set()
    with open(bed_report_merge, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[0]
            if "," in entry[3]:
                len_list = entry[8].split(",")
                keep_idx = len_list.index(max(len_list))
                ids = entry[9].split(",")
                final_id = ids[keep_idx]
                for index in ids:
                    if index != final_id:
                        remove_ids.add(index)

    ## step five: write final report
    data_new = []
    for entry in data:
        te_id = entry["ID"]
        if te_id not in remove_ids:
            data_new.append(entry)

    ## step six: clean up
    # os.remove(bed_report_tmp)
    # os.remove(bed_report_sort)
    # os.remove(bed_report_merge)
    # os.remove(json_report_tmp)

    json_report = out + "/" + "liftover_report.json"
    with open(json_report, "w") as output:
        json.dump(data_new, output, indent=4, sort_keys=False)

    # write non-reference lifted annotations in BED format
    bed_report = out + "/" + "liftover_nonref.bed"
    with open(bed_report, "w") as output:
        for item in data_new:
            if item["num_hits"] == 1:
                info = item["report"]
                chrom = info["chrom"]
                start = info["start"]
                end = info["end"]
                family = info["family"]
                score = "."
                strand = info["strand"]
                out_line = "\t".join(
                    [chrom, str(start), str(end), family, score, strand]
                )
                output.write(out_line + "\n")

    # generate summary
    nonref_total = 0
    ref_total = 0
    unlift_total = 0
    nonref_comments = dict()
    ref_comments = dict()
    unlift_comments = dict()
    for item in data_new:
        info = item["report"]
        if info["type"] == "non-reference":
            nonref_total = nonref_total + 1
            if "comment" in info.keys():
                if info["comment"] not in nonref_comments.keys():
                    nonref_comments[info["comment"]] = 1
                else:
                    nonref_comments[info["comment"]] = (
                        nonref_comments[info["comment"]] + 1
                    )
        elif info["type"] == "reference":
            ref_total = ref_total + 1
            if "comment" in info.keys():
                if info["comment"] not in ref_comments.keys():
                    ref_comments[info["comment"]] = 1
                else:
                    ref_comments[info["comment"]] = ref_comments[info["comment"]] + 1
        elif info["type"] == "unlifted":
            unlift_total = unlift_total + 1
            if "comment" in info.keys():
                if info["comment"] not in unlift_comments.keys():
                    unlift_comments[info["comment"]] = 1
                else:
                    unlift_comments[info["comment"]] = (
                        unlift_comments[info["comment"]] + 1
                    )
    summary_data = {
        "non-reference": {"total": nonref_total, "comments": nonref_comments},
        "reference": {"total": ref_total, "comments": ref_comments},
        "unlifted": {"total": unlift_total, "comments": unlift_comments},
    }

    summary_report = out + "/" + "liftover_summary.json"
    with open(summary_report, "w") as output:
        json.dump(summary_data, output, indent=4, sort_keys=False)

    # clean tmp files
    if not keep_files:
        print("Clean intermediate files...")
        shutil.rmtree(tmp_dir)
        shutil.rmtree(input_json_dir)

    print("Liftover finished!")
    return json_report


def main():
    args = get_args(program_version=__version__)

    # TODO: create soft link for fa1 and fa2?
    json_report = liftover(
        fasta1=args.fasta1,
        fasta2=args.fasta2,
        bed1=args.bed1,
        bed2=args.bed2,
        preset=args.preset,
        flank_len=args.flank_len,
        flank_gap_max=args.flank_gap_max,
        flank_overlap_max=args.flank_overlap_max,
        out=args.out,
        threads=args.threads,
        keep_files=args.keep_files,
        # single_flank=args.single_flank,
        different_contig_name=args.different_contig_name,
        telr_mode=args.telr_mode,
    )


if __name__ == "__main__":
    main()

import sys
import os
import json
import subprocess
import shutil
from multiprocessing import Pool
from STELR_utility.py import string_to_bool, get_subseq, check_exist

def extract_genome_seqs(ref, ref_size, te_json, flank_len, output_file):
    """
    Extract subsequence from fasta based on coordinates
    """
    with open(ref_size, "r") as ref_size_file:
        ref_size_dict = json.load(ref_size_file)
    ref_size = ref_size_dict[chrom]
    with open(te_json, "r") as te_json_file:
        te_dict = json.load(te_json_file)
        chrom = te_dict["chrom"]
        start = te_dict["start"]
        end = te_dict["end"]
    
    if end > ref_size or start < 0:
        pass#touch file in smk
    else:
        direction = {True:"5p",False:"3p"}["5p" in os.path.basename(output_file)]
        flank_len = int(flank_len)
        start, end = {"5p":start-flank_len+1,"3p":end}[direction], {"5p":start,"3p":end+flank_len}[direction]

        with open(output_file, "w") as output:
            output.write(f"{chrom}\t{start}\t{end}\n")
    #bedtools call moved to smk file


def paf_to_bed(paf, bed, contig_name, different_contig_name):
    """
    convert PAF to sorted BED format
    """
    if "3p" in os.path.basename(paf): different_contig_name = True
    else: different_contig_name = string_to_bool(different_contig_name)
    if different_contig_name:
        filter = "_".join(contig_name.split("_")[:-2])
    else:
        filter = None
    with open(paf, "r") as input, open(bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chrom = entry[5]
            if filter is None:
                start = entry[7]
                end = entry[8]
                name = entry[0]
                score = entry[11]
                strand = entry[4]
                out_line = "\t".join([chrom, start, end, name, score, strand])
                output.write(out_line + "\n")
            elif chrom == filter:
                start = entry[7]
                end = entry[8]
                name = entry[0]
                score = entry[11]
                strand = entry[4]
                out_line = "\t".join([chrom, start, end, name, score, strand])
                output.write(f"{out_line}\n")

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
    return int(start), int(end), gap


def create_bed(chrom, start, end, family, strand, filename):
    with open(filename, "w") as output:
        out_line = "\t".join([chrom, str(start), str(end), family, ".", strand])
        output.write(out_line)

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


def get_paf_info(paf, output_file):
    """
    Extract information from minimap2 PAF file
    """
    paf_info = dict()
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
    with open(output_file, "w") as output:
        json.dump(paf_info, output)


def get_genome_size(genome_index, size_file):
    genome_size_dict = dict()
    with open(genome_index, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            genome_size_dict[entry[0]] = int(entry[1])
    with open(size_file, "w") as output:
        json.dump(genome_size_dict, output)

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

def make_json(bed_input, json_output):
    with open(bed_input, "r") as input, open(json_output, "w") as output:
        entry = input[0].replace("\n", "").split("\t")
        chrom = entry[0]
        start = int(entry[1])
        end = int(entry[2])
        family = entry[3]
        strand = entry[5]

        annotation = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "family": family,
            "strand": strand#,
            #"fasta1": fasta1,
            #"fasta2": fasta2,
            #"out_dir": tmp_dir,
            #"flank_len": flank_len,
            #"flank_gap_max": flank_gap_max,
            #"flank_overlap_max": flank_overlap_max,
            #"bed2": bed2,
            #"preset": preset,
            # "single_flank": single_flank,
            #"different_contig_name": different_contig_name,
            #"telr_mode": telr_mode,
        }
        json.dump(annotation, output)

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
    #Already completed:
        # generate genome1 index file if not present
        genome1_index = fasta1 + ".fai"
        if not check_exist(genome1_index):
            build_index(fasta1)

        # generate genome2 index file if not present
        genome2_index = fasta2 + ".fai"
        if not check_exist(genome2_index):
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

            prefix = "_".join([chrom, str(start), str(end)])
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

def bed_to_json(overlap, info_5p, info_3p, json_file):
    with open(info_5p, "r") as info:
        align_5p_flank_qcs = json.load(info)
    with open(info_3p, "r") as info:
        align_3p_flank_qcs = json.load(info)
    overlap_dict = dict()
    with open(overlap, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")

            chrom_5p = entry[0]
            chrom_3p = entry[6]
            if (chrom_5p == chrom_3p and chrom_3p != "."):  # to make sure the entry exists
                align_5p_flank_id = "_".join([entry[3], entry[0], entry[1], entry[2]])
                align_3p_flank_id = "_".join([entry[9], entry[6], entry[7], entry[8]])

                start_5p = int(entry[1])
                end_5p = int(entry[2])
                start_3p = int(entry[7])
                end_3p = int(entry[8])

                id_chrom = chrom_5p.replace("(","").replace(")","").replace("_","").replace(",","").replace(":","")
                entry_id = f"({id_chrom}:{start_5p}-{end_5p},{start_3p}-{end_3p})"

                overlap_dict[entry_id] = {
                    "chrom_5p":chrom_5p,
                    "chrom_3p":chrom_3p,
                    "flank_strand":entry[5],
                    "mapp_quality_5p":int(entry[4]),
                    "mapp_quality_3p":int(entry[10]),
                    "align_5p_flank_qc":align_5p_flank_qcs[align_5p_flank_id],
                    "align_3p_flank_qc":align_3p_flank_qcs[align_3p_flank_id],
                    "start_5p":start_5p,
                    "end_5p":end_5p,
                    "start_3p":start_3p,
                    "end_3p":end_3p
                }
    if len(overlap_dict) > 0:
        with open(json_file, "w") as output:
            json.dump(overlap_dict, output)

def make_report(overlap_json, overlap_id, te_json, ref_bed, ref, flank_overlap_max, flank_gap_max, report_file):
    with open(te_json, "r") as te_json_file:
        te_dict = json.load(te_json_file)
        family = te_dict["family"]
        strand = te_dict["strand"]
        te_length = te_dict["end"] - te_dict["start"]
    reference_tes = True

    with open(overlap_json, "r") as overlap:
        overlap_dict = json.load(overlap)[overlap_id]
        chrom_genome2_5p = overlap_dict["chrom_5p"]
        chrom_genome2_3p = overlap_dict["chrom_3p"]
        flank_strand = overlap_dict["flank_strand"]
        mapp_quality_5p = overlap_dict["mapp_quality_5p"]
        mapp_quality_3p = overlap_dict["mapp_quality_3p"]
        align_5p_flank_qc = overlap_dict["align_5p_flank_qc"]
        align_3p_flank_qc = overlap_dict["align_3p_flank_qc"]
        start_5p = overlap_dict["start_5p"]
        end_5p = overlap_dict["end_5p"]
        start_3p = overlap_dict["start_3p"]
        end_3p = overlap_dict["end_3p"]

    lift_chrom = chrom_genome2_5p
    # get final liftover coordiantes
    lift_start, lift_end, lift_gap = get_coord(start_5p, end_5p, start_3p, end_3p, flank_strand)
    # figure out the strand of the lifted annotation
    if flank_strand == strand: lift_strand = "+"
    else: lift_strand = "-"
    # check the distance between flank alignment and nearast TE annotation
    distance_5p = check_nearby_ref(
        lift_chrom,
        start_5p,
        end_5p,
        family,
        lift_strand,
        ref_bed
    )
    distance_3p = check_nearby_ref(
        lift_chrom,
        start_3p,
        end_3p,
        family,
        lift_strand,
        ref_bed
    )
    # report the flanking sequence alignments
    align_5p_coord = f"{chrom_genome2_5p}:{start_5p}-{end_5p}"
    align_3p_coord = f"{chrom_genome2_3p}:{start_3p}-{end_3p}"
    lift_entry = {
        "type": None,
        "family": family,
        "chrom": lift_chrom,
        "start": lift_start,
        "end": lift_end,
        "strand": lift_strand,
        "gap": lift_gap,
        "TSD_length": None,
        "TSD_sequence": None,
        "5p_flank_align_coord": align_5p_coord,
        "5p_flank_mapping_quality": mapp_quality_5p,
        "5p_flank_num_residue_matches": align_5p_flank_qc["num_residue_matches"],
        "5p_flank_alignment_block_length": align_5p_flank_qc["alignment_block_length"],
        "5p_flank_sequence_identity": align_5p_flank_qc["sequence_identity"],
        "3p_flank_align_coord": align_3p_coord,
        "3p_flank_mapping_quality": mapp_quality_3p,
        "3p_flank_num_residue_matches": align_3p_flank_qc["num_residue_matches"],
        "3p_flank_alignment_block_length": align_3p_flank_qc["alignment_block_length"],
        "3p_flank_sequence_identity": align_3p_flank_qc["sequence_identity"],
        "distance_5p_flank_ref_te": distance_5p,
        "distance_3p_flank_ref_te": distance_3p,
        "comment": None,
    }

    # if the overlap between flank alignments is greater than 50bp, don't report
    reported = False
    if lift_gap < -flank_overlap_max:
        pass
    elif lift_gap >= -flank_overlap_max and lift_gap <= flank_gap_max:
        # if the gap between flank alignments is smaller than 50bp
        # report as reference if 1) there is a reference 2 TE in between two flanks or 2) the gap size and original TE size is similar or 3) the gap size is bigger than the size of the TE.
        # report as non-reference otherwise.
        if (
            is_reference(distance_5p, distance_3p, lift_gap)
            or check_nums_similar(lift_gap, te_length)
            or (lift_gap >= te_length)
        ):
            lift_entry["type"] = "reference"
            lift_entry["comment"] = "overlap/gap size between 3p and 5p flanks within threshold, include genome2 TE in between"
        else:
            lift_entry["type"] = "non-reference"
            lift_entry["comment"] = "overlap/gap size between 3p and 5p flanks within threshold"
            # get TSD length and sequence
            if lift_gap == 0:
                lift_entry["TSD_length"] = 0
                lift_entry["TSD_sequence"] = None
            if lift_gap < 0:
                lift_entry["TSD_length"] = -lift_gap
                lift_entry["TSD_sequence"] = get_subseq(ref, lift_chrom, lift_start, lift_end)
        reported = True
    else:  # if the gap between flanks is greater than 50bp
        if lift_gap > flank_gap_max and lift_gap <= 0.5 * te_length:
            # if the gap size is smaller than half of the original TE size
            # report as reference if 1) there is a reference 2 TE in between two flanks
            # report as non-reference otherwise.
            if is_reference(distance_5p, distance_3p, lift_gap):
                lift_entry["type"] = "reference"
                lift_entry["comment"] = "flanks gap size less than half of TE annotation, include genome2 TE in between"
            else:
                lift_entry["type"] = "non-reference"
                lift_entry["comment"] = "flanks gap size exceeds threshold but less than half of TE annotation, no genome2 TE in between"
                num_hits = num_hits + 1
            lift_entries["report"].append(lift_entry)
            reported = True
        elif lift_gap >= 0.5 * te_length and lift_gap <= 20000:
            # if the gap size is greater than half of the original TE size and smaller than 20kb
            lift_entry["type"] = "reference"
            if is_reference(distance_5p, distance_3p, lift_gap):
                lift_entry["comment"] = "flanks gap size greater than half of TE annotation, include genome2 TE in between"  # TODO: check same family?
            else:
                lift_entry["comment"] = "flanks gap size greater than half of TE annotation, no genome2 TE in between"
            lift_entries["report"].append(lift_entry)
            reported = True
    if reported:
        with open(report_file, "w") as output:
            json.dump(lift_entry, output)

def choose_report(*input_files):
    # more than one report per annotation, need to filter
    # first step is to choose the best reference annotation by comparing gap size with TE length
    reports = []
    flanks = {"5p":{"bed_file":None},"3p":{"bed_file":None}}
    te_dict = {}
    for file in input_files:
        if "_flank" in os.path.basename(file):
            flank = {True:"5p",False:"3p"}["5p" in os.path.basename(file)]
            for file_type in ["bed", "info"]
                if file_type in os.path.basename(file):      
                    flanks[flank][f"{file_type}_file"] = file
        elif "_best_report.json" in file:
            out_file = file
        elif check_exist(file):
            with open(file, "r") as input:
                if os.path.basename(file) == "00_annotation.json": 
                    te_dict = json.load(input)
                else: reports.append(json.load(input))
    strand = te_dict["strand"]
    best_report = {}
    reported = True
    if len(reports) > 1:
        best_ref_entry = dict()
        best_nonref_entry = dict()
        for report in reports:
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
        if reported:
            if best_ref_entry and best_nonref_entry:
                # if both reference and non-reference liftover can be found, report the non-reference one
                # lift_entries["report"].append(best_ref_entry)
                best_report = best_nonref_entry
            elif best_ref_entry:
                best_report = best_ref_entry
            elif best_nonref_entry:
                best_report = best_nonref_entry
            else:
                reported = False
    elif len(lift_entries["report"]) == 1:
        lift_entries["report"] = lift_entries["report"][0]
    else: reported = False

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
        for flank in flanks:
            flanks[flank]["alignment_coords"] = []
            if check_exist(flanks[flank]["bed_file"]):
                with open(flanks[flank]["bed_file"], "r") as input:
                    for line in input:
                        entry = line.replace("\n", "").split("\t")
                        flanks[flank]["alignment_coords"].append(entry[0] + ":" + entry[1] + "-" + entry[2])

            if len(flanks[flank]["alignment_coords"]) == 1:
                lift_entry[f"{flank}_flank_align_coord"] = flanks[flank]["alignment_coords"][0]
            elif len(flanks[flank]["alignment_coords"]) > 1:
                lift_entry[f"{flank}_flank_align_coord"] = flanks[flank]["alignment_coords"]

        # if single flank mode is turned on, inspect single flanks and report as lifted
        # if single_flank:
        if len(flanks["5p"]["alignment_coords"]) == 1 and len(flanks["3p"]["alignment_coords"]) == 0:
            lift_entry = single_flank_liftover(flanks, "5p", lift_entry, strand)
        elif len(flanks["5p"]["alignment_coords"]) == 0 and len(flanks["3p"]["alignment_coords"]) == 1:
            lift_entry = single_flank_liftover(flanks, "3p", lift_entry, strand)
        

    # write
    with open(out_file, "w") as output:
        json.dump(liftover_report, output)

def is_reference(distance_5p, distance_3p, lift_gap):
    return (
        distance_5p is not None
        and distance_5p >= 0
        and distance_5p <= lift_gap
        and distance_3p is not None
        and distance_3p <= 0
        and -distance_3p <= lift_gap
    )

def single_flank_liftover(flank_info, flank, lift_entry, strand, ref_bed):
    flank_info = flank_info[flank]

    with open(flank_info["info_file"], "r") as input:
        align_flank_qcs = json.load(input)

    with open(flank_info["bed_file"], "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            flank_chrom = entry[0]
            flank_start = int(entry[1])
            flank_end = int(entry[2])
            flank_mapping_quality = int(entry[4])
            flank_strand = entry[5]
            align_flank_id = "_".join(
                [entry[3], entry[0], entry[1], entry[2]]
            )
    if flank_strand == strand:
        lift_strand = "+"
    else:
        lift_strand = "-"
    if flank_strand == "+":
        lift_start = lift_end = {"5p":flank_end,"3p":flank_start}[flank]
    else:
        lift_start = lift_end = {"5p":flank_start,"3p":flank_end}[flank]

    align_flank_qc = align_flank_qcs[align_flank_id]

    lift_entry["chrom"] = flank_chrom
    lift_entry["start"] = int(lift_start)
    lift_entry["end"] = int(lift_end)
    lift_entry["mapp_quality_5p"] = flank_mapping_quality
    lift_entry["strand"] = lift_strand
    lift_entry["5p_flank_num_residue_matches"] = align_flank_qc["num_residue_matches"]
    lift_entry["5p_flank_alignment_block_length"] = align_flank_qc["alignment_block_length"]
    lift_entry["5p_flank_sequence_identity"] = align_flank_qc["sequence_identity"]

    distance = check_nearby_ref(
        flank_chrom,
        flank_start,
        flank_end,
        family,
        lift_strand,
        ref_bed
    )
    lift_entry[f"distance_{flank}_flank_ref_te"] = distance
    # if distance between flank and ref is small, report as ref, otherwise as non-ref
    if distance is not None and abs(distance) <= 5:
        lift_entry["type"] = "reference"
        lift_entry[
            "comment"
        ] = "only one flank aligned, flank alignment adjacent to reference TE"
    else:
        lift_entry["type"] = "non-reference"
        lift_entry[
            "comment"
        ] = "only one flank aligned, flank alignment not adjacent to reference TE"
    
    return lift_entry

check_nearby_ref(chrom, start_query, end_query, family, strand, ref_bed, threshold=5000):
    """
    Check if flanking seqeunce alignments have nearby TE annotations in genome 2
    """
    distance = None

    if check_exist(ref_bed):
        bed = f"{chrom}\t{start}\t{end}\t{family}\t.\t{strand}".replace("'","")
        overlap = subprocess.call(f"echo '{bed}' | bedtools closest -a - -b {ref_bed} -d -D ref -k 5", shell=True, capture_output=True, text=True).stdout

        for line in overlap.split("\n"):
            entry = line.split("\t")
            chrom2 = entry[6]
            family2 = entry[9]
            strand2 = entry[11]
            if chrom == chrom2 and family == family2 and strand == strand2:
                distance_new = int(entry[12])
                if distance is None:
                    distance = distance_new
                else:
                    distance = absmin(distance, distance_new)
    if distance is not None:
        if abs(distance) > threshold:
            distance = None
    return distance

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
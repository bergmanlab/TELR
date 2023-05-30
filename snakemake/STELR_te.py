import os
import sys
import subprocess
import json
from Bio import SeqIO
import re
import logging
import time
import statistics
from multiprocessing import Pool
from TELR_utility import (
    check_exist,
    mkdir,
    format_time,
    get_cmd_output,
    get_rev_comp_sequence,
    get_contig_name,
    get_contig_length,
    read_vcf
)
#from telr.TELR_liftover import liftover
#from telr.TELR_assembly import prep_assembly_inputs

def get_vcf_seq(contig, vcf_parsed, sequence_file):
    sequence = read_vcf(vcf_parsed, contig, column=7, single_contig=True)
    with open(sequence_file, "w") as output:
        output.write(f">{contig}\n")
        output.write(sequence)

def vcf_alignment_filter(intersect_file, output_file):
    with open(intersect_file, "r") as input, open(output_file, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            # the overlap between VCF sequence alignment and TE-contig alignment has to be over 10bp
            if int(entry[12]) > 10:
                out_line = "\t".join(
                    [entry[0], entry[1], entry[2], entry[3], entry[4], entry[5]]
                )
                output.write(out_line + "\n")

def annotate_contig(merge_output, annotation_file):
    te_dir = annotation_file[:annotation_file.rindex("/")]
    mkdir(te_dir)
    with open(merged_output, "r") as input, open(annotation_file, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            contig_te_start = entry[1]
            contig_te_end = entry[2]
            contig_te_family = entry[3]
            contig_te_strand = entry[4]
            if contig_te_strand != "+" and contig_te_strand != "-":
                contig_te_strand = "."
            out_line = "\t".join(
                [
                    contig_name,
                    contig_te_start,
                    contig_te_end,
                    contig_te_family,
                    ".",
                    contig_te_strand
                ]
            )
            output.write(out_line + "\n")
            te_name = f"te_{contig_te_start}_{contig_te_end}"
            mkdir(f"{te_dir}/{te_name}")
            with open(f"{te_dir}/{te_name}/00_annotation.bed", "w") as te_annotation:
                te_annotation.write(out_line)

def rm_parse_merge(input_file, output_file):
    with open(input_file, "r") as input, open(output_file, "w") as output:
        for line in input:
            if "##" not in line:
                entry = line.replace("\n", "").split("\t")
                contig_name = entry[0].rsplit(":", 1)[0]
                start = entry[0].rsplit(":", 1)[1].split("-")[0]
                end = entry[0].rsplit(":", 1)[1].split("-")[1]
                # contigs = entry[0].replace(':', '-').split("-")
                family = re.sub('Target "Motif:|".*', "", entry[8])
                strand = entry[6]
                score = entry[5]
                out_line = "\t".join(
                    [contig_name, start, end, family, score, strand]
                )
                output.write(out_line + "\n")

def rm_reannotate(rm_raw, original_bed, output_file):
    # replace contig_te_annotation family with ones from RM
    contig_te_annotation_new = contig_te_annotation_sorted.replace(
        "bed", "family_reannotated.bed"
    )
    contig_rm_family_dict = dict()
    with open(rm_raw, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            family = entry[3]
            contig_rm_family_dict[contig_name] = family

    with open(output_file, "w") as output, open(original_bed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            contig_te_start = entry[1]
            contig_te_end = entry[2]
            if contig_name in contig_rm_family_dict:
                contig_te_family = contig_rm_family_dict[contig_name]
                contig_te_strand = entry[5]
                out_line = "\t".join(
                    [
                        contig_name,
                        contig_te_start,
                        contig_te_end,
                        contig_te_family,
                        ".",
                        contig_te_strand,
                    ]
                )
                output.write(out_line + "\n")

#delete this?
def seq2contig(seq, contig, out):
    with open(out, "a") as output:
        subprocess.call(
            ["minimap2", "-cx", "map-pb", "--secondary=no", contig, seq], stdout=output
        )  # only retain primary alignment

def gff3tobed(gff, bed):
    # check GFF3 format
    repeats_in_ref = True
    with open(gff, "r") as input:
        for line in input:
            if "No repetitive sequences detected" in line:
                repeats_in_ref = False
            elif "#" not in line:
                if "Target=" not in line:
                    print(
                        "Incorrect GFF3 format, please check README for expected format, exiting..."
                    )
                    logging.exception(
                        "Incorrect GFF3 format, please check README for expected format, exiting..."
                    )
                    sys.exit(1)
                break
    if repeats_in_ref:
        with open(bed, "w") as output, open(gff, "r") as input:
            for line in input:
                if "#" not in line:
                    entry = line.replace("\n", "").split("\t")
                    info = entry[8].split(";")
                    for item in info:
                        if "Target=" in item:
                            family = item.replace("Target=", "")
                    out_line = "\t".join(
                        [entry[0], str(int(entry[3]) - 1), entry[4], family, ".", entry[6]]
                    )
                    output.write(out_line + "\n")

def parse_rm_out(rm_gff, rm_out, gff3):
    if os.stat(rm_gff).st_size != 0:
        with open(gff3, "w") as output, open(rm_gff, "r") as input:
            for line in input:
                if "RepeatMasker" in line:
                    entry = line.replace("\n", "").split("\t")
                    family = entry[8].split(" ")[1]
                    family = re.sub('"Motif:', "", family)
                    family = re.sub('"', "", family)
                    out_line = "\t".join(
                        [
                            entry[0],
                            "RepeatMasker",
                            "dispersed_repeat",
                            entry[3],
                            entry[4],
                            entry[5],
                            entry[6],
                            entry[7],
                            f"Target={family}",
                        ]
                    )
                    output.write(f"{out_line}\n")
    else:
        with open(gff3, "w") as output, open(rm_out, "r") as input:
            input = "".join([line for line in input])
            if "There were no repetitive sequences detected" in input:
                output.write("No repetitive sequences detected")
            else:
                raise Exception("Repeatmasking failed, exiting...")

def get_te_flank_ratio(cov_dict):
    if cov_dict["te"] and cov_dict["flank"]:
        if cov_dict["flank"] == 0:
            return None
        else:
            ratio = cov_dict["te"] / cov_dict["flank"]
            if ratio > 1.5:
                return None
            else:
                return ratio
    else:
        return None

def estimate_coverage(te_5p, te_3p, flank_5p, flank_3p, frequency_file):
    # read contig annotation to dict
    depths = {
        "5p":{
            "te":median_cov(te_5p),
            "flank":median_cov(flank_5p)
        }
        "flank":{
            "te":median_cov(te_3p),
            "flank":median_cov(flank_3p)
        }
    }
    if depths["3p"]["te"] is None:
        depths["3p"]["te"] = depths["5p"]["te"]
    with open(frequency_file, "w") as output:
        json.dump(depths, output)

def get_af(fwd_freq, rev_freq, out_file):
    allele_freqs = {
        "fwd": fwd_freq,
        "rev": rev_freq
    }
    for direction in allele_freqs:
        with open(allele_freqs[direction], "r") as input:
            allele_freqs[direction] = json.load(input)
    
    taf_5p = get_te_flank_ratio(allele_freqs["fwd"]["5p"])
    taf_3p = get_te_flank_ratio(allele_freqs["rev"]["5p"])
    if taf_5p and taf_3p:
        if abs(taf_5p - taf_3p) <= 0.3:
            freq = (taf_5p + taf_3p) / 2
        else:
            freq = None
    elif taf_5p:
        freq = taf_5p
    elif taf_3p:
        freq = taf_3p
    else:
        freq = None
    if freq:
        if freq > 1:
            freq = 1
        allele_freqs["freq"] = round(freq, 3)
    else: allele_freqs["freq"] = None

    with open(outfile, "w") as output:
        json.dump(allele_freqs, output)

def get_start_end(te, contig):
    start = int(te.split("_")[-2])
    end = int(te.split("_")[-1])
    if "revcomp" in contig:
        length = get_contig_length(contig)
        start, end = length - end, length - start
    return start, end    

def estimate_te_depth(bam, contig, te, te_interval_size, te_offset, output_5p, output_3p):
    start, end = get_start_end(te, contig)
    contig_name = get_contig_name(contig)
    te_interval_size, te_offset = int(te_interval_size), int(te_offset)
    commands = {
        "5p":f"samtools depth -aa -r {contig_name}:{start}-{end} {bam} > {output_5p}",
        "3p":f"touch {output_3p}"
    }
    if te_interval_size:
        if start + te_offset + te_interval_size < end:
            commands["5p"] = f"samtools depth -aa -r {contig_name}:{start+te_offset}-{end+te_offset+te_interval_size} {bam} > {output_5p}"
            commands["3p"] = f"samtools depth -aa -r {contig_name}:{end-te_interval_size-te_offset}-{end-te_offset} {bam} > {output_3p}"
    
    subprocess.run(commands["5p"], shell=True)
    subprocess.run(commands["3p"], shell=True)

def estimate_flank_depth(bam, contig, te, flank_len, flank_offset, output_5p, output_3p):
    start, end = get_start_end(te, contig)
    contig_length = get_contig_length(contig)
    contig_name = get_contig_name(contig)
    flank_len, flank_offset = int(flank_len), int(flank_offset)
    commands = {
        "5p":f"touch {output_5p}",
        "3p":f"touch {output_3p}"
    }
    if start - flank_len - flank_offset >= 0:
        commands["5p"] = f"samtools depth -aa -r {contig_name}:{start-flank_len-offset}-{start-offset} {bam} > {output_5p}"
    if end + flank_len + flank_offset <= contig_length:
        commands["3p"] = f"samtools depth -aa -r {contig_name}:{end+flank_offset}-{end+flank_len+flank_offset} {bam} > {output_3p}"
    
    subprocess.run(commands["5p"], shell=True)
    subprocess.run(commands["3p"], shell=True)

def get_median_cov(depth):
    covs = []
    if os.stat(depth).st_size != 0:
        with open(depth, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                covs.append(int(entry[2]))
        median_cov = statistics.median(covs)
    else: median_cov = None
    return median_cov

def find_te(
    contigs_fa,
    reference,
    contig_te_bed,
    ref_te_bed,
    out,
    gap,
    overlap,
    flank_len,
    # single_flank,
    different_contig_name,
    keep_files,
    thread,
):
    """
    Identify non-reference TE insertions in the reference genome using assembled contigs
    """
    # default parameters
    presets = "asm10"

    # lift over
    logging.info("Map contigs to reference...")

    json_report = liftover(
        fasta1=contigs_fa,
        fasta2=reference,
        bed1=contig_te_bed,
        bed2=ref_te_bed,
        preset=presets,
        flank_len=flank_len,
        flank_gap_max=gap,
        flank_overlap_max=overlap,
        out=out,
        threads=thread,
        keep_files=keep_files,
        # single_flank=single_flank,
        different_contig_name=different_contig_name,
        telr_mode=True,
    )

    return json_report

def create_fa(header, seq, out):
    with open(out, "w") as output:
        output.write(">" + header + "\n")
        output.write(seq)

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
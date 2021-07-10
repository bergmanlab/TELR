from TELR_utility import check_exist
import os
import pandas as pd
import logging
import json
from Bio import SeqIO
from datetime import date
import subprocess


def generate_output(
    liftover_report_path,
    te_freq_dict,
    te_fa,
    vcf_parsed,
    contig_te_annotation,
    contig_fa,
    out,
    sample_name,
    ref,
):
    logging.info("Write output...")

    # load merged contigs
    contig_length_dict = dict()
    with open(contig_fa, "r") as handle:
        contig_records = SeqIO.parse(handle, "fasta")
        for record in contig_records:
            contig_name = record.id
            contig_desc = record.description.split(" ")[1]
            contig_length = int(contig_desc.replace("len=", ""))
            contig_length_dict[contig_name] = contig_length

    # load contig te annotation
    contig_te_strand_dict = dict()
    with open(contig_te_annotation, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            te_strand = entry[5]
            if te_strand != "+" and te_strand != "-":
                contig_te_strand_dict[contig_name] = "."
            else:
                contig_te_strand_dict[contig_name] = te_strand

    # load info from VCF file
    sniffles_info = dict()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = (
                line.replace("\n", "").replace(" ", "").split("\t")
            )  # TODO: some how there are spaces in the parsed vcf file
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            gt = entry[10]
            ref_count = entry[11]
            alt_count = entry[12]
            sniffles_info[contig_name] = {
                "gt": gt,
                "alt_count": alt_count,
                "ref_count": ref_count,
            }

    # load TE sequence fasta
    te_seqs_dict = dict()
    with open(te_fa, "r") as input:
        for record in SeqIO.parse(input, "fasta"):
            ins_name = record.id
            te_seqs_dict[ins_name] = record.seq

    # load liftover report
    with open(liftover_report_path) as f:
        liftover_report = json.load(f)

    # initiate final report
    final_report = []
    final_report_expanded = []

    contig_ids = set()
    for item in liftover_report:
        insertion_report = {
            "type": None,
            "ID": None,
            "chrom": None,
            "start": None,
            "end": None,
            "family": None,
            "strand": None,
            "support": None,
            "tsd_length": None,
            "tsd_sequence": None,
            "te_sequence": None,
            "gt": None,
            "num_sv_reads": None,
            "num_ref_reads": None,
            "allele_frequency": None,
        }

        insertion_report_expanded = {
            "type": None,
            "ID": None,
            "chrom": None,
            "start": None,
            "end": None,
            "family": None,
            "strand": None,
            "support": None,
            "tsd_length": None,
            "tsd_sequence": None,
            "te_sequence": None,
            "gt": None,
            "num_sv_reads": None,
            "num_ref_reads": None,
            "allele_frequency": None,
            "gap_between_flank": None,
            "te_length": None,
            "contig_id": None,
            "contig_length": None,
            "contig_te_start": None,
            "contig_te_end": None,
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
        }
        report_info = item["report"]
        ins_name = item["genome1_coord"]
        if report_info["type"] == "non-reference":
            # fill simple record
            insertion_report["type"] = report_info["type"]
            insertion_report["chrom"] = report_info["chrom"]
            insertion_report["start"] = report_info["start"]
            insertion_report["end"] = report_info["end"]
            insertion_report["family"] = report_info["family"]
            insertion_report["ID"] = "_".join(
                [
                    report_info["chrom"],
                    str(report_info["start"]),
                    str(report_info["end"]),
                    report_info["family"],
                ]
            )
            insertion_report["strand"] = report_info["strand"]
            insertion_report["tsd_length"] = report_info["TSD_length"]
            insertion_report["tsd_sequence"] = report_info["TSD_sequence"]

            contig_id = ins_name.split(":")[0]
            contig_ids.add(contig_id)
            contig_te_strand = contig_te_strand_dict[contig_id]

            if contig_te_strand == "+" or contig_te_strand == ".":
                insertion_report["te_sequence"] = str(te_seqs_dict[ins_name])
            else:
                insertion_report["te_sequence"] = str(
                    te_seqs_dict[ins_name].reverse_complement()
                )

            insertion_report["gt"] = sniffles_info[contig_id]["gt"]
            insertion_report["num_sv_reads"] = sniffles_info[contig_id]["alt_count"]
            insertion_report["num_ref_reads"] = sniffles_info[contig_id]["ref_count"]
            insertion_report["allele_frequency"] = te_freq_dict[contig_id]

            # fill expaned report
            insertion_report_expanded["contig_length"] = contig_length_dict[contig_id]
            insertion_report_expanded["gap_between_flank"] = report_info["gap"]
            if report_info["gap"] is not None:
                insertion_report["support"] = "single_side"
            else:
                insertion_report["support"] = "both_sides"

            insertion_report_expanded["te_length"] = item["te_length"]

            insertion_report_expanded["contig_id"] = contig_id
            insertion_report_expanded["te_length"] = len(
                insertion_report["te_sequence"]
            )

            contig_te_coord = ins_name.split(":")[1]
            insertion_report_expanded["contig_te_start"] = int(
                contig_te_coord.split("-")[0]
            )
            insertion_report_expanded["contig_te_end"] = int(
                contig_te_coord.split("-")[1]
            )

            # fill expanded record
            insertion_report_expanded["5p_flank_align_coord"] = report_info[
                "5p_flank_align_coord"
            ]
            insertion_report_expanded["5p_flank_mapping_quality"] = report_info[
                "5p_flank_mapping_quality"
            ]
            insertion_report_expanded["5p_flank_num_residue_matches"] = report_info[
                "5p_flank_num_residue_matches"
            ]
            insertion_report_expanded["5p_flank_alignment_block_length"] = report_info[
                "5p_flank_alignment_block_length"
            ]
            insertion_report_expanded["5p_flank_sequence_identity"] = report_info[
                "5p_flank_sequence_identity"
            ]
            insertion_report_expanded["3p_flank_align_coord"] = report_info[
                "3p_flank_align_coord"
            ]
            insertion_report_expanded["3p_flank_mapping_quality"] = report_info[
                "3p_flank_mapping_quality"
            ]
            insertion_report_expanded["3p_flank_num_residue_matches"] = report_info[
                "3p_flank_num_residue_matches"
            ]
            insertion_report_expanded["3p_flank_alignment_block_length"] = report_info[
                "3p_flank_alignment_block_length"
            ]
            insertion_report_expanded["3p_flank_sequence_identity"] = report_info[
                "3p_flank_sequence_identity"
            ]

            final_report.append(insertion_report)
            # print(final_report)
            insertion_report_expanded.update(insertion_report)
            final_report_expanded.append(insertion_report_expanded)
            # print(final_report)

    # write in JSON format
    insertion_report_path = out + "/" + sample_name + ".telr.json"
    with open(insertion_report_path, "w") as output:
        json.dump(final_report, output, indent=4, sort_keys=False)

    insertion_report_expanded_path = out + "/" + sample_name + ".telr.expanded.json"
    with open(insertion_report_expanded_path, "w") as output:
        json.dump(final_report_expanded, output, indent=4, sort_keys=False)

    # write TE sequences in fasta
    te_fa_path = out + "/" + sample_name + ".telr.te.fasta"
    with open(te_fa_path, "w") as output:
        for item in final_report:
            record.id = (
                item["chrom"]
                + "_"
                + str(item["start"])
                + "_"
                + str(item["end"])
                + "#"
                + item["family"]
            )
            output.write(">" + record.id + "\n" + item["te_sequence"] + "\n")

    # write contig sequences in fasta
    contig_fa_path = out + "/" + sample_name + ".telr.contig.fasta"
    with open(contig_fa_path, "w") as output, open(contig_fa, "r") as handle:
        contig_records = SeqIO.parse(handle, "fasta")
        for record in contig_records:
            contig_name = record.id
            if contig_name in contig_ids:
                SeqIO.write(record, output, "fasta")

    # write in VCF format
    vcf_out = os.path.join(out, sample_name + ".telr.vcf")
    ref_info = get_contig_info(ref)
    write_vcf(final_report, ref, ref_info, vcf_out)

    # write in BED format
    bed_out = os.path.join(out, sample_name + ".telr.bed")
    write_bed(final_report, bed_out)


def write_bed(final_report, bed):
    with open(bed, "w") as output:
        for item in final_report:
            chrom = item["chrom"]
            start = item["start"]
            end = item["end"]
            family = item["family"]
            score = "."
            strand = item["strand"]
            out_line = "\t".join([chrom, str(start), str(end), family, score, strand])
            output.write(out_line + "\n")


def write_vcf(input, ref, ref_info, out_vcf):
    df = pd.DataFrame(input)
    if not df.empty:
        df["ID"] = df.index
        df["start"] = df["start"] + 1
        df["REF"] = "N"
        df["QUAL"] = "."
        df["FILTER"] = "PASS"
        df["FORMAT"] = "GT:DR:DV"
        df["gt"] = df["gt"] + ":" + df["num_sv_reads"] + ":" + df["num_ref_reads"]
        df["INFO"] = df.apply(
            lambda x: "SVTYPE=INS"
            + ";END="
            + str(x.end)
            + ";FAMILY="
            + str(x.family)
            + ";STRANDS="
            + str(x.strand)
            + ";SUPPORT_TYPE="
            + str(x.support)
            + ";RE="
            + str(x.num_sv_reads)
            + ";AF="
            + str(x.allele_frequency)
            + ";TSD_LEN="
            + str(x.tsd_length)
            + ";TSD_SEQ="
            + str(x.tsd_sequence),
            axis=1,
        )

        df = df[
            [
                "chrom",
                "start",
                "ID",
                "REF",
                "te_sequence",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "gt",
            ]
        ]
    with open(out_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.1" + "\n")
        vcf.write("##fileDate={}".format(date.today()) + "\n")
        vcf.write("##source=TELR" + "\n")
        vcf.write("##reference=" + ref + "\n")
        vcf.write("\n".join(ref_info) + "\n")
        vcf.write(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structure variant">'
            + "\n"
        )
        vcf.write(
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structure variant">'
            + "\n"
        )
        vcf.write(
            '##INFO=<ID=STRANDS,Number=A,Type=String,Description="Strand orientation">'
            + "\n"
        )
        vcf.write(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">' + "\n"
        )
        vcf.write(
            '##INFO=<ID=FAMILY,Number=1,Type=String,Description="TE family">' + "\n"
        )
        vcf.write(
            '##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">' + "\n"
        )
        vcf.write(
            '##INFO=<ID=SUPPORT_TYPE,Number=1,Type=String,Description="Type of TE flank alignment to reference genome (1: single flank; 2: two flanks with gap; 3: two flanks with overlap)">'
            + "\n"
        )
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
        vcf.write(
            '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">'
            + "\n"
        )
        vcf.write(
            '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">'
            + "\n"
        )
        vcf.write(
            "#"
            + "\t".join(
                [
                    "CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "SAMPLE",
                ]
            )
            + "\n"
        )
    if not df.empty:
        df.to_csv(out_vcf, sep="\t", mode="a", index=False, header=False)


def get_contig_info(reference):
    with open(reference, "r") as input:
        subprocess.call(["samtools", "faidx", reference])
        index_file = reference + ".fai"
    contig_info = []
    with open(index_file, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_info.append("##contig=<ID={},length={}>".format(entry[0], entry[1]))
    return contig_info

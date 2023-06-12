import os
import sys
import pandas as pd
import logging
import json
from Bio import SeqIO
from datetime import date
import subprocess
from STELR_utility import check_exist

def write_vcf(input, ref, ref_index, out_vcf):
    ref_info = get_contig_info(ref_index)
    df = pd.DataFrame(input)
    if not df.empty:
        df["ID"] = df.index
        df["start"] = df["start"] + 1
        df["REF"] = "N"
        df["QUAL"] = "."
        df["FILTER"] = "PASS"
        df["FORMAT"] = "GT:DR:DV"
        df["gt"] = df["genotype"] + ":" + df["num_sv_reads"] + ":" + df["num_ref_reads"]
        df["INFO"] = df.apply(
            lambda x: f"SVTYPE=INS;END={x.end};FAMILY={x.family};STRANDS={x.strand};SUPPORT_TYPE={x.support};RE={x.num_sv_reads};AF={x.allele_frequency};TSD_LEN={x.tsd_length};TSD_SEQ={x.tsd_sequence}",
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
        df = df.fillna("NA")
    with open(out_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.1\n")
        vcf.write("##fileDate={}".format(date.today()) + "\n")
        vcf.write("##source=TELR\n")
        vcf.write(f"##reference={ref}\n")
        vcf.write("\n".join(ref_info) + "\n")
        vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structure variant">\n')
        vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structure variant">\n')
        vcf.write('##INFO=<ID=STRANDS,Number=A,Type=String,Description="Strand orientation">\n')
        vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        vcf.write('##INFO=<ID=FAMILY,Number=1,Type=String,Description="TE family">\n')
        vcf.write('##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">\n')
        vcf.write('##INFO=<ID=SUPPORT_TYPE,Number=1,Type=String,Description="single_side or both_sides">\n')
        vcf.write('##INFO=<ID=TSD_LEN,Number=1,Type=String,Description="Length of the TSD sequence if available">\n')
        vcf.write('##INFO=<ID=TSD_SEQ,Number=1,Type=String,Description="TSD sequence if available">\n')
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">\n')
        vcf.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        
    if not df.empty:
        df.to_csv(out_vcf, sep="\t", mode="a", index=False, header=False)


def get_contig_info(reference_index):
    contig_info = []
    with open(reference_index, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_info.append("##contig=<ID={},length={}>".format(entry[0], entry[1]))
    return contig_info

def make_json_output(liftover_file, af_file, vcf_parsed_file, annotation_file, contig_file, json_output):
    with open(liftover_file, "r") as data:
        liftover_report = json.load(data)
    if not liftover_report["type"] == "non-reference":
        quit()
    with open(af_file, "r") as data:
        af_dict = json.load(data)
    with open(vcf_parsed_file, "r") as data:
        sniffles_info = [line for line in data][0].replace("\n","").split("\t")
        sniffles_info = {
            "genotype":sniffles_info[10],
            "alt_count":sniffles_info[12],
            "ref_count":sniffles_info[11]
        }
    sequence = subprocess.run(f"bedtools getfasta -fi '{contig_file}' -bed '{annotation_file}' -s", shell=True, capture_output=True, text=True).stdout.split("\n")[1]
    contig_name = liftover_file.split("contigs/")[1].split("/tes/")[0]
    te_name = liftover_file.split("/tes/")[1].split("/")[0]
    te_start, te_end = [int(index) for index in te_name.replace("te_","").split("_")]
    contig_length = 0
    with open(contig_file, "r") as data:
        next(data)
        for line in data:
            contig_length += len(line.replace("\n",""))

    full_report = {
        "type": "non-reference",
        "ID": f"{liftover_report['chrom']}_{liftover_report['start']}_{liftover_report['end']}_{liftover_report['family']}",
        "chrom": liftover_report["chrom"],
        "start": liftover_report["start"],
        "end": liftover_report["end"],
        "family": liftover_report["family"],
        "strand": liftover_report["strand"],
        "support": "single_side",#will be processed further
        "tsd_length": liftover_report["TSD_length"],
        "tsd_sequence": liftover_report["TSD_sequence"],#will be processed further
        "te_sequence": sequence,
        "genotype": sniffles_info["genotype"],
        "num_sv_reads": sniffles_info["alt_count"],
        "num_ref_reads": sniffles_info["ref_count"],
        "allele_frequency": af_dict["freq"],
        "gap_between_flank": liftover_report["gap"],
        "te_length": len(sequence),
        "contig_id": contig_name,
        "contig_length": contig_length,
        "contig_te_start": te_start,
        "contig_te_end": te_end,
        "5p_flank_align_coord": liftover_report["5p_flank_align_coord"],
        "5p_flank_mapping_quality": liftover_report["5p_flank_mapping_quality"],
        "5p_flank_num_residue_matches": liftover_report["5p_flank_num_residue_matches"],
        "5p_flank_alignment_block_length": liftover_report["5p_flank_alignment_block_length"],
        "5p_flank_sequence_identity": liftover_report["5p_flank_sequence_identity"],
        "3p_flank_align_coord": liftover_report["3p_flank_align_coord"],
        "3p_flank_mapping_quality": liftover_report["3p_flank_mapping_quality"],
        "3p_flank_num_residue_matches": liftover_report["3p_flank_num_residue_matches"],
        "3p_flank_alignment_block_length": liftover_report["3p_flank_alignment_block_length"],
        "3p_flank_sequence_identity": liftover_report["3p_flank_sequence_identity"],
        "te_5p_cov": af_dict["fwd"]["5p"]["te"],
        "te_3p_cov": af_dict["fwd"]["3p"]["te"],
        "flank_5p_cov": af_dict["fwd"]["5p"]["flank"],
        "flank_3p_cov": af_dict["fwd"]["3p"]["flank"],
        "te_5p_cov_rc": af_dict["rev"]["5p"]["te"],
        "te_3p_cov_rc": af_dict["rev"]["3p"]["te"],
        "flank_5p_cov_rc": af_dict["rev"]["5p"]["flank"],
        "flank_3p_cov_rc": af_dict["rev"]["3p"]["flank"]
    }

    if full_report["tsd_sequence"]:
        full_report["tsd_sequence"] = full_report["tsd_sequence"].upper()
    
    if full_report["5p_flank_align_coord"] and full_report["3p_flank_align_coord"]:
        full_report["support"] = "both_sides"
    
    basic_report = ["type","ID","chrom","start","end","family","strand","support","tsd_length","tsd_sequence","te_sequence","genotype","num_sv_reads","num_ref_reads","allele_frequency"]
    basic_report = {key:full_report[key] for key in basic_report}

    output_information = {
        "expanded_json":full_report,
        "json":basic_report,
        "contig_path":contig_file,
        "te_fasta":f">{full_report['ID']}\n{full_report['te_sequence']}\n",
        "contig_name":contig_name,
        "bed_out":f"{full_report['chrom']}\t{full_report['start']}\t{full_report['end']}\t{full_report['family']}\t.\t{full_report['strand']}\n"
    }

    with open(json_output, "w") as output:
        json.dump(output_information, output)
    
def write_output(contig_fa_outfile, te_fa_outfile, bed_outfile, json_outfile, expanded_json_outfile, vcf_outfile, reference, reference_index, *json_files):
    contigs = {}
    for file in json_files:
        with open(file, "r") as data:
            te_info = json.load(data)
            contig_name = te_info["contig_name"]
            te_name = te_info["json"]["ID"]
            if not contig_name in contigs:
                contigs[contig_name] = {"contig_path":te_info["contig_path"],"te_list":{}}
            te_dict = {
                "expanded_json":te_info["expanded_json"],
                "json":te_info["json"],
                "te_fasta":te_info["te_fasta"],
                "bed_out":te_info["bed_out"]
            }
            contigs[contig_name]["te_list"][te_name] = te_dict
    
    json_output = []
    json_expanded_output = []
    if check_exist(contig_fa_outfile):
        os.remove(contig_fa_outfile)
    with open(te_fa_outfile, "w") as te_fa, open(bed_outfile, "w") as bed_out:
        for contig in contigs:
            subprocess.run(f"cat '{contigs[contig]['contig_path']}' >> {contig_fa_outfile}", shell=True)
            te_list = contigs[contig]["te_list"]
            for te in te_list:
                outputs = te_list[te]
                json_output.append(outputs["json"])
                json_expanded_output.append(outputs["expanded_json"])
                te_fa.write(outputs["te_fasta"])
                bed_out.write(outputs["bed_out"])
    with open(json_outfile, "w") as output:
        json.dump(json_output, output, indent=4, sort_keys=False)
    with open(expanded_json_outfile, "w") as output:
        json.dump(json_expanded_output, output, indent=4, sort_keys=False)

    write_vcf(json_output, reference, reference_index, vcf_outfile)


if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
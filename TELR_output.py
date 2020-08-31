import os
import pandas as pd
import logging
import json
from Bio import SeqIO
from datetime import date


def generate_output(meta, te_fa, vcf_parsed, out, sample_name):
    logging.info("Write output...")
    # convert meta to dict
    ins_dict = dict()
    for item in meta:
        ins_dict[item['ins_name']] = item

    # generate TE sequence fasta
    final_te_seqs = out+"/"+sample_name+".final.fa"
    if os.path.isfile(final_te_seqs):
        os.remove(final_te_seqs)

    with open(te_fa, "r") as input, open(final_te_seqs, "a") as output:
        for record in SeqIO.parse(input, "fasta"):
            ins_name = record.id
            if ins_name in ins_dict:
                chr = ins_dict[ins_name]['chr']
                start = ins_dict[ins_name]['start']
                end = ins_dict[ins_name]['end']
                family = ins_dict[ins_name]['family']
                te_strand = ins_dict[ins_name]['te_strand']
                record.id = chr+"_"+str(start)+"_"+str(end)+"#"+family

                if te_strand == "+" or te_strand == ".":
                    te_seq = str(record.seq)
                else:
                    te_seq = str(record.seq.reverse_complement())
                output.write(">"+record.id+"\n"+te_seq+"\n")
                ins_dict[ins_name]['sequence'] = te_seq

    # add additional info to meta
    contig_meta = dict()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            gt = entry[10]
            ref_count = entry[11]
            alt_count = entry[12]
            contig_meta[contig_name] = {
                "gt": ':'.join([gt, ref_count, alt_count]),
                "alt_count": alt_count
            }

    # write meta data in json format
    for item in meta:
        item['gt'] = contig_meta[item['ID']]['gt']
        item['alt_count'] = contig_meta[item['ID']]['alt_count']
        del item['ins_name']
        del item['te_strand']
    report_json = out + "/" + sample_name + ".final.json"
    with open(report_json, 'w') as output:
        json.dump(meta, output, indent=4, sort_keys=False)

    # write in VCF format
    vcf_out = os.path.join(out, sample_name+'.final.vcf')
    write_vcf(meta, vcf_out)


def write_vcf(input_dict, out_vcf):
    df = pd.DataFrame(input_dict)

    df['ID'] = df.index
    df['REF'] = "N"
    df['QUAL'] = 60
    df['FILTER'] = "PASS"
    df['FORMAT'] = "GT:DR:DV"
    df['INFO'] = df.apply(lambda x: "SVTYPE=INS" + ";END=" + str(x.end) + ";FAMILY=" + str(x.family) +
                          ";STRANDS=" + str(x.strand) + ";SUPPORT_TYPE=" + str(x.support_type) + ";RE=" + str(x.alt_count) + ";AF=" + str(x.frequency), axis=1)

    df = df[['chr', 'start', 'ID', 'REF',
             'sequence', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'gt']]
    with open(out_vcf, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.1" + '\n')
        vcf.write("##fileDate={}".format(date.today()) + '\n')
        vcf.write("##source=TELR" + '\n')
        vcf.write(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">" + '\n')
        vcf.write(
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" + '\n')
        vcf.write(
            "##INFO=<ID=STRANDS,Number=A,Type=String,Description=\"Strand orientation\">" + '\n')
        vcf.write(
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" + '\n')
        vcf.write(
            "##INFO=<ID=FAMILY,Number=1,Type=String,Description=\"TE family\">" + '\n')
        vcf.write(
            "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">" + '\n')
        vcf.write("##INFO=<ID=SUPPORT_TYPE,Number=1,Type=String,Description=\"Type of TE flank alignment to reference genome (1: single flank; 2: two flanks with gap; 3: two flanks with overlap)\">" + '\n')
        vcf.write(
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" + '\n')
        vcf.write(
            "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference reads\">" + '\n')
        vcf.write(
            "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant reads\">" + '\n')
        vcf.write("#" + '\t'.join(["CHROM", "POS", "ID", "REF", "ALT",
                                   "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]) + '\n')
    df.to_csv(out_vcf, sep="\t", mode='a', index=False, header=False)

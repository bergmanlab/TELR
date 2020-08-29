import os
import subprocess
from Bio import SeqIO
import json
import re

from liftover import annotation_liftover


def annotate_contig(asm_dir, TE_library, vcf, out, sample_name, thread, presets):
    print("Generating contig annotation...")
    if presets == "ont":
        presets_minimap2 = "map-ont"
    else:
        presets_minimap2 = "map-pb"

    # merge all contigs into a single file
    merge_contigs = out+"/"+sample_name+".contigs.fa"
    contig_list = out+"/"+sample_name+".contigs.list"
    # print ("Generate merged contig file...")
    with open(merge_contigs, "w") as output_contigs, open(contig_list, "w") as output_list:
        for file in os.listdir(asm_dir):
            if ".cns.fa" in file and os.stat(asm_dir+"/"+file).st_size > 0:
                contig_name = file.replace('.cns.fa', '')
                with open(asm_dir+"/"+file, "r") as handle:
                    records = SeqIO.parse(handle, "fasta")
                    for record in records:
                        if record.id == 'ctg1':
                            record.id = contig_name
                            record.description = "len="+str(len(record.seq))
                            SeqIO.write(record, output_contigs, 'fasta')
                            output_list.write(contig_name+"\n")

    # map sequence to contigs
    seq2contig_out = out+"/"+"seq2contig.paf"
    if os.path.isfile(seq2contig_out):
        os.remove(seq2contig_out)

    # TODO: consider that some contigs might not exist
    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            vcf_seq = entry[7]
            sv_len = entry[3]
            query = out+"/"+contig_name+".seq.fa"
            create_fa(contig_name, vcf_seq, query)
            subject = out+"/"+contig_name+".fa"
            with open(subject, "w") as output:
                try:
                    subprocess.check_output(
                        ["samtools", "faidx", merge_contigs, contig_name], stderr=subprocess.DEVNULL)
                except subprocess.CalledProcessError as e:
                    print(contig_name + ":contig assembly doesn't exist")
                    continue
                else:
                    subprocess.call(
                        ["samtools", "faidx", merge_contigs, contig_name], stdout=output)
            # subject=asm_dir+"/"+contig_name+".cns.fa"
            if os.path.isfile(subject):
                with open(seq2contig_out, "a") as output:
                    mm2_output = subprocess.check_output(
                        ["minimap2", "-cx", presets_minimap2, "--secondary=no", "-v", "0", subject, query])
                    mm2_output_parsed = mm2_output.decode("utf-8")
                    if mm2_output_parsed == "":
                        print(contig_name + ':VCF sequence can not map to contig')
                    subprocess.call(["minimap2", "-cx", presets_minimap2,
                                     "--secondary=no", "-v", "0", subject, query], stdout=output)
            os.remove(query)
            os.remove(subject)
    seq2contig_bed = out+"/"+"seq2contig.bed"
    # covert to bed format
    with open(seq2contig_out, "r") as input, open(seq2contig_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            seq_id = entry[0]
            contig_id = entry[5]
            bed_line = "\t".join(
                [entry[0], entry[7], entry[8], entry[5], entry[11], entry[4]])
            output.write(bed_line+"\n")
    print("Done\n")

    # map TE library to contigs using minimap2
    # TE-contig alignment
    te2contig_out = out+"/"+sample_name+".te2contig.paf"
    print("Generating TE-contig alignment...")
    if os.path.isfile(te2contig_out):
        os.remove(te2contig_out)
    with open(contig_list, "r") as input:
        for line in input:
            contig_name = line.replace('\n', '')
            contig = out+"/"+contig_name+".fa"
            with open(contig, "w") as output:
                subprocess.call(
                    ["samtools", "faidx", merge_contigs, contig_name], stdout=output)
            # map TE library to contig using minimap2 map-pb -p 0.8 -c
            with open(te2contig_out, "a") as output:
                subprocess.call(["minimap2", "-cx", presets_minimap2, contig,
                                 TE_library, "-v", "0", "-t", str(thread)], stdout=output)
            # remove contig file
            os.remove(contig)
    # convert to bed format
    te2contig_bed = out+"/"+sample_name+".te2contig.bed"
    with open(te2contig_out, "r") as input, open(te2contig_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            bed_line = "\t".join(
                [entry[5], entry[7], entry[8], entry[0], entry[11], entry[4]])
            output.write(bed_line+"\n")
    print("Done\n")

    # Use VCF sequence alignment to filter minimap2 TE-contig alignment
    te2contig_filter_raw = out+"/"+sample_name+".te2contig_filter.tsv"
    command = "bedtools intersect -a "+te2contig_bed+" -b "+seq2contig_bed+" -wao"
    # print ("Filter TE-contig alignment...")
    with open(te2contig_filter_raw, "w") as output:
        subprocess.call(command, shell=True, stdout=output)
    # print ("Done")

    # filter and merge
    # get rid of -1 and make it into bed format
    te2contig_filter_tmp_bed = out+"/"+sample_name+".te2contig_filter.tmp.bed"
    with open(te2contig_filter_raw, "r") as input, open(te2contig_filter_tmp_bed, "w") as output:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            # the overlap between VCF sequence alignment and TE-contig alingment has to be over 10bp
            if int(entry[12]) > 10:
                out_line = "\t".join(
                    [entry[0], entry[1], entry[2], entry[3], entry[4], entry[5]])
                output.write(out_line+"\n")
    # sort
    te2contig_filter_tmp_sort_bed = out+"/" + \
        sample_name+".te2contig_filter.tmp.sort.bed"
    command = "bedtools sort -i "+te2contig_filter_tmp_bed
    with open(te2contig_filter_tmp_sort_bed, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # merge
    te2contig_filter_bed = out+"/"+sample_name+".te2contig_filter.bed"
    command = "bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim \"|\" -i " + \
        te2contig_filter_tmp_sort_bed
    with open(te2contig_filter_bed, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # remove tmp files
    os.remove(te2contig_filter_raw)
    os.remove(te2contig_filter_tmp_bed)
    os.remove(te2contig_filter_tmp_sort_bed)

    # extract sequence and RM
    te_fa = out+"/"+sample_name+".te.fa"
    with open(te_fa, "w") as output:
        subprocess.call(["bedtools", "getfasta", "-fi", merge_contigs,
                         "-bed", te2contig_filter_bed], stdout=output)
    te_rm_out = te_fa+".out.gff"
    print("Generating TE sequences repeatmasking output...")
    subprocess.call(["RepeatMasker", "-dir", out, "-gff", "-s", "-nolow", "-no_is",
                     "-xsmall", "-e", "ncbi", "-lib", TE_library, "-pa", str(thread), te_fa])
    print("Done\n")

    ## parse and merge
    te2contig_rm = out+"/"+sample_name+".te2contig_rm.bed"
    print("Repeatmask TE sequences...")
    with open(te_rm_out, "r") as input, open(te2contig_rm, "w") as output:
        for line in input:
            if "##" not in line:
                entry = line.replace('\n', '').split("\t")
                contig_name = entry[0].rsplit(':', 1)[0]
                start = entry[0].rsplit(':', 1)[1].split("-")[0]
                end = entry[0].rsplit(':', 1)[1].split("-")[1]
                # contigs = entry[0].replace(':', '-').split("-")
                family = re.sub('Target \"Motif:|\".*', '', entry[8])
                strand = entry[6]
                score = entry[5]
                out_line = "\t".join(
                    [contig_name, start, end, family, score, strand])
                output.write(out_line+"\n")
    print("Done\n")

    te2contig_rm_merge = out+"/"+sample_name+".te2contig_rm.merge.bed"
    command = "bedtools merge -c 4,6 -o distinct -delim \"|\" -i "+te2contig_rm
    with open(te2contig_rm_merge, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # build frequency dict
    te_freq = dict()
    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            freq = entry[5]
            te_freq[contig_name] = freq

    return te2contig_filter_bed, te2contig_rm_merge, te_freq, te_fa, merge_contigs


def seq2contig(seq, contig, out):
    with open(out, "a") as output:
        subprocess.call(["minimap2", "-cx", "map-pb", "--secondary=no",
                         contig, seq], stdout=output)  # only retain primary alignment


def find_te(contigs, ref, te_contigs_annotation, family_annotation, te_freq, te_fa, out, sample_name, gap, overlap, presets):
    if presets == "ont":
        presets_minimap2 = "map-ont"
    else:
        presets_minimap2 = "map-pb"

    # lift over
    report_meta = annotation_liftover(fasta1=contigs, fasta2=ref, bed=te_contigs_annotation, sample_name=sample_name, out_dir=out,
                                      preset=presets_minimap2, overlap=overlap, gap=gap, flank_len=500, family_rm=family_annotation, freq=te_freq)

    # convert meta to dict
    ins_dict = dict()
    for item in report_meta:
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

    # write meta data in json format
    for item in report_meta:
        del item['ins_name']
        del item['te_strand']
    report_json = out + "/" + sample_name + ".final.json"
    with open(report_json, 'w') as output:
        json.dump(report_meta, output, indent=4, sort_keys=False)


def create_fa(header, seq, out):
    with open(out, "w") as output:
        output.write(">"+header+"\n")
        output.write(seq)

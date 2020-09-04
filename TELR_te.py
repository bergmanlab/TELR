import os
import subprocess
from Bio import SeqIO
import json
import re
import logging
from TELR_utility import mkdir, create_loci_set
from TELR_liftover import annotation_liftover
from TELR_output import generate_output


def annotate_contig(
    contig_dir, te_library, vcf_parsed, out, sample_name, thread, presets, loci_eval
):
    logging.info("Annotate contigs...")
    if presets == "ont":
        presets = "map-ont"
    else:
        presets = "map-pb"

    all_loci = create_loci_set(vcf_parsed)
    assembly_passed_loci = set()
    merge_contigs = os.path.join(out, sample_name + ".contigs.fa")
    with open(merge_contigs, "w") as output:
        for locus in all_loci:
            assembly = os.path.join(contig_dir, locus + ".cns.fa")
            if os.path.isfile(assembly) and os.stat(assembly).st_size > 0:
                assembly_passed_loci.add(locus)
                with open(assembly, "r") as handle:
                    records = SeqIO.parse(handle, "fasta")
                    for record in records:
                        if record.id == "ctg1":
                            record.id = locus
                            record.description = "len=" + str(len(record.seq))
                            SeqIO.write(record, output, "fasta")

    # report assembly failed loci
    with open(loci_eval, "a") as output:
        for locus in all_loci:
            if locus not in assembly_passed_loci:
                output.write("\t".join([locus, "assembly failed"]) + "\n")

    # map sequence to contigs
    seq2contig_out = os.path.join(out, "seq2contig.paf")
    if os.path.isfile(seq2contig_out):
        os.remove(seq2contig_out)

    # TODO: consider that some contigs might not exist
    seq2contig_passed_loci = set()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            if contig_name in assembly_passed_loci:
                vcf_seq = entry[7]
                query = os.path.join(out, contig_name + ".seq.fa")
                create_fa(contig_name, vcf_seq, query)
                subject = os.path.join(out, contig_name + ".contig.fa")
                with open(subject, "w") as output:
                    try:
                        subprocess.check_output(
                            ["samtools", "faidx", merge_contigs, contig_name],
                            stderr=subprocess.DEVNULL,
                        )
                    except subprocess.CalledProcessError:
                        print(contig_name + ":contig assembly doesn't exist")
                        continue
                    else:
                        subprocess.call(
                            ["samtools", "faidx", merge_contigs, contig_name],
                            stdout=output,
                        )
                seq2contig_output = subprocess.check_output(
                    [
                        "minimap2",
                        "-cx",
                        presets,
                        "--secondary=no",
                        "-v",
                        "0",
                        subject,
                        query,
                    ]
                )
                seq2contig_output = seq2contig_output.decode("utf-8")
                if seq2contig_output != "":
                    seq2contig_passed_loci.add(contig_name)
                    with open(seq2contig_out, "a") as output:
                        output.write(seq2contig_output)
                os.remove(query)
                os.remove(subject)
    seq2contig_bed = os.path.join(out, "seq2contig.bed")
    # covert to bed format
    with open(seq2contig_out, "r") as input, open(seq2contig_bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            bed_line = "\t".join(
                [entry[0], entry[7], entry[8], entry[5], entry[11], entry[4]]
            )
            output.write(bed_line + "\n")

    # report ins-contig failed loci
    with open(loci_eval, "a") as output:
        for locus in assembly_passed_loci:
            if locus not in seq2contig_passed_loci:
                output.write("\t".join([locus, "seq-contig failed"]) + "\n")

    # map TE library to contigs using minimap2
    # TE-contig alignment
    te2contig_out = os.path.join(out, sample_name + ".te2contig.paf")
    if os.path.isfile(te2contig_out):
        os.remove(te2contig_out)
    for locus in seq2contig_passed_loci:
        contig_fa = os.path.join(out, locus + ".fa")
        with open(contig_fa, "w") as output:
            subprocess.call(["samtools", "faidx", merge_contigs, locus], stdout=output)
        # map TE library to contig using minimap2
        with open(te2contig_out, "a") as output:
            subprocess.call(
                [
                    "minimap2",
                    "-cx",
                    presets,
                    contig_fa,
                    te_library,
                    "-v",
                    "0",
                    "-t",
                    str(thread),
                ],
                stdout=output,
            )
        os.remove(contig_fa)
    # convert to bed format
    te2contig_bed = os.path.join(out, sample_name + ".te2contig.bed")
    with open(te2contig_out, "r") as input, open(te2contig_bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            bed_line = "\t".join(
                [entry[5], entry[7], entry[8], entry[0], entry[11], entry[4]]
            )
            output.write(bed_line + "\n")

    # Use VCF sequence alignment to filter minimap2 TE-contig alignment
    te2contig_filter_raw = os.path.join(out, sample_name + ".te2contig_filter.tsv")
    with open(te2contig_filter_raw, "w") as output:
        subprocess.call(
            [
                "bedtools",
                "intersect",
                "-a",
                te2contig_bed,
                "-b",
                seq2contig_bed,
                "-wao",
            ],
            stdout=output,
        )

    # filter and merge
    # get rid of -1 and make it into bed format
    te2contig_filter_tmp_bed = os.path.join(
        out, sample_name + ".te2contig_filter.tmp.bed"
    )
    with open(te2contig_filter_raw, "r") as input, open(
        te2contig_filter_tmp_bed, "w"
    ) as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            # the overlap between VCF sequence alignment and TE-contig alignment has to be over 10bp
            if int(entry[12]) > 10:
                out_line = "\t".join(
                    [entry[0], entry[1], entry[2], entry[3], entry[4], entry[5]]
                )
                output.write(out_line + "\n")
    # sort
    te2contig_filter_tmp_sort_bed = (
        out + "/" + sample_name + ".te2contig_filter.tmp.sort.bed"
    )
    command = "bedtools sort -i " + te2contig_filter_tmp_bed
    with open(te2contig_filter_tmp_sort_bed, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # merge
    te2contig_filter_bed = out + "/" + sample_name + ".te2contig_filter.bed"
    command = (
        'bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim "|" -i '
        + te2contig_filter_tmp_sort_bed
    )
    with open(te2contig_filter_bed, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # remove tmp files
    os.remove(te2contig_out)
    os.remove(seq2contig_bed)
    os.remove(te2contig_filter_raw)
    os.remove(te2contig_filter_tmp_bed)
    os.remove(te2contig_filter_tmp_sort_bed)

    # extract sequence and RM
    te_fa = out + "/" + sample_name + ".te.fa"
    with open(te_fa, "w") as output:
        subprocess.call(
            [
                "bedtools",
                "getfasta",
                "-fi",
                merge_contigs,
                "-bed",
                te2contig_filter_bed,
            ],
            stdout=output,
        )
    repeatmasker_dir = os.path.join(out, "contig_te_repeatmask")
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
                te_fa,
            ]
        )
        contig_te_repeatmasked = os.path.join(
            repeatmasker_dir, os.path.basename(te_fa) + ".out.gff"
        )
        open(contig_te_repeatmasked, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking contig TE sequences failed, exiting...")
        sys.exit(1)

    ## parse and merge
    te2contig_rm = out + "/" + sample_name + ".te2contig_rm.bed"
    with open(contig_te_repeatmasked, "r") as input, open(te2contig_rm, "w") as output:
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
                out_line = "\t".join([contig_name, start, end, family, score, strand])
                output.write(out_line + "\n")
    print("Done\n")

    te2contig_rm_merge = out + "/" + sample_name + ".te2contig_rm.merge.bed"
    command = 'bedtools merge -c 4,6 -o distinct -delim "|" -i ' + te2contig_rm
    with open(te2contig_rm_merge, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # build frequency dict
    te_freq = dict()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            freq = entry[5]
            te_freq[contig_name] = freq

    return te2contig_filter_bed, te2contig_rm_merge, te_freq, te_fa, merge_contigs


def seq2contig(seq, contig, out):
    with open(out, "a") as output:
        subprocess.call(
            ["minimap2", "-cx", "map-pb", "--secondary=no", contig, seq], stdout=output
        )  # only retain primary alignment


# def find_te(contigs, ref, te_contigs_annotation, family_annotation, te_freq, te_fa, out, sample_name, gap, overlap, presets):
def find_te(
    contig_dir,
    vcf_parsed,
    ref,
    te_library,
    out,
    sample_name,
    thread,
    gap,
    overlap,
    presets,
    loci_eval,
):
    """
    Identify non-reference TE insertions in the reference genome using assembled contigs
    """
    if presets == "ont":
        presets = "map-ont"
    else:
        presets = "map-pb"

    # contig annotation
    (
        te_contigs_annotation,
        family_annotation,
        te_freq,
        te_fa,
        merge_contigs,
    ) = annotate_contig(
        contig_dir, te_library, vcf_parsed, out, sample_name, thread, presets, loci_eval
    )

    # lift over
    logging.info("Map contigs to reference...")
    report_meta = annotation_liftover(
        fasta1=merge_contigs,
        fasta2=ref,
        bed=te_contigs_annotation,
        sample_name=sample_name,
        out_dir=out,
        preset=presets,
        overlap=overlap,
        gap=gap,
        flank_len=500,
        family_rm=family_annotation,
        freq=te_freq,
    )

    generate_output(report_meta, te_fa, vcf_parsed, out, sample_name)


def create_fa(header, seq, out):
    with open(out, "w") as output:
        output.write(">" + header + "\n")
        output.write(seq)

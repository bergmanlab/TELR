import os
import subprocess
from Bio import SeqIO
import json
import re
import logging
import pysam
import statistics
from TELR_utility import mkdir, create_loci_set
from TELR_liftover import annotation_liftover
from TELR_output import generate_output
from TELR_assembly import prep_assembly


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
    with open(merge_contigs, "w") as MERGE:
        for locus in all_loci:
            assembly = os.path.join(contig_dir, locus + ".cns.fa")
            new_assembly = os.path.join(contig_dir, locus + ".cns.ctg1.fa")
            if os.path.isfile(assembly) and os.stat(assembly).st_size > 0:
                assembly_passed_loci.add(locus)
                with open(assembly, "r") as handle:
                    records = SeqIO.parse(handle, "fasta")
                    for record in records:
                        if record.id == "ctg1":
                            record.id = locus
                            record.description = "len=" + str(len(record.seq))
                            SeqIO.write(record, MERGE, "fasta")
                            with open(new_assembly, "w") as CTG1:
                                SeqIO.write(record, CTG1, "fasta")

    # report assembly failed loci
    with open(loci_eval, "a") as output:
        for locus in all_loci:
            if locus not in assembly_passed_loci:
                output.write("\t".join([locus, "Contig assembly failed"]) + "\n")

    # map sequence to contigs
    seq2contig_out = os.path.join(out, "seq2contig.paf")
    if os.path.isfile(seq2contig_out):
        os.remove(seq2contig_out)

    # TODO: consider that some contigs might not exist
    seq2contig_passed_loci = set()
    seq2contig_dir = os.path.join(out, "seq2contig")
    seq2contig = os.path.join(out, "seq2contig.paf")
    mkdir(seq2contig_dir)
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            if contig_name in assembly_passed_loci:
                vcf_seq = entry[7]
                query = os.path.join(seq2contig_dir, contig_name + ".seq.fa")
                create_fa(contig_name, vcf_seq, query)
                subject = os.path.join(seq2contig_dir, contig_name + ".contig.fa")
                with open(subject, "w") as output:
                    try:
                        subprocess.call(
                            ["samtools", "faidx", merge_contigs, contig_name],
                            stdout=output,
                        )
                    except subprocess.CalledProcessError:
                        print(contig_name + ":contig assembly doesn't exist")
                        continue
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
                    with open(seq2contig, "a") as output:
                        output.write(seq2contig_output)
                os.remove(query)
                os.remove(subject)
    os.rmdir(seq2contig_dir)
    # covert to bed format
    seq2contig_bed = os.path.join(out, "seq2contig.bed")
    with open(seq2contig, "r") as input, open(seq2contig_bed, "w") as output:
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
                output.write(
                    "\t".join(
                        [locus, "Sniffles VCF sequence not mapped to assembled contig"]
                    )
                    + "\n"
                )

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

    # find out what's filtered out
    seq_mm2_overlap_loci = set()
    with open(te2contig_filter_tmp_sort_bed, "r") as input:
        for line in input:
            seq_mm2_overlap_loci.add(line.split("\t")[0])
    # seq_mm2_overlap_loci = create_loci_set(te2contig_filter_tmp_sort_bed)
    with open(loci_eval, "a") as output:
        for locus in seq2contig_passed_loci:
            if locus not in seq_mm2_overlap_loci:
                output.write(
                    "\t".join([locus, "VCF sequence doesn't overlap contig annotation"])
                    + "\n"
                )

    # merge
    contig_te_annotation = out + "/" + sample_name + ".te2contig_filter.bed"
    command = (
        'bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim "|" -i '
        + te2contig_filter_tmp_sort_bed
    )
    with open(contig_te_annotation, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    # seq_mm2_overlap_merge_loci = create_loci_set(contig_te_annotation)

    # remove tmp files
    os.remove(seq2contig)
    os.remove(te2contig_bed)
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
                contig_te_annotation,
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

    contig_rm_annotation = out + "/" + sample_name + ".te2contig_rm.merge.bed"
    command = 'bedtools merge -c 4,6 -o distinct -delim "|" -i ' + te2contig_rm
    with open(contig_rm_annotation, "w") as output:
        subprocess.call(command, shell=True, stdout=output)
    os.remove(te2contig_rm)

    # seq_mm2_overlap_merge_rm_loci = create_loci_set(te2contig_rm_merge)
    # with open(loci_eval, "a") as output:
    #     for locus in seq_mm2_overlap_merge_loci:
    #         if locus not in seq_mm2_overlap_merge_rm_loci:
    #             print(locus, "contig seq RM failed")

    # build frequency dict
    te_freq = dict()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            freq = entry[5]
            te_freq[contig_name] = freq

    return contig_te_annotation, contig_rm_annotation, te_freq, te_fa, merge_contigs


def seq2contig(seq, contig, out):
    with open(out, "a") as output:
        subprocess.call(
            ["minimap2", "-cx", "map-pb", "--secondary=no", contig, seq], stdout=output
        )  # only retain primary alignment


def minimap2bam(query, subject, presets, thread=1):
    sam = query.replace(".reads.fa", ".realign.sam")
    with open(sam, "w") as output:
        subprocess.call(
            ["minimap2", "-a", "-x", presets, "-v", "0", subject, query], stdout=output
        )
    bam = query.replace(".reads.fa", ".realign.bam")
    with open(bam, "w") as output:
        subprocess.call(["samtools", "view", "-bS", sam], stdout=output)
    sorted_bam = query.replace(".reads.fa", ".realign.sort.bam")
    subprocess.call(["samtools", "sort", "-@", str(thread), "-o", sorted_bam, bam])
    subprocess.call(["samtools", "index", "-@", str(thread), sorted_bam])

    os.remove(sam)
    os.remove(bam)

    return sorted_bam


def get_af(
    out,
    sample_name,
    bam,
    raw_reads,
    contig_te_annotation,
    contig_dir,
    vcf_parsed,
    presets,
    thread,
):
    logging.info("Estimating allele frequency")
    if presets == "ont":
        presets = "map-ont"
    else:
        presets = "map-pb"

    # prepare reads
    telr_reads_dir = os.path.join(out, "telr_reads")
    prep_assembly(
        vcf_parsed, out, sample_name, bam, raw_reads, telr_reads_dir, read_type="all"
    )

    # read contig annotation to dict
    contig_te_dict = dict()
    with open(contig_te_annotation, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            start = entry[1]
            end = entry[2]
            contig_te_dict[contig_name] = int(start), int(end)

    # re-align reads to assembly
    k = 0
    vcf_parsed_freq = vcf_parsed + ".freq"
    with open(vcf_parsed, "r") as input, open(vcf_parsed_freq, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            # rename all reads around breakpoints
            telr_reads = telr_reads_dir + "/contig" + str(k)
            telr_reads_rename = telr_reads_dir + "/" + contig_name + ".reads.fa"
            os.rename(telr_reads, telr_reads_rename)
            k = k + 1

            contig = os.path.join(contig_dir, contig_name + ".cns.ctg1.fa")

            # get contig length
            contig_length = 0
            if os.path.isfile(contig):
                with open(contig, "r") as handle:
                    records = SeqIO.parse(handle, "fasta")
                    for record in records:
                        contig_length = len(record.seq)
            else:
                print(contig_name + " no assembly")
                continue

            # map reads to assembly
            raw_reads = os.path.join(telr_reads_dir, contig_name + ".reads.fa")
            # bam = raw_reads.replace(".reads.fa", ".realign.sort.bam")
            bam = minimap2bam(raw_reads, contig, presets, thread)

            # get average
            if contig_name in contig_te_dict:
                start, end = contig_te_dict[contig_name]
                # get TE coverage
                te_cov = get_median_cov(bam, contig_name, start, end)
                # get flanking coverage
                flank_len = 1000
                offset = 200
                left_flank_present = True
                right_flank_present = True
                left_flank_cov = 0
                right_flank_cov = 0
                if start - flank_len - offset >= 0:
                    left_flank_cov = get_median_cov(
                        bam, contig_name, start - flank_len - offset, start - offset
                    )
                else:
                    left_flank_present = False

                if end + flank_len + offset <= contig_length:
                    right_flank_cov = get_median_cov(
                        bam, contig_name, end + offset, end + flank_len + offset
                    )
                else:
                    right_flank_present = False

                if left_flank_present and right_flank_present:
                    flank_cov = (left_flank_cov + right_flank_cov) / 2
                elif left_flank_present:
                    flank_cov = left_flank_cov
                elif right_flank_present:
                    flank_cov = right_flank_cov
                else:
                    print(contig_name + " no flanks")
                    continue
                out_line = "\t".join(
                    entry
                    + [
                        str(te_cov),
                        str(left_flank_cov),
                        str(right_flank_cov),
                        str(flank_cov),
                    ]
                )
                output.write(out_line + "\n")
            else:
                print(contig_name + " not in contig_te")

    # get frequency
    te_freq = dict()
    with open(vcf_parsed_freq, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            te_cov = float(entry[13])
            flank_cov = float(entry[16])
            freq = round(te_cov / flank_cov, 2)
            if freq > 1:
                freq = 1
            te_freq[contig_name] = freq
    return te_freq

def get_median_cov(bam, chr, start, end):
    commands = (
        "samtools depth -aa -r " + chr + ":" + str(start) + "-" + str(end) + " " + bam
    )
    depth = bam + ".depth"
    with open(depth, "w") as output:
        subprocess.call(commands, shell=True, stdout=output)
    covs = []
    with open(depth, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            covs.append(int(entry[2]))
    median_cov = statistics.median(covs)
    os.remove(depth)
    return median_cov


# def find_te(contigs, ref, te_contigs_annotation, family_annotation, te_freq, te_fa, out, sample_name, gap, overlap, presets):
def find_te(
    contig_te_annotation,
    contig_family_annotation,
    te_freq,
    merge_contigs,
    ref,
    out,
    sample_name,
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

    # lift over
    logging.info("Map contigs to reference...")
    report_meta = annotation_liftover(
        fasta1=merge_contigs,
        fasta2=ref,
        bed=contig_te_annotation,
        sample_name=sample_name,
        out_dir=out,
        preset=presets,
        overlap=overlap,
        gap=gap,
        flank_len=500,
        family_rm=contig_family_annotation,
        freq=te_freq,
        loci_eval=loci_eval,
    )

    return report_meta


def create_fa(header, seq, out):
    with open(out, "w") as output:
        output.write(">" + header + "\n")
        output.write(seq)

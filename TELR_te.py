import os
import sys
import subprocess
from Bio import SeqIO
import re
import logging
import time
import statistics
from multiprocessing import Pool
from TELR_utility import mkdir, create_loci_set, format_time

# from TELR_liftover import annotation_liftover
from TELR_liftover import liftover
from TELR_assembly import prep_assembly


def annotate_contig(
    contig_dir,
    te_library,
    vcf_parsed,
    out,
    sample_name,
    thread,
    presets,
    minimap2_family,
    loci_eval,
):
    logging.info("Annotate contigs...")
    if presets == "pacbio":
        presets = "map-pb"
    else:
        presets = "map-ont"

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
                        if record.id == "ctg1" or record.id == "contig_1":
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
                    with open(seq2contig_out, "a") as output:
                        output.write(seq2contig_output)
                os.remove(query)
                os.remove(subject)
    os.rmdir(seq2contig_dir)
    # covert to bed format
    seq2contig_bed = os.path.join(out, "seq2contig.bed")
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
    # sort # TODO: package this part, hide variables
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
    contig_te_annotation_tmp = out + "/" + sample_name + ".te2contig_filter.bed.tmp"
    command = (
        'bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim "|" -i '
        + te2contig_filter_tmp_sort_bed
    )
    with open(contig_te_annotation_tmp, "w") as output:
        subprocess.call(command, shell=True, stdout=output)

    contig_te_annotation = out + "/" + sample_name + ".te2contig_filter.bed"
    with open(contig_te_annotation_tmp, "r") as input, open(
        contig_te_annotation, "w"
    ) as output:
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
                    contig_te_strand,
                ]
            )
            output.write(out_line + "\n")

    # seq_mm2_overlap_merge_loci = create_loci_set(contig_te_annotation)

    # remove tmp files
    os.remove(te2contig_bed)
    os.remove(te2contig_out)
    os.remove(seq2contig_bed)
    os.remove(te2contig_filter_raw)
    os.remove(te2contig_filter_tmp_bed)
    os.remove(te2contig_filter_tmp_sort_bed)

    # extract sequence and RM
    if "+" in sample_name:
        sample_name_replace = sample_name.replace("+", "plus")
    else:
        sample_name_replace = sample_name
    te_fa = out + "/" + sample_name_replace + ".te.fa"
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

    if not minimap2_family:
        print("Use repeatmasker to annotate contig TE families instead of minimap2")
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
        with open(contig_te_repeatmasked, "r") as input, open(
            te2contig_rm, "w"
        ) as output:
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
        print("Done\n")

        contig_rm_annotation = out + "/" + sample_name + ".te2contig_rm.merge.bed"
        command = 'bedtools merge -c 4,6 -o distinct -delim "|" -i ' + te2contig_rm
        with open(contig_rm_annotation, "w") as output:
            subprocess.call(command, shell=True, stdout=output)
        # os.remove(te2contig_rm)

        # replace contig_te_annotation family with ones from RM
        contig_te_annotation_new = contig_te_annotation.replace(
            "bed", "family_reannotated.bed"
        )
        contig_rm_family_dict = dict()
        with open(contig_rm_annotation, "r") as input:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                contig_name = entry[0]
                family = entry[3]
                contig_rm_family_dict[contig_name] = family

        with open(contig_te_annotation_new, "w") as output, open(
            contig_te_annotation, "r"
        ) as input:
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

        contig_te_annotation = contig_te_annotation_new

    # build frequency dict
    te_freq = dict()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            freq = entry[5]
            te_freq[contig_name] = freq

    return contig_te_annotation, te_freq, te_fa, merge_contigs


def seq2contig(seq, contig, out):
    with open(out, "a") as output:
        subprocess.call(
            ["minimap2", "-cx", "map-pb", "--secondary=no", contig, seq], stdout=output
        )  # only retain primary alignment


def repeatmask(ref, library, outdir, thread):
    mkdir(outdir)
    try:
        subprocess.call(
            [
                "RepeatMasker",
                "-dir",
                outdir,
                "-gff",
                "-s",
                "-nolow",
                "-no_is",
                "-e",
                "ncbi",
                "-lib",
                library,
                "-pa",
                str(thread),
                ref,
            ]
        )
        ref_rm = os.path.join(outdir, os.path.basename(ref) + ".masked")
        gff = os.path.join(outdir, os.path.basename(ref) + ".out.gff")
        gff3 = os.path.join(outdir, os.path.basename(ref) + ".out.gff3")
        if not os.path.isfile(ref_rm):
            ref_rm_out = os.path.join(outdir, os.path.basename(ref) + ".out")
            with open(ref_rm_out, "r") as input:
                for line in input:
                    if "There were no repetitive sequences detected" in line:
                        print("No repetitive sequences detected")
                        ref_rm = ref
                        gff = None
                        gff3 = None
                    else:
                        raise Exception("Repeatmasking failed, exiting...")
        else:
            parse_rm_out(gff, gff3)
            open(ref_rm, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking failed, exiting...")
        sys.exit(1)
    return ref_rm, gff3


def gff3tobed(gff, bed):
    # check GFF3 format
    with open(gff, "r") as input:
        for line in input:
            if "#" not in line:
                if "Target=" not in line:
                    print(
                        "Incorrect GFF3 format, please check README for expected format, exiting..."
                    )
                    logging.exception(
                        "Incorrect GFF3 format, please check README for expected format, exiting..."
                    )
                    sys.exit(1)
                break
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


def parse_rm_out(rm_gff, gff3):
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
                        "Target=" + family,
                    ]
                )
                output.write(out_line + "\n")


def realignment(args):
    query = args[0]
    subject = args[1]
    presets = args[2]
    thread = 1

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


def get_flank_cov(
    bam, contig_name, contig_length, start, end, flank_len=100, offset=200
):
    left_flank_present = True
    right_flank_present = True
    left_flank_cov = 0
    right_flank_cov = 0
    flank_cov = 0
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
        print(contig_name + " no flank present")

    return left_flank_cov, right_flank_cov, flank_cov


def get_contig_length(contig):
    if os.path.isfile(contig):
        with open(contig, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                contig_length = len(record.seq)
                return contig_length
    else:
        print("no contig " + contig)


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
    logging.info("Estimating allele frequency...")
    start_time = time.time()
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
    realign_pa_list = []
    with open(vcf_parsed, "r") as input:
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
            if not os.path.isfile(contig):
                print(contig_name + " no assembly")
                continue

            # prepare reads for assembly
            raw_reads = os.path.join(telr_reads_dir, contig_name + ".reads.fa")

            align_pa = [raw_reads, contig, presets]
            realign_pa_list.append(align_pa)

    # run realignment in parallel
    logging.info("Perform local realignment...")
    start_time = time.time()
    try:
        pool = Pool(processes=thread)
        pool.map(realignment, realign_pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("Local realignment failed, exiting...")
        sys.exit(1)

    proc_time = time.time() - start_time
    logging.info("Local realignment finished in " + format_time(proc_time))

    # analyze realignment and estimate coverage
    vcf_parsed_freq = vcf_parsed + ".freq"
    flank_len = 100
    offset = 200
    with open(vcf_parsed, "r") as input, open(vcf_parsed_freq, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            bam = os.path.join(telr_reads_dir, contig_name + ".realign.sort.bam")
            if os.path.isfile(bam):
                if contig_name in contig_te_dict:
                    start, end = contig_te_dict[contig_name]
                    # get contig size
                    contig = os.path.join(contig_dir, contig_name + ".cns.ctg1.fa")
                    contig_length = get_contig_length(contig)
                    # get TE coverage
                    te_cov = get_median_cov(bam, contig_name, start, end)
                    left_flank_cov, right_flank_cov, flank_cov = get_flank_cov(
                        bam, contig_name, contig_length, start, end, flank_len, offset
                    )
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

    # get frequency
    te_freq = dict()
    with open(vcf_parsed_freq, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            te_cov = float(entry[14])
            flank_cov = float(entry[17])
            if flank_cov == 0:
                freq = 1
                print(contig_name + " zero flank coverage for")
            else:
                freq = round(te_cov / flank_cov, 2)
            if freq > 1:
                freq = 1
            te_freq[contig_name] = freq
    proc_time = time.time() - start_time
    logging.info("Allele frequency estimation finished in " + format_time(proc_time))
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

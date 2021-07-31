import os
import sys
import subprocess
from Bio import SeqIO
import re
import logging
import time
import statistics
from multiprocessing import Pool
from telr.TELR_utility import (
    check_exist,
    mkdir,
    format_time,
    get_cmd_output,
    get_rev_comp_sequence,
)
from telr.TELR_liftover import liftover
from telr.TELR_assembly import prep_assembly_inputs


def annotate_contig(
    contigs,
    assembly_passed_loci,
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
        minimap2_presets = "map-pb"
    else:
        minimap2_presets = "map-ont"

    # map sequence to contigs
    vcf_seq2contig_out = os.path.join(out, "seq2contig.paf")
    # if os.path.isfile(vcf_seq2contig_out):
    #     os.remove(vcf_seq2contig_out)

    # TODO: consider that some contigs might not exist
    seq2contig_passed_loci = set()
    vcf_seq2contig_dir = os.path.join(out, "vcf_seq2contig")
    mkdir(vcf_seq2contig_dir)
    with open(vcf_parsed, "r") as input, open(vcf_seq2contig_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            if contig_name in assembly_passed_loci:
                vcf_seq = entry[7]
                query = os.path.join(vcf_seq2contig_dir, contig_name + ".seq.fa")
                create_fa(contig_name, vcf_seq, query)
                subject = os.path.join(
                    vcf_seq2contig_dir, contig_name + ".contig.fa"
                )  ## TODO: this can be replaced
                with open(subject, "w") as subject_output_handle:
                    try:
                        subprocess.call(
                            ["samtools", "faidx", contigs, contig_name],
                            stdout=subject_output_handle,
                        )
                    except subprocess.CalledProcessError:
                        print(contig_name + ":contig assembly doesn't exist")
                        continue
                cmd = [
                    "minimap2",
                    "-cx",
                    minimap2_presets,
                    "--secondary=no",
                    "-v",
                    "0",
                    subject,
                    query,
                ]
                vcf_seq2contig_output = get_cmd_output(cmd)
                if vcf_seq2contig_output != "":
                    output.write(vcf_seq2contig_output)
                    seq2contig_passed_loci.add(contig_name)
                    # with open(vcf_seq2contig_out, "a") as output:
                os.remove(query)
                os.remove(subject)
    os.rmdir(vcf_seq2contig_dir)

    # covert to bed format
    seq2contig_bed = os.path.join(out, "seq2contig.bed")
    with open(vcf_seq2contig_out, "r") as input, open(seq2contig_bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            bed_line = "\t".join(
                [entry[0], entry[7], entry[8], entry[5], entry[11], entry[4]]
            )
            output.write(bed_line + "\n")

    # # report ins-contig failed loci
    # with open(loci_eval, "a") as output:
    #     for locus in assembly_passed_loci:
    #         if locus not in seq2contig_passed_loci:
    #             output.write(
    #                 "\t".join(
    #                     [locus, "Sniffles VCF sequence not mapped to assembled contig"]
    #                 )
    #                 + "\n"
    #             )

    # map TE library to contigs using minimap2
    # TE-contig alignment
    te2contig_out = os.path.join(out, sample_name + ".te2contig.paf")
    if os.path.isfile(te2contig_out):
        os.remove(te2contig_out)
    for locus in seq2contig_passed_loci:
        contig_fa = os.path.join(out, locus + ".fa")
        with open(contig_fa, "w") as output:
            subprocess.call(["samtools", "faidx", contigs, locus], stdout=output)
        # map TE library to contig using minimap2
        with open(te2contig_out, "a") as output:
            subprocess.call(
                [
                    "minimap2",
                    "-cx",
                    minimap2_presets,
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
                contigs,
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

    return contig_te_annotation, te_fa


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
    prefix = args[3]
    thread = 1

    sam = prefix + ".sam"
    with open(sam, "w") as output:
        subprocess.call(
            ["minimap2", "-a", "-x", presets, "-v", "0", subject, query], stdout=output
        )
    bam = prefix + ".realign.bam"
    with open(bam, "w") as output:
        subprocess.call(["samtools", "view", "-bS", sam], stdout=output)
    sorted_bam = prefix + ".realign.sort.bam"
    subprocess.call(["samtools", "sort", "-@", str(thread), "-o", sorted_bam, bam])
    subprocess.call(["samtools", "index", "-@", str(thread), sorted_bam])

    os.remove(sam)
    os.remove(bam)


def get_flank_cov(bam, contig_name, contig_length, start, end, flank_len, offset):
    """
    Get the coverage of the flanking regions of a TE
    """
    # left_flank_present = True
    # right_flank_present = True
    left_flank_cov = None
    right_flank_cov = None
    # flank_cov = 0
    if start - flank_len - offset >= 0:
        left_flank_cov = get_median_cov(
            bam, contig_name, start - flank_len - offset, start - offset
        )
    # else:
    #     left_flank_present = False

    if end + flank_len + offset <= contig_length:
        right_flank_cov = get_median_cov(
            bam, contig_name, end + offset, end + flank_len + offset
        )
    # else:
    #     right_flank_present = False

    # if left_flank_present and right_flank_present:
    #     flank_cov = (left_flank_cov + right_flank_cov) / 2
    # elif left_flank_present:
    #     flank_cov = left_flank_cov
    # elif right_flank_present:
    #     flank_cov = right_flank_cov
    # else:
    #     print(contig_name + " no flank present")

    return left_flank_cov, right_flank_cov


def get_contig_length(contig):
    if os.path.isfile(contig):
        with open(contig, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                contig_length = len(record.seq)
                return contig_length
    else:
        print("no contig " + contig)


def get_te_flank_ratio(te_cov, flank_cov):
    if te_cov and flank_cov:
        if flank_cov == 0:
            return None
        else:
            ratio = te_cov / flank_cov
            if ratio > 1.5:
                return None
            else:
                return ratio
    else:
        return None


def get_af(
    out,
    sample_name,
    bam,
    raw_reads,
    contig_te_annotation,
    contig_dir,
    vcf_parsed,
    flank_intervel_size,
    flank_offset,
    te_interval_size,
    te_offset,
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
    prep_assembly_inputs(
        vcf_parsed, out, sample_name, bam, raw_reads, telr_reads_dir, read_type="all"
    )

    # re-align reads to assembly, both forward and reverse
    k = 0
    realign_pa1_list = []
    realign_pa2_list = []
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
            if not os.path.isfile(contig):
                print(contig_name + " no assembly")
                continue
            contig_rev_comp = os.path.join(
                contig_dir, contig_name + ".cns.ctg1.revcomp.fa"
            )
            get_rev_comp_sequence(contig, contig_rev_comp)

            # prepare reads for assembly
            raw_reads = os.path.join(telr_reads_dir, contig_name + ".reads.fa")

            prefix1 = os.path.join(telr_reads_dir, contig_name)
            align_pa1 = [raw_reads, contig, presets, prefix1]
            realign_pa1_list.append(align_pa1)

            prefix2 = os.path.join(telr_reads_dir, contig_name + ".revcomp")
            align_pa2 = [raw_reads, contig_rev_comp, presets, prefix2]
            realign_pa2_list.append(align_pa2)

    # run realignment in parallel
    logging.info("Perform local realignment...")
    start_time = time.time()
    try:
        pool = Pool(processes=thread)
        pool.map(realignment, realign_pa1_list)
        pool.map(realignment, realign_pa2_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("Local realignment failed, exiting...")
        sys.exit(1)
    proc_time = time.time() - start_time
    logging.info("Local realignment finished in " + format_time(proc_time))

    # read contig annotation to dict
    contig_te_coord_dict = dict()
    with open(contig_te_annotation, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = entry[0]
            # get coordinates of TEs on reverse complement
            contig = os.path.join(contig_dir, contig_name + ".cns.ctg1.fa")
            if check_exist(contig):
                start = int(entry[1])
                end = int(entry[2])
                contig_length = get_contig_length(contig)
                start_revcomp = contig_length - end
                end_revcomp = contig_length - start
                contig_te_coord_dict[contig_name] = {
                    "fw": (start, end),
                    "rc": (start_revcomp, end_revcomp),
                }
            else:
                continue

    # analyze realignment and estimate coverage
    vcf_parsed_freq = vcf_parsed + ".freq"
    with open(vcf_parsed, "r") as input, open(vcf_parsed_freq, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            bam = os.path.join(telr_reads_dir, contig_name + ".realign.sort.bam")
            if os.path.isfile(bam):
                if contig_name in contig_te_coord_dict:
                    start, end = contig_te_coord_dict[contig_name]["fw"]
                    # get contig size
                    contig = os.path.join(contig_dir, contig_name + ".cns.ctg1.fa")
                    contig_length = get_contig_length(contig)
                    # get TE locus coverage
                    te_5p_cov, te_3p_cov = get_te_cov(
                        bam, contig_name, start, end, te_interval_size, te_offset
                    )
                    # get flanking coverage
                    flank_5p_cov, flank_3p_cov = get_flank_cov(
                        bam,
                        contig_name,
                        contig_length,
                        start,
                        end,
                        flank_intervel_size,
                        flank_offset,
                    )
                    out_line = "\t".join(
                        entry
                        + [
                            str(te_5p_cov),
                            str(te_3p_cov),
                            str(flank_5p_cov),
                            str(flank_3p_cov),
                        ]
                    )
                    output.write(out_line + "\n")

    # analyze realignment and estimate coverage, on rev-comp contigs
    vcf_parsed_freq_revcomp = vcf_parsed + ".revcomp.freq"
    with open(vcf_parsed, "r") as input, open(vcf_parsed_freq_revcomp, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            bam = os.path.join(
                telr_reads_dir, contig_name + ".revcomp.realign.sort.bam"
            )
            if os.path.isfile(bam):
                if contig_name in contig_te_coord_dict:
                    start, end = contig_te_coord_dict[contig_name]["rc"]
                    # get contig size
                    contig = os.path.join(
                        contig_dir, contig_name + ".cns.ctg1.revcomp.fa"
                    )
                    contig_length = get_contig_length(contig)
                    # get TE locus coverage
                    te_5p_cov, te_3p_cov = get_te_cov(
                        bam, contig_name, start, end, te_interval_size, te_offset
                    )
                    # get flanking coverage
                    flank_5p_cov, flank_3p_cov = get_flank_cov(
                        bam,
                        contig_name,
                        contig_length,
                        start,
                        end,
                        flank_intervel_size,
                        flank_offset,
                    )
                    out_line = "\t".join(
                        entry
                        + [
                            str(te_5p_cov),
                            str(te_3p_cov),
                            str(flank_5p_cov),
                            str(flank_3p_cov),
                        ]
                    )
                    output.write(out_line + "\n")

    # get frequency
    te_freq = dict()
    with open(vcf_parsed_freq, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            te_freq[contig_name] = dict()
            if entry[14] != "None":
                te_5p_cov = float(entry[14])
            else:
                te_5p_cov = None
            if entry[15] != "None":
                te_3p_cov = float(entry[15])
            else:
                te_3p_cov = None
            if entry[16] != "None":
                flank_5p_cov = float(entry[16])
            else:
                flank_5p_cov = None
            if entry[17] != "None":
                flank_3p_cov = float(entry[17])
            else:
                flank_3p_cov = None
            te_freq[contig_name]["te_5p_cov"] = te_5p_cov
            te_freq[contig_name]["te_3p_cov"] = te_3p_cov
            te_freq[contig_name]["flank_5p_cov"] = flank_5p_cov
            te_freq[contig_name]["flank_3p_cov"] = flank_3p_cov

    # get TAF on reverse-complemented contigs
    with open(vcf_parsed_freq_revcomp, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            if entry[14] != "None":
                te_5p_cov = float(entry[14])
            else:
                te_5p_cov = None
            if entry[15] != "None":
                te_3p_cov = float(entry[15])
            else:
                te_3p_cov = None
            if entry[16] != "None":
                flank_5p_cov = float(entry[16])
            else:
                flank_5p_cov = None
            if entry[17] != "None":
                flank_3p_cov = float(entry[17])
            else:
                flank_3p_cov = None
            te_freq[contig_name]["te_5p_cov_rc"] = te_5p_cov
            te_freq[contig_name]["te_3p_cov_rc"] = te_3p_cov
            te_freq[contig_name]["flank_5p_cov_rc"] = flank_5p_cov
            te_freq[contig_name]["flank_3p_cov_rc"] = flank_3p_cov
            taf_5p = get_te_flank_ratio(
                te_freq[contig_name]["te_5p_cov"],
                te_freq[contig_name]["flank_5p_cov"],
            )
            taf_3p = get_te_flank_ratio(
                te_freq[contig_name]["te_5p_cov_rc"],
                te_freq[contig_name]["flank_5p_cov_rc"],
            )
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
            if freq:
                te_freq[contig_name]["freq"] = round(freq, 3)
            else:
                te_freq[contig_name]["freq"] = None
    proc_time = time.time() - start_time
    logging.info("Allele frequency estimation finished in " + format_time(proc_time))
    return te_freq


def get_te_cov(bam, contig_name, start, end, te_interval_size, te_offset):
    """
    Get TE locus coverage
    """
    te_5p_cov = None
    te_3p_cov = None
    whole_te_locus_cov = False
    if te_interval_size:
        if start + te_offset + te_interval_size < end:
            te_5p_cov = get_median_cov(
                bam,
                contig_name,
                start + te_offset,
                start + te_offset + te_interval_size,
            )
            te_3p_cov = get_median_cov(
                bam, contig_name, end - te_interval_size - te_offset, end - te_offset
            )
        else:
            whole_te_locus_cov = True
    else:
        whole_te_locus_cov = True

    if whole_te_locus_cov:
        te_5p_cov = get_median_cov(bam, contig_name, start, end)
        te_3p_cov = te_5p_cov
    return te_5p_cov, te_3p_cov


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

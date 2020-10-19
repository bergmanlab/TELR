import sys
import os
import subprocess
import time
import logging
from Bio import SeqIO
from multiprocessing import Process, Pool
import pysam
from TELR_utility import mkdir, check_exist, format_time


def local_assembly(
    contig_dir, vcf_parsed, out, sample_name, bam, raw_reads, thread, presets, polish
):
    """Perform local assembly using reads from parsed VCF file"""

    # Prepare reads used for local assembly and polishing
    sv_reads_dir = os.path.join(out, "sv_reads")
    
    prep_assembly(
        vcf_parsed, out, sample_name, bam, raw_reads, sv_reads_dir, read_type="sv"
    )

    mkdir(contig_dir)

    if presets == "ont":
        presets_wtdbg2 = "ont"
        presets_minimap2 = "map-ont"
    else:
        presets_wtdbg2 = "rs"
        presets_minimap2 = "map-pb"

    k = 0
    asm_pa_list = []
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            # rename variant reads
            sv_reads = sv_reads_dir + "/contig" + str(k)
            sv_reads_rename = sv_reads_dir + "/" + contig_name + ".reads.fa"
            os.rename(sv_reads, sv_reads_rename)
            thread_asm = 1
            asm_pa = [
                sv_reads_rename,
                contig_dir,
                contig_name,
                thread_asm,
                presets_wtdbg2,
                presets_minimap2,
                polish,
            ]
            asm_pa_list.append(asm_pa)
            k = k + 1
    # run assembly in parallel
    logging.info("Perform local assembly of non-reference TE loci...")
    start_time = time.time()
    try:
        pool = Pool(processes=thread)
        pool.map(run_wtdbg2, asm_pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("Local assembly failed, exiting...")
        sys.exit(1)

    proc_time = time.time() - start_time
    logging.info("Local assembly finished in " + format_time(proc_time))


def prep_assembly(
    vcf_parsed, out, sample_name, bam, raw_reads, reads_dir, read_type="sv"
):
    """Prepare reads for local assembly"""
    logging.info("Prepare reads for local assembly")

    if read_type == "sv":
        # extract read IDs
        read_ids = os.path.join(out, sample_name + ".id")
        with open(vcf_parsed, "r") as input, open(read_ids, "w") as output:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                read_list = entry[8].split(",")
                for read in read_list:
                    output.write(read + "\n")
    else:
        window = 50
        samfile = pysam.AlignmentFile(bam, "rb")
        read_ids = os.path.join(out, sample_name + ".id")
        vcf_parsed_new = vcf_parsed + ".new"
        with open(vcf_parsed, "r") as input, open(read_ids, "w") as output, open(
            vcf_parsed_new, "w"
        ) as VCF:
            for line in input:
                entry = line.replace("\n", "").split("\t")

                # get sniffles read list
                read_list = entry[8].split(",")
                reads_sniffles = set(read_list)

                ins_chr = entry[0]
                ins_breakpoint = round((int(entry[1]) + int(entry[2])) / 2)
                start = ins_breakpoint - window
                end = ins_breakpoint + window
                reads = set()
                # coverage = 0
                for read in samfile.fetch(ins_chr, start, end):
                    reads.add(read.query_name)
                for read in reads:
                    output.write(read + "\n")
                # # get breakpoint coverage
                # coverage = 0
                # for pileupcolumn in samfile.pileup(
                #     ins_chr, ins_breakpoint - 20, ins_breakpoint + 20, truncate=True
                # ):
                #     if pileupcolumn.n > coverage:
                #         coverage = pileupcolumn.n

                # # compare sniffles reads
                # num_sniffles_unique = len(reads_sniffles.difference(reads))
                # out_line = (
                #     line.replace("\n", "")
                #     + "\t"
                #     + str(len(reads))
                #     + "\t"
                #     + str(num_sniffles_unique)
                # )
                # VCF.write(out_line + "\n")

                # write
                out_line = line.replace("\n", "") + "\t" + str(len(reads))
                VCF.write(out_line + "\n")
                vcf_parsed = vcf_parsed_new

    # generate unique ID list
    read_ids_unique = read_ids + ".unique"
    command = "cat " + read_ids + " | sort | uniq"
    with open(read_ids_unique, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # filter raw reads using read list
    subset_fa = os.path.join(out, sample_name + ".subset.fa")
    command = "seqtk subseq " + raw_reads + " " + read_ids_unique + " | seqtk seq -a"
    with open(subset_fa, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # reorder reads
    subset_fa_reorder = out + "/" + sample_name + ".subset.reorder.fa"
    extract_reads(subset_fa, read_ids, subset_fa_reorder)

    # separate reads into multiple files, using csplit
    mkdir(reads_dir)
    csplit_prefix = reads_dir + "/contig"
    m = []
    k = 1
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if read_type == "sv":
                k = k + 2 * (len(entry[8].split(",")))
            else:
                k = k + 2 * int(entry[13])
            m.append(k)
    if len(m) == 1:
        subprocess.call(["cp", subset_fa_reorder, reads_dir + "/contig0"])
    elif len(m) == 0:
        print("No insertion detected, exiting...")
    else:
        m = m[:-1]
        index = " ".join(str(i) for i in m)
        command = (
            "csplit -s -f " + csplit_prefix + " -n 1 " + subset_fa_reorder + " " + index
        )
        subprocess.call(command, shell=True)

    # remove tmp files
    os.remove(read_ids)
    os.remove(read_ids_unique)
    os.remove(subset_fa)
    os.remove(subset_fa_reorder)


def run_wtdbg2(args):
    sv_reads = args[0]
    asm_dir = args[1]
    contig_name = args[2]
    thread = args[3]
    presets_wtdbg2 = args[4]
    presets_minimap2 = args[5]
    polish = args[6]

    prefix = sv_reads.replace(".reads.fa", "")
    try:
        subprocess.run(
            [
                "wtdbg2",
                "-x",
                presets_wtdbg2,
                "-q",
                "-AS",
                "1",
                "-g",
                "30k",
                "-t",
                str(thread),
                "-i",
                sv_reads,
                "-fo",
                prefix,
            ],
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        print("fail to build contig layout for contig: " + contig_name)
        return

    contig_layout = prefix + ".ctg.lay.gz"
    if check_exist(contig_layout):
        # derive consensus
        cns_thread = str(min(thread, 4))
        consensus = prefix + ".raw.fa"
        try:
            subprocess.run(
                [
                    "wtpoa-cns",
                    "-q",
                    "-t",
                    cns_thread,
                    "-i",
                    contig_layout,
                    "-fo",
                    consensus,
                ],
                timeout=300,
            )
        except subprocess.TimeoutExpired:
            print("fail to assemble contig: " + contig_name)
            return
        if check_exist(consensus):
            consensus_final = asm_dir + "/" + contig_name + ".cns.fa"
            if polish > 0:
                polish_contig(
                    prefix,
                    sv_reads,
                    consensus,
                    consensus_final,
                    thread,
                    presets_minimap2,
                    polish,
                    contig_name,
                )
            else:
                # move raw consensus to asm dir
                os.rename(consensus, consensus_final)
        else:
            print("Initial assembly fail for " + prefix + "\n")
    else:
        print("Build contig layout fail for " + prefix + "\n")


def polish_contig(
    prefix, reads, raw_contig, polished_contig, thread, preset, round, contig_name
):
    # polish consensus
    polish_thread = str(min(thread, 4))
    bam = raw_contig + ".bam"
    k = 0
    while True:
        tmp_contig = raw_contig + ".tmp"
        command = (
            "minimap2 -t "
            + polish_thread
            + " -ax "
            + preset
            + " -r2k "
            + raw_contig
            + " "
            + reads
            + " | samtools sort -@"
            + polish_thread
            + " > "
            + bam
        )
        try:
            subprocess.run(
                command,
                shell=True,
                timeout=300,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except subprocess.TimeoutExpired:
            print("fail to map reads to contig: " + contig_name)
            return
        command = (
            "samtools view -F0x900 "
            + bam
            + " | wtpoa-cns -t "
            + polish_thread
            + " -d "
            + raw_contig
            + " -i - -fo "
            + tmp_contig
        )
        try:
            subprocess.run(
                command,
                shell=True,
                timeout=300,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except subprocess.TimeoutExpired:
            print("fail to polish contig: " + contig_name)
            return
        raw_contig = tmp_contig
        k = k + 1
        if k >= round:
            break
    if os.path.isfile(tmp_contig) and os.stat(tmp_contig).st_size != 0:
        os.rename(tmp_contig, polished_contig)
    else:
        print("polishing failed for " + contig_name + "\n")
    return


def extract_reads(reads, list, out):
    """Extract reads from fasta using read ID list"""
    record_dict = SeqIO.index(reads, "fasta")
    with open(out, "wb") as output_handle, open(list, "r") as ID:
        for entry in ID:
            entry = entry.replace("\n", "")
            output_handle.write(record_dict.get_raw(entry))

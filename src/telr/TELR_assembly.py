import sys
import os
import subprocess
import shutil
import time
import logging
from Bio import SeqIO
from multiprocessing import Pool
import pysam
from telr.TELR_utility import mkdir, check_exist, format_time


def get_local_contigs(
    assembler,
    polisher,
    contig_dir,
    vcf_parsed,
    out,
    sample_name,
    bam,
    raw_reads,
    thread,
    presets,
    polish_iterations,
):
    """Perform local assembly using reads from parsed VCF file in parallel"""

    # Prepare reads used for local assembly and polishing
    sv_reads_dir = os.path.join(out, "sv_reads")

    try:
        prep_assembly_inputs(
            vcf_parsed, out, sample_name, bam, raw_reads, sv_reads_dir, read_type="sv"
        )
    except Exception as e:
        print(e)
        print("Prepare local assembly input data failed, exiting...")
        sys.exit(1)

    mkdir(contig_dir)

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
                presets,
                assembler,
                polisher,
                polish_iterations,
            ]
            asm_pa_list.append(asm_pa)
            k = k + 1
    # run assembly in parallel
    logging.info("Perform local assembly of non-reference TE loci...")
    start_time = time.time()

    try:
        pool = Pool(processes=thread)
        contig_list = pool.map(run_assembly_polishing, asm_pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("Local assembly failed, exiting...")
        sys.exit(1)

    proc_time = time.time() - start_time

    # merge all contigs
    assembly_passed_loci = set()
    merged_contigs = os.path.join(out, sample_name + ".contigs.fa")
    with open(merged_contigs, "w") as merged_output_handle:
        for contig in contig_list:
            if check_exist(contig):
                contig_name = os.path.basename(contig).replace(".cns.fa", "")
                assembly_passed_loci.add(contig_name)
                parsed_contig = os.path.join(contig_dir, contig_name + ".cns.ctg1.fa")
                with open(contig, "r") as input:
                    records = SeqIO.parse(input, "fasta")
                    for record in records:
                        if record.id == "ctg1" or record.id == "contig_1":
                            record.id = contig_name
                            record.description = "len=" + str(len(record.seq))
                            SeqIO.write(record, merged_output_handle, "fasta")
                            with open(parsed_contig, "w") as parsed_output_handle:
                                SeqIO.write(record, parsed_output_handle, "fasta")

    logging.info("Local assembly finished in " + format_time(proc_time))
    return merged_contigs, assembly_passed_loci


def run_assembly_polishing(args):
    reads = args[0]
    asm_dir = args[1]
    contig_name = args[2]
    thread = args[3]
    presets = args[4]
    assembler = args[5]
    polisher = args[6]
    polish_iterations = args[7]

    # run assembly
    if assembler == "wtdbg2":
        asm_cns = run_wtdbg2_assembly(reads, asm_dir, contig_name, thread, presets)
    else:
        asm_cns = run_flye_assembly(reads, asm_dir, contig_name, thread, presets)

    if not check_exist(asm_cns):
        print("assembly failed")
        return None

    # run polishing
    if polish_iterations > 0:
        if polisher == "wtdbg2":
            asm_cns = run_wtdbg2_polishing(
                asm_cns, reads, thread, polish_iterations, presets
            )
        else:
            asm_cns = run_flye_polishing(
                asm_cns, reads, asm_dir, contig_name, thread, polish_iterations, presets
            )

    if check_exist(asm_cns):
        return asm_cns
    else:
        return None


def run_flye_polishing(
    asm_cns, reads, asm_dir, contig_name, thread, polish_iterations, presets
):
    """Run Flye polishing"""
    if presets == "pacbio":
        presets_flye = "--pacbio-raw"
    else:
        presets_flye = "--nano-raw"

    tmp_out_dir = os.path.join(asm_dir, contig_name)
    mkdir(tmp_out_dir)
    try:
        subprocess.call(
            [
                "flye",
                "--polish-target",
                asm_cns,
                presets_flye,
                reads,
                "--out-dir",
                tmp_out_dir,
                "--thread",
                str(thread),
                "--iterations",
                str(polish_iterations),
            ]
        )
    except Exception as e:
        print(e)
        print("Polishing failed, exiting...")
        return None

    # rename contig file
    polished_contig = os.path.join(
        tmp_out_dir, "polished_" + str(polish_iterations) + ".fasta"
    )
    if check_exist(polished_contig):
        os.rename(polished_contig, asm_cns)
        shutil.rmtree(tmp_out_dir)
        return asm_cns
    else:
        return None


def run_wtdbg2_polishing(asm_cns, reads, threads, polish_iterations, presets):
    """Run wtdbg2 polishing"""
    if presets == "pacbio":
        presets_minimap2 = "map-pb"
    else:
        presets_minimap2 = "map-ont"

    # polish consensus
    threads = str(min(threads, 4))
    bam = asm_cns + ".bam"
    k = 0
    while True:
        # align reads to contigs

        command = (
            "minimap2 -t "
            + threads
            + " -ax "
            + presets_minimap2
            + " -r2k "
            + asm_cns
            + " "
            + reads
            + " | samtools sort -@"
            + threads
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
            print("fail to map reads to contig: " + asm_cns)
            return

        # run wtpoa-cns to get polished contig
        cns_tmp = asm_cns + ".tmp"
        command = (
            "samtools view -F0x900 "
            + bam
            + " | wtpoa-cns -t "
            + threads
            + " -d "
            + asm_cns
            + " -i - -fo "
            + cns_tmp
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
            print("fail to polish contig: " + asm_cns)
            return
        if check_exist(cns_tmp):
            os.rename(cns_tmp, asm_cns)
            os.remove(bam)
        else:
            break
        k = k + 1
        if k >= polish_iterations:
            break
    if check_exist(asm_cns):
        return asm_cns
    else:
        print("polishing failed for " + asm_cns + "\n")
    return None


def run_flye_assembly(sv_reads, asm_dir, contig_name, thread, presets):
    """Run Flye assembly"""
    if presets == "pacbio":
        presets_flye = "--pacbio-raw"
    else:
        presets_flye = "--nano-raw"

    tmp_out_dir = os.path.join(asm_dir, contig_name)
    mkdir(tmp_out_dir)
    try:
        subprocess.call(
            [
                "flye",
                presets_flye,
                sv_reads,
                "--out-dir",
                tmp_out_dir,
                "--thread",
                str(thread),
                "--iterations",
                "0",
            ]
        )
    except Exception as e:
        print(e)
        print("Assembly failed, exiting...")
        return
    # rename contigs
    contig_path = os.path.join(tmp_out_dir, "assembly.fasta")
    contig_path_new = os.path.join(asm_dir, contig_name + ".cns.fa")
    if check_exist(contig_path):
        os.rename(contig_path, contig_path_new)
        # remove tmp files
        shutil.rmtree(tmp_out_dir)
        return contig_path_new
    else:
        print("assembly failed")
        return None


def run_wtdbg2_assembly(sv_reads, asm_dir, contig_name, thread, presets):
    """Run wtdbg2 assembly"""
    if presets == "pacbio":
        presets_wtdbg2 = "rs"
    else:
        presets_wtdbg2 = "ont"
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
    except Exception as e:
        print(e)
        print("wtdbg2 failed, exiting...")
        return None

    # derive consensus
    contig_layout = prefix + ".ctg.lay.gz"
    if check_exist(contig_layout):
        cns_thread = str(min(thread, 4))
        consensus = prefix + ".cns.fa"
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
            return None

    if check_exist(consensus):
        consensus_rename = os.path.join(asm_dir, contig_name + ".cns.fa")
        os.rename(consensus, consensus_rename)
        return consensus_rename
    else:
        return None


def prep_assembly_inputs(
    vcf_parsed, out, sample_name, bam, raw_reads, reads_dir, read_type="sv"
):
    """Prepare reads for local assembly"""
    # logging.info("Prepare reads for local assembly")

    if read_type == "sv":  # TODO: figure out what this does
        # extract read IDs
        read_ids = os.path.join(out, f"{sample_name}.id")
        with open(vcf_parsed, "r") as input, open(read_ids, "w") as output:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                read_list = entry[8].split(",")
                for read in read_list:
                    output.write(read + "\n")
    else:  # TODO: think about using this for assembly, filter for cigar reads
        window = 1000
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
                if start < 0:
                    start = 0
                end = ins_breakpoint + window
                reads = set()
                # coverage = 0
                for read in samfile.fetch(ins_chr, start, end):
                    reads.add(read.query_name)
                for read in reads:
                    output.write(read + "\n")

                # write
                out_line = line.replace("\n", "") + "\t" + str(len(reads))
                VCF.write(out_line + "\n")
                vcf_parsed = vcf_parsed_new

    # generate unique ID list
    read_ids_unique = f"{read_ids}.unique"
    command = f"cat {read_ids} | sort | uniq"
    with open(read_ids_unique, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # filter raw reads using read list
    subset_fa = os.path.join(out, f"{sample_name}.subset.fa")
    command = f"seqtk subseq {raw_reads} {read_ids_unique} | seqtk seq -a"
    with open(subset_fa, "w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # reorder reads
    subset_fa_reorder = f"{out}/{sample_name}.subset.reorder.fa"
    extract_reads(subset_fa, read_ids, subset_fa_reorder)

    # separate reads into multiple files, using csplit
    mkdir(reads_dir)
    csplit_prefix = f"{reads_dir}/contig"
    m = []
    k = 1
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if read_type == "sv":
                k += 2 * (len(entry[8].split(",")))
            else:
                k += 2 * int(entry[14])
            m.append(k)
    if len(m) == 1:
        subprocess.call(["cp", subset_fa_reorder, f"{reads_dir}/contig0"])
    elif len(m) == 0:
        print("No insertion detected, exiting...")
    else:
        m = m[:-1]
        index = " ".join(str(i) for i in m)
        command = f"csplit -s -f {csplit_prefix} -n 1 {subset_fa_reorder} {index}"
        subprocess.call(command, shell=True)

    # remove tmp files
    os.remove(read_ids)
    os.remove(read_ids_unique)
    os.remove(subset_fa)
    os.remove(subset_fa_reorder)


def extract_reads(reads, list, out):
    """Extract reads from fasta using read ID list"""
    record_dict = SeqIO.index(reads, "fasta") #SeqIO.index() pulls reads in a memory efficient way/without opening the whole file.
    with open(out, "wb") as output_handle, open(list, "r") as ID:
        for entry in ID:
            entry = entry.replace("\n", "")
            output_handle.write(record_dict.get_raw(entry))

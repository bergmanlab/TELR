import sys
import os
import subprocess
import shutil
import time
import logging
from Bio import SeqIO
from multiprocessing import Pool
import pysam
from TELR_utility import mkdir, check_exist, format_time, get_contig_name


def parse_assembled_contig(input_contig, parsed_contig):
    if check_exist(input_contig):
        contig_name = os.path.basename(input_contig).replace(".cns.fa", "")
        with open(input_contig, "r") as input, open(parsed_contig, "w") as parsed_output_handle:
            records = SeqIO.parse(input, "fasta")
            for record in records:
                if record.id == "ctg1" or record.id == "contig_1":
                    record.id = contig_name
                    record.description = "len=" + str(len(record.seq))
                    SeqIO.write(record, parsed_output_handle, "fasta")

def run_flye_polishing(initial_assembly, assembled_consensus, reads, asm_dir, contig_name, thread, polish_iterations, presets):
    """Run Flye polishing"""
    presets_flye = {"pacbio":"--pacbio-raw","ont":"--nano-raw"}[presets]
    tmp_out_dir = os.path.join(asm_dir, contig_name)
    mkdir(tmp_out_dir)
    try:
        subprocess.call(
            [
                "flye",
                "--polish-target",
                initial_assembly,
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
        sys.exit(1)

    # rename contig file
    polished_contig = os.path.join(tmp_out_dir, f"polished_{polish_iterations}.fasta")
    if check_exist(polished_contig):
        os.rename(polished_contig, assembled_consensus)
        shutil.rmtree(tmp_out_dir)
        sys.exit(0)
    else:
        sys.exit(1)


def run_wtdbg2_polishing(initial_assembly, assembled_consensus, reads, asm_dir, contig_name, threads, polish_iterations, presets):
    """Run wtdbg2 polishing"""
    presets_minimap2 = {"pacbio":"map-pb","ont":"map-ont"}[presets]
    polish_iterations = int(polish_iterations)
    # polish consensus
    threads = str(min(int(threads), 4))
    bam = f"{initial_assembly}.bam"
    iteration = 0
    while True:
        # align reads to contigs
        command = f"minimap2 -t {threads} -ax {presets_minimap2} -r2k {initial_assembly} {reads} | samtools sort -@{threads} > {bam}"
        try:
            subprocess.run(
                command,
                shell=True,
                timeout=300,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except subprocess.TimeoutExpired:
            print(f"fail to map reads to contig: {initial_assembly}")
            sys.exit(1)

        # run wtpoa-cns to get polished contig
        cns_tmp = f"{initial_assembly}.tmp"
        command = f"samtools view -F0x900 {bam} | wtpoa-cns -t {threads} -d {initial_assembly} -i - -fo {cns_tmp}"
        try:
            subprocess.run(
                command,
                shell=True,
                timeout=300,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except subprocess.TimeoutExpired:
            print(f"fail to polish contig: {initial_assembly}")
            sys.exit(1)
        if check_exist(cns_tmp):
            os.rename(cns_tmp, assembled_consensus)
            os.remove(bam)
        else:
            break
        iteration = iteration + 1
        if iteration >= polish_iterations:
            break
    if check_exist(assembled_consensus):
        sys.exit(0)
    else:
        print(f"polishing failed for {initial_assembly}\n")
        sys.exit(1)


def run_flye_assembly(sv_reads, asm_dir, contig_name, thread, presets, output):
    mkdir(asm_dir, verbose = False)
    """Run Flye assembly"""
    presets_flye = {"pacbio":"--pacbio-raw","ont":"--nano-raw"}[presets]

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
        sys.exit(1)
    # rename contigs
    contig_path = os.path.join(tmp_out_dir, "assembly.fasta")
    if check_exist(contig_path):
        os.rename(contig_path, output)
        # remove tmp files
        shutil.rmtree(tmp_out_dir)
        sys.exit(0)
    else:
        print("assembly failed")
        sys.exit(1)


def run_wtdbg2_assembly(sv_reads, asm_dir, contig_name, thread, presets, output):
    mkdir(asm_dir, verbose = False)
    """Run wtdbg2 assembly"""
    presets_wtdbg2 = {"pacbio":"rs","ont":"ont"}[presets]
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
        print(f"fail to build contig layout for contig: {contig_name}")
        sys.exit(1)
    except Exception as e:
        print(e)
        print("wtdbg2 failed, exiting...")
        sys.exit(1)

    # derive consensus
    contig_layout = f"{prefix}.ctg.lay.gz"
    if check_exist(contig_layout):
        cns_thread = str(min(int(thread), 4))
        consensus = f"{prefix}.cns.fa"
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
            print(f"fail to assemble contig: {contig_name}")
            sys.exit(1)

    if check_exist(consensus):
        os.rename(consensus, output)
        sys.exit(0)
    else:
        sys.exit(1)

def write_read_IDs(vcf_parsed, read_ids):
    with open(vcf_parsed, "r") as input, open(read_ids, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            read_list = entry[8].split(",")
            for read in read_list:
                output.write(read + "\n")

def separate_read_files(reads_dir, vcf_parsed, subset_fa_reorder):
    mkdir(reads_dir)
    csplit_prefix = reads_dir + "/contig"
    m = []
    k = 1
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            k = k + 2 * (len(entry[8].split(",")))
            m.append(k)
    if len(m) == 1:
        subprocess.call(["cp", subset_fa_reorder, reads_dir + "/contig0"])
    elif len(m) == 0:
        print("No insertion detected, exiting...")
    else:
        m = m[:-1]
        index = " ".join(str(i) for i in m)
        command = f"csplit -s -f {csplit_prefix} -n 1 {subset_fa_reorder} {index}"
        subprocess.call(command, shell=True)

def extract_reads(reads, id_list, out):
    """Extract reads from fasta using read ID list"""
    record_dict = SeqIO.index(reads, "fasta")
    if type(id_list) in [list, set]:
        ID = id_list
    else:
        ID = open(id_list, "r")
    with open(out, "wb") as output_handle:
        for entry in ID:
            entry = entry.replace("\n", "")
            output_handle.write(record_dict.get_raw(entry))

def make_contig_file(vcf_parsed, contig_name, contig_file, reads):
    mkdir(os.path.split(contig_file)[0], verbose = False)
    contig_name = contig_name.split("_")
    contig_name = f"{'_'.join(contig_name[:-2])}\t{contig_name[-2]}\t{contig_name[-1]}"
    with open(vcf_parsed, "r") as input:
        sequences = [i for i in input if contig_name in i]
        sequences = set(sequences[0].split("\t")[8].split(",")) # get all of the read names for the matching contig
    extract_reads(reads, sequences, contig_file)



if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
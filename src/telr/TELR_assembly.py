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
    files,
    assembler,
    polisher,
    out,
    thread,
    presets,
    polish_iterations,
):
    """Perform local assembly using reads from parsed VCF file in parallel"""

    # Prepare reads used for local assembly and polishing
    files.__dict__[out].dir("sv_reads_dir", "sv_reads")

    try:
        prep_assembly_inputs(files, out)
    except Exception as e:
        print(e)
        print("Prepare local assembly input data failed, exiting...")
        sys.exit(1)

    files.contig.make()

    k = 0
    asm_pa_list = []
    with files.vcf_parsed.open() as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            # rename variant reads
            sv_reads = f"{files.sv_reads_dir.path}/{contig_name}.reads.fa"
            os.rename(f"{files.sv_reads_dir.path}/contig{k}", sv_reads)
            thread_asm = 1
            asm_pa = [
                sv_reads,
                files.contig.path,
                contig_name,
                thread_asm,
                presets,
                assembler,
                polisher,
                polish_iterations,
            ]
            asm_pa_list.append(asm_pa)
            k += 1
    # run assembly in parallel
    logging.info("Perform local assembly of non-reference TE loci...")
    start_time = time.perf_counter()

    try:
        pool = Pool(processes=thread)
        contig_list = pool.map(run_assembly_polishing, asm_pa_list)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("Local assembly failed, exiting...")
        sys.exit(1)

    proc_time = time.perf_counter() - start_time

    # merge all contigs
    files.set("assembly_passed_loci")
    files.set("assembly_passed_loci_merged")
    files.add("merged_contigs",out, ".contigs.fa")
    with files.merged_contigs.open("w") as merged_output_handle:
        for contig in contig_list:
            if check_exist(contig):
                contig_name = os.path.basename(contig).replace(".cns.fa", "")
                files.set_file("assembly_parsed_loci","sv_reads_dir",contig_name,".cns.fa")
                parsed_contig = files.set_file("assembly_parsed_loci_merged","contig",contig_name,".cns.ctg1.fa")
                with open(contig, "r") as input: #contig == output from assembly/polishing process
                    records = SeqIO.parse(input, "fasta")
                    for record in records:
                        if record.id == "ctg1" or record.id == "contig_1":
                            record.id = contig_name
                            record.description = f"len={len(record.seq)}"
                            SeqIO.write(record, merged_output_handle, "fasta")
                            with parsed_contig.open("w") as parsed_output_handle:
                                SeqIO.write(record, parsed_output_handle, "fasta")

    logging.info(f"Local assembly finished in {format_time(proc_time)}")


def run_assembly_polishing(args):
    '''Function called from multiprocessing function for threaded assembly. 
    Makes decisions about other functions within TELR_assembly.py to call 
    on the data based on input parameters, and calls them.'''
    reads = args[0]
    asm_dir = args[1]
    contig_name = args[2]
    thread = args[3]
    presets = args[4]
    assembler = {"wtdbg2":run_wtdbg2_assembly,"flye":run_flye_assembly}[args[5]] #potential to throw error if invalid input
    polisher = {"wtdbg2":run_wtdbg2_polishing,"flye":run_flye_polishing}[args[6]]
    polish_iterations = args[7]

    # run assembly
    asm_cns = assembler(reads, asm_dir, contig_name, thread, presets)

    if not check_exist(asm_cns):
        print("assembly failed")
        return None

    # run polishing
    if polish_iterations > 0:
        asm_cns = polisher(asm_cns, reads, asm_dir, contig_name, thread, polish_iterations, presets)

    if check_exist(asm_cns):
        return asm_cns
    else:
        return None


def run_flye_polishing(asm_cns, reads, asm_dir, contig_name, thread, polish_iterations, presets):
    """Run Flye polishing"""
    presets_flye = {"pacbio":"--pacbio-raw","ont":"--nano-raw"}[presets]
    
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
                str(polish_iterations)
            ]
        )
    except Exception as e:
        print(e)
        print("Polishing failed, exiting...")
        return None

    # rename contig file
    polished_contig = os.path.join(tmp_out_dir, f"polished_{polish_iterations}.fasta")
    if check_exist(polished_contig):
        os.rename(polished_contig, asm_cns)
        shutil.rmtree(tmp_out_dir)
        return asm_cns
    else:
        return None


def run_wtdbg2_polishing(asm_cns, reads, asm_dir, contig_name, threads, polish_iterations, presets):
    """Run wtdbg2 polishing"""
    presets_minimap2 = {"pacbio":"map-pb","ont":"map-ont"}[presets]

    # polish consensus
    threads = str(min(threads, 4))
    bam = f"{asm_cns}.bam"
    k = 0
    while True:
        # align reads to contigs
        command = f"minimap2 -t {threads} -ax {presets_minimap2} -r2k {asm_cns} {reads} | samtools sort -@{threads} > {bam}"
        try:
            subprocess.run(
                command,
                shell=True,
                timeout=300,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except subprocess.TimeoutExpired:
            print(f"fail to map reads to contig: {asm_cns}")
            return

        # run wtpoa-cns to get polished contig
        cns_tmp = f"{asm_cns}.tmp"
        command = f"samtools view -F0x900 {bam} | wtpoa-cns -t {threads} -d {asm_cns} -i - -fo {cns_tmp}"
        try:
            subprocess.run(
                command,
                shell=True,
                timeout=300,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except subprocess.TimeoutExpired:
            print(f"fail to polish contig: {asm_cns}")
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
        print(f"polishing failed for {asm_cns}\n")
    return None


def run_flye_assembly(sv_reads, asm_dir, contig_name, thread, presets):
    """Run Flye assembly"""
    tmp_out_dir = os.path.join(asm_dir, contig_name)
    mkdir(tmp_out_dir)
    try:
        subprocess.call(
            [
                "flye",
                {"pacbio":"--pacbio-raw","ont":"--nano-raw"}[presets],
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
    contig_path_new = os.path.join(asm_dir, f"{contig_name}.cns.fa")
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
    prefix = sv_reads.replace(".reads.fa", "")
    try:
        subprocess.run(
            [
                "wtdbg2",
                "-x",
                {"pacbio":"rs","ont":"ont"}[presets],
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
                prefix
            ],
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        print(f"fail to build contig layout for contig: {contig_name}")
        return
    except Exception as e:
        print(e)
        print("wtdbg2 failed, exiting...")
        return None

    # derive consensus
    contig_layout = f"{prefix}.ctg.lay.gz"
    if check_exist(contig_layout):
        consensus = f"{prefix}.cns.fa"
        try:
            subprocess.run(
                [
                    "wtpoa-cns",
                    "-q",
                    "-t",
                    str(min(thread, 4)),
                    "-i",
                    contig_layout,
                    "-fo",
                    consensus,
                ],
                timeout=300,
            )
        except subprocess.TimeoutExpired:
            print(f"fail to assemble contig: {contig_name}")
            return None

    if check_exist(consensus):
        consensus_rename = os.path.join(asm_dir, f"{contig_name}.cns.fa")
        os.rename(consensus, consensus_rename)
        return consensus_rename
    else:
        return None


def prep_assembly_inputs(
    files, out, read_type="sv"
):
    """Prepare reads for local assembly"""
    # logging.info("Prepare reads for local assembly")

    if read_type == "sv":  # TODO: figure out what this does
        # extract read IDs
        files.add("read_ids",out,".id")
        with files.vcf_parsed.open() as input, files.read_ids.open("w") as output:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                read_list = entry[8].split(",")
                for read in read_list:
                    output.write(read + "\n")
    else:  # TODO: think about using this for assembly, filter for cigar reads
        window = 1000
        samfile = pysam.AlignmentFile(files.bam.path, "rb")
        files.add("read_ids",out,".id")
        files.vcf_parsed.extend("vcf_parsed_new",".new", file_format = "vcf")
        with files.vcf_parsed.open() as input, files.read_ids.open("w") as output, files.vcf_parsed_new.open("w") as VCF:
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
                files.vcf_parsed_new.rename("vcf_parsed")

    # generate unique ID list
    files.read_ids.extend("read_ids_unique",".unique",file_format="txt")
    command = f"cat {files.read_ids.path} | sort | uniq"
    with files.read_ids_unique.open("w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # filter raw reads using read list
    files.add("subset_fa",out,".subset.fa")
    command = f"seqtk subseq {files.reads.path} {files.read_ids_unique.path} | seqtk seq -a"
    with files.subset_fa.open("w") as output:
        subprocess.call(command, stdout=output, shell=True)

    # reorder reads
    files.add("subset_fa_reorder",out,".subset.reorder.fa")
    extract_reads(files.subset_fa, files.read_ids, files.subset_fa_reorder)

    # separate reads into multiple files, using csplit
    files.mkdir("sv_reads_dir")
    csplit_indices = []
    k = 1
    with files.vcf_parsed.open() as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if read_type == "sv":
                k += 2 * (len(entry[8].split(",")))
            else:
                k += 2 * int(entry[14])
            csplit_indices.append(k)
    if len(csplit_indices) == 1:
        subprocess.call(["cp", files.subset_fa_reorder.path, files.directories["sv_reads_dir"] + "/contig0"])
    elif len(csplit_indices) == 0:
        print("No insertion detected, exiting...")
    else:
        index = " ".join(str(i) for i in csplit_indices[:-1])
        command = "csplit -s -f " + files.directories["sv_reads_dir"] + f"/contig -n 1 {files.subset_fa_reorder.path} {index}"
        subprocess.call(command, shell=True)

    # remove tmp files
    files.read_ids.remove()
    files.read_ids_unique.remove()
    files.subset_fa.remove()
    files.subset_fa_reorder.remove()


def extract_reads(reads, list, out):
    """Extract reads from fasta using read ID list"""
    record_dict = SeqIO.index(reads.path, "fasta") #SeqIO.index() pulls reads in a memory efficient way/without opening the whole file.
    with out.open("wb") as output_handle, list.open() as ID:
        for entry in ID:
            entry = entry.replace("\n", "")
            output_handle.write(record_dict.get_raw(entry))

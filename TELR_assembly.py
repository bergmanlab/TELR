import sys
import os
import subprocess
from multiprocessing import Process, Pool
from TELR_utility import mkdir, check_exist

def local_assembly(contig_reads_dir, vcf, out, raw_reads, thread, presets, polish):
    contig_assembly_dir=out+"/"+"contig_assembly"
    mkdir(contig_assembly_dir)

    if presets == "ont":
        presets_wtdbg2 = "ont"
        presets_minimap2 = "map-ont"
    else:
        presets_wtdbg2 = "rs"
        presets_minimap2 = "map-pb"

    print ("Assemble contigs...")
    k = 0
    asm_pa_list=[]
    with open(vcf, "r") as input:
        for line in input:
            entry = line.replace('\n', '').split("\t")
            contig_name = "_".join([entry[0], entry[1], entry[2]])
            contig_reads=contig_reads_dir + "/contig" + str(k)
            # rename contig reads
            contig_reads_rename=contig_reads_dir + "/" + contig_name + ".reads.fa"
            os.rename(contig_reads, contig_reads_rename)
            thread_asm=1
            asm_pa = [contig_reads_rename, contig_assembly_dir, contig_name, thread_asm, presets_wtdbg2, presets_minimap2, polish]
            asm_pa_list.append(asm_pa)
            k = k + 1
    # run assembly in parallel
    pool = Pool(processes=thread)
    pool.map(run_wtdbg2, asm_pa_list)
    pool.close()
    pool.join()
    print ("Done\n")
    return(contig_assembly_dir)

def run_wtdbg2(args):
    reads = args[0]
    asm_dir = args[1]
    contig_name = args[2]
    thread = args[3]
    presets_wtdbg2 = args[4]
    presets_minimap2 = args[5]
    polish = args[6]

    prefix = reads.replace('.reads.fa', '')
    command = "wtdbg2 -x "+presets_wtdbg2+" -q -AS 1 -g 30k -t "+str(thread)+" -i "+reads+" -fo "+prefix
    try:
        subprocess.run(command, shell = True, timeout = 300)
    except subprocess.TimeoutExpired:
        print("fail to build contig layout for contig: " + contig_name)
        return

    contig_layout=prefix+".ctg.lay.gz"
    if check_exist(contig_layout):
        # derive consensus
        cns_thread = str(min(thread, 4))
        contig_raw=prefix+".ctg.lay.gz"
        consensus = prefix + ".raw.fa"
        command = "wtpoa-cns -q -t "+cns_thread+" -i "+contig_layout+" -fo "+consensus
        try:
            subprocess.run(command, shell = True, timeout = 300)
        except subprocess.TimeoutExpired:
            print("fail to assemble contig: " + contig_name)
            return
        if check_exist(consensus):
            consensus_final = asm_dir + "/" + contig_name + ".cns.fa"
            if polish > 0:
                polish_contig(prefix, reads, consensus, consensus_final, thread, presets_minimap2, polish, contig_name)
            else:
                # move raw consensus to asm dir
                os.rename(consensus, consensus_final)
        else:
            print("Initial assembly fail for "+prefix+"\n")
    else:
        print("Build contig layout fail for "+prefix+"\n")

def polish_contig(prefix, reads, raw_contig, polished_contig, thread, preset, round, contig_name):
    # polish consensus
    polish_thread = str(min(thread, 4))
    bam = raw_contig + ".bam"
    k = 0
    while True:
        tmp_contig = raw_contig + ".tmp"
        command = "minimap2 -t " + polish_thread + " -ax " + preset + " -r2k " + raw_contig + " " + reads + " | samtools sort -@" + polish_thread + " > " + bam
        try:
            subprocess.run(command, shell = True, timeout = 300, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
        except subprocess.TimeoutExpired:
            print("fail to map reads to contig: " + contig_name)
            return
        command = "samtools view -F0x900 " + bam + " | wtpoa-cns -t " + polish_thread + " -d " + raw_contig + " -i - -fo " + tmp_contig
        try:
            subprocess.run(command, shell = True, timeout = 300, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
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
        print("polishing failed for "+contig_name+"\n")
    return
    
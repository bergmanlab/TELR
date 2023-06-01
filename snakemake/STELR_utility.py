import os
import sys
from datetime import datetime, timedelta
import subprocess
from Bio import SeqIO


def rm_file(file):
    if os.path.exists(file):
        os.remove(file)


def get_contig_name(some_input):
    if type(some_input) is list:
        return "_".join([some_input[0], some_input[1], some_input[2]])
    elif type(some_input) is str:
        if "contigs/" in some_input:
            some_input = some_input[some_input.index("contigs/")+8:]
            return some_input[:some_input.index("/")] 


def mkdir(dir, verbose = True):
    if os.path.isdir(dir):
        if verbose:
            print(f"Directory {dir} exists")
        return
    try:
        os.mkdir(dir)
    except OSError:
        print(f"Creation of the directory {dir} failed")
    else:
        if verbose:
            print(f"Successfully created the directory {dir}")


def check_exist(file):
    if file:
        if os.path.isfile(file) and os.stat(file).st_size != 0:
            return True
        else:
            return False
    else:
        return False


def format_time(time):
    d = datetime(1, 1, 1) + timedelta(seconds=time)
    if d.hour == 0 and d.minute == 0:
        return f"{d.second} seconds"
    elif d.hour == 0 and d.minute != 0:
        return f"{d.minute} minutes {d.second} seconds"
    else:
        return f"{d.hour} hours {d.minute} minutes {d.second} seconds"


def create_loci_set(vcf_parsed):
    all_loci = set()
    with open(vcf_parsed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            all_loci.add("_".join(entry[0:3]))
    return all_loci


def report_bad_loci(set_raw, set_filtered, report, message):
    with open(report, "a") as output:
        for locus in set_raw:
            if locus not in set_filtered:
                output.write("\t".join([locus, message]) + "\n")


def get_cmd_output(cmd_list):
    """get output from subprocess"""
    output = subprocess.check_output(cmd_list)
    output = output.decode("utf-8")
    return output


def get_rev_comp_sequence(fasta_in, fasta_out):
    """get reverse complement of a sequence"""
    with open(fasta_out, "w") as output:
        for record in SeqIO.parse(fasta_in, "fasta"):
            output.write(f">{record.id}\n{record.seq.reverse_complement()}\n")

def export_env(file):
    """export conda environment"""
    file_tmp = file + ".tmp"
    cmd_list = ["conda", "env", "export", "--name", "TELR", "--file", file_tmp]
    subprocess.call(cmd_list)
    with open(file, "w") as output, open(file_tmp, "r") as input:
        for line in input:
            if (
                not line.startswith("prefix:")
                and "- pip:" not in line
                and "- telr==" not in line
            ):
                output.write(line)
    os.remove(file_tmp)

def read_vcf(vcf_file, contig_name = False, column = False, single_contig = False):
    if contig_name is False:
        pass #add functionality later if needed
    else:
        contig_name = contig_name.split("_")
        contig_name = f"{'_'.join(contig_name[:-2])}\t{contig_name[-2]}\t{contig_name[-1]}"
        with open(vcf_file, "r") as input:
            if single_contig: matching_row = input[0]
            else: matching_row = [i for i in input if contig_name in i][0]
            matching_row = matching_row.replace("\n","").split("\t")
        if type(column) is int: 
            return matching_row[column]
        elif type(column) in [list, set]:
            return {column_number:matching_row[column_number] for column_number in column}
        elif column == False:
            return matching_row

def minimap2bed(minimap_file, bed_file):#0 and 5 swapped, recheck later
    with open(minimap_file, "r") as input, open(bed_file, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            bed_line = "\t".join(
                [entry[5], entry[7], entry[8], entry[0], entry[11], entry[4]]
            )
            output.write(bed_line + "\n")

def get_contig_length(contig):
    if os.path.isfile(contig):
        with open(contig, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                contig_length = len(record.seq)
                return contig_length
    else:
        print("no contig " + contig)

def string_to_bool(bool_string, default=False)
    bool_string = str(bool_string).lower()
    bool_dict = {"true":True,"false":False}
    if bool_string in bool_dict: return bool_dict[bool_string]
    else: return default

if __name__ == '__main__':
    globals()[sys.argv[1]](*sys.argv[2:])
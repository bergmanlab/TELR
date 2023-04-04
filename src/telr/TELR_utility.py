import os
from datetime import datetime, timedelta
import subprocess
from Bio import SeqIO

class file_container:
    def __init__(self, sample_name, logging=False):
        self.sample_name = sample_name
        self.quiet = True
        if(logging):
            self.quiet = False
            self.file_log = ""
            self.frame = "main()"

    def add(self, key, directory, extension, **kwargs):
        file_name = self.sample_name
        if "new_dir" in kwargs:
            self.directories.update(kwargs.pop("new_dir"))
        if("file_name" in kwargs): file_name = kwargs.pop("file_name")
        self.prev_frame = self.frame
        if("frame" in kwargs): self.frame = kwargs.pop("frame")
        self.__dict__[key] = file(
            self, key, directory,
            f"{file_name}{extension}",
            self.quiet,
            **kwargs
        )
        self.frame = self.prev_frame
    
    def dir(self, key, path):
        self.__dict__[key] = directory(self, path)
    
    def extend(self, file, key, extension, **kwargs):
        filename = file.name
        if "file_name" in kwargs: filename = kwargs["file_name"]
        directory = file.directory
        if "new_dir" in kwargs:
            directory = kwargs.pop("new_dir")
        self.add(key, directory, extension, file_name = filename, **kwargs)
        return self.__dict__[key]

    def input_file(self, key, path, **kwargs):
        filename = os.path.basename(path)
        directory = path[:path.rindex("/")]
        if "input" in self.__dict__:
            if self.input == directory:
                directory = "input"
            else:
                self.directories[f"input_{key}"] = directory
                directory = f"input_{key}"
        else: 
            self.directories["input"] = directory
            directory = "input"
        self.add(key, directory, "", file_name = filename, **kwargs)
    
    def mkdir(self, directory):
        if type(directory) is dict:
            self.directories.update(directory)
            directory = directory[next(iter(directory))]
        mkdir(self.directories[directory])
    
    def set(self, name):
        self.__dict__[name] = {}
    
    def set_file(self, set_id, directory, name, extension, **kwargs):
        new_file = file(self, directory, f"{name}{extension}", **kwargs)
        self.__dict__[set_id][name] = new_file
        self.__dict__[f"{set_id}_{name}"] = new_file
        return new_file
    
    def log(self, file, action):
        if not isinstance(self.file_log, file):
            self.file_log = file(self, "out", "file_log.tsv", quiet=True)
            with self.file_log.open("a") as output:
                output.write("file key\tfile path\taction\tframe\n")
        with self.file_log.open("a") as output:
            output.write("\t".join([file.key,file.path,action,self.frame])+"\n")
    


class directory:
    def __init__(self, container, path):
        self.container = container
        self.path = path
    
    def make(self):
        mkdir(self.path)
    
    def dir(self, key, name):
        self.container.__dict__[key] = os.path.join(self.path, name)
    
    def exists(self):
        return os.path.isdir(self.path)

    def remove(self):
        if self.exists():
            os.rmdir(self.path)

class file:
    def __init__(self, container, key, directory, name, quiet, **kwargs):
        self.container = container
        self.key = key
        self.directory = directory
        self.name = name
        self.quiet = quiet
        self.path = os.path.join(self.container.__dict__[directory].path, self.name)
        self.file_format = self.name[self.name.rindex(".")+1:]
        self.__dict__.update(kwargs)
        if not self.quiet: self.log("initialized")
    
    def add(self, **kwargs):
        self.__dict__.update(kwargs)
        if not self.quiet: self.log("metadata updated")
    
    def extend(self, key, extension, **kwargs):
        return self.container.extend(self, key, extension, **kwargs)
    
    def exists(self):
        return os.path.isfile(self.path)
    
    def open(self, options="r"):
        if not self.quiet: self.log(f"opened {options}")
        return open(self.path, options)
    
    def remove(self):
        if self.exists():
            os.remove(self.path)
            if not self.quiet: self.log("removed")

    def file_path(self, new_file):
        if new_file in self.container.__dict__:
            new_file = self.container.__dict__[new_file]
        if isinstance(new_file, file):
            new_file = new_file.path
        return new_file
    
    def rename(self, new_file):
        new_file = self.file_path(new_file)
        os.rename(self.path, new_file)
        self.file_format = f"renamed to {new_file}"
        if not self.quiet: self.log(self.file_format)
    
    def awk(self, new_file, column_array, ofs = "\t", index_correction = 1):
        '''Not called (at the moment) but meant as a replacement for lines 77-84 of TELR_te.py'''
        columns = "$" + ", $".join([str(column+index_correction) for column in column_array])
        new_file = self.file_path(new_file)
        command = f"awk -v OFS='{ofs}' '{{print {columns}}}' {self.path} > {new_file}"
        subprocess.run(command, shell=True)
        if not self.quiet: self.log(f"awk called, written to {new_file}")
    
    def log(self, action):
        self.container.log(self, action)



def rm_file(file):
    if os.path.exists(file):
        os.remove(file)


def mkdir(dir):
    if os.path.isdir(dir):
        print("Directory %s exists" % dir)
        return
    try:
        os.mkdir(dir)
    except OSError:
        print("Creation of the directory %s failed" % dir)
    else:
        print("Successfully created the directory %s " % dir)


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
        return "%d seconds" % (d.second)
    elif d.hour == 0 and d.minute != 0:
        return "%d minutes %d seconds" % (d.minute, d.second)
    else:
        return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)


def create_loci_set(vcf_parsed):
    #create set off all loci in vcf_parsed file as chr_start_end
    all_loci = set()
    with vcf_parsed.open() as input:
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
            output.write(
                ">" + record.id + "\n" + str(record.seq.reverse_complement()) + "\n"
            )


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

def name_contig(line):
    line = line.replace("\n", "").split("\t")
    return "_".join([line[0], line[1], line[2]])
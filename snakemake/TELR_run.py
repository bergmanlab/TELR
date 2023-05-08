import sys
import os
import random
import argparse
import subprocess
import logging
import json

def main():
    if len(sys.argv) > 2:
        config = get_args()
        config["tmp_dir"] = os.path.join(config["out"], "intermediate_files")
        config["verbose"] = False # add as option later
        if_verbose = verbose(config["verbose"])
        mkdir(if_verbose, config["tmp_dir"])
        mkdir(if_verbose, os.path.join(config["tmp_dir"], "input"))
        config["sample_name"] = os.path.splitext(os.path.basename(config["reads"]))[0]
        run_id = make_run_config(if_verbose, config)
    else: 
        run_id = sys.argv[1] #untested
    process_input_files(
        {"reads":config["reads"],
         "reference":config["reference"],
         "library":config["library"]},
         config["sample_name"]
         )
    run_workflow(config, run_id)

def make_run_config(if_verbose, config):

    # config["conda"] = os.path.join(os.path.dirname(os.path.abspath(__file__)),"envs/telr.yaml")

    snake_dir = os.path.join(config["out"], "snakemake") #make directory for snakemake
    mkdir(if_verbose, snake_dir)
    mkdir(if_verbose, f"{snake_dir}/config") #make directory for config file

    run_id = random.randint(1000000,9999999) # generate a random run ID
    run_config = f"{snake_dir}/config/config_{run_id}.json" # path to config file
    with open(run_config, "w") as conf:
        json.dump(config, conf, indent=4) #write config file as json
    
    return run_id

def run_workflow(config, run_id):
    #telr_dir = os.path.dirname(os.path.abspath(__file__)) #remember to implement later!
    command = [
        "snakemake", #"--use-conda",#"--conda-prefix",f"{telr_dir}/envs",
        "--configfile", f"{os.path.join(config['out'], 'snakemake')}/config/config_{run_id}.json",
        "--cores", str(config["thread"])#, "--quiet"
    ]
    subprocess.call(command)

def process_input_files(input_file_paths, sample_name):
    def file_extension_of(file):
        return os.path.splittext(input_file_paths[file])[1]
    accepted_formats_for = { #valid file extensions for each type of input
        "reads":[".fasta",".fastq",".fa",".fq",".bam"],
        "library":[".fasta",".fastq",".fa",".fq"],
        "reference":[".fasta",".fastq",".fa",".fq"]
        }
    for file in ["reads","reference","library"]:
        try: open(input_file_paths[file], "r")
        except Exception as e:
            print(e)
            logging.exception(f"Cannot open input file: {input_file_paths[file]}")
            sys.exit(1)
        if not file_extension_of(file) in accepted_formats_for[file]:
                print(f"Input {file} format not recognized, exiting...")
                logging.error("Input format not recognized")
                sys.exit(1)
        else: 
            symlink(input_file_paths[file], f"intermediate_files/input/{file}-{sample_name}{file_extension_of(file)}")

def symlink(input, output): #create a symbolic link at the output location referencing the input path.
    if os.path.islink(output):
        os.remove(output)
    try:
        os.symlink(input, output)
    except Exception:
        logging.exception(f"Create symbolic link for {input} failed")
        sys.exit(1)

def mkdir(if_verbose, dir):
    if os.path.isdir(dir):
        if_verbose.print(f"Directory {dir} exists") # verbose object print
        return
    try:
        os.mkdir(dir)
    except OSError:
        print(f"Creation of the directory {dir} failed")
    else:
        if_verbose.print(f"Successfully created the directory {dir}")

def get_args():
    parser = argparse.ArgumentParser(
        description="Program for detecting non-reference TEs in long read data"
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    # required
    required.add_argument(
        "-i",
        "--reads",
        type=str,
        help="reads in fasta/fastq format or read alignments in bam format",
        required=True,
    )
    required.add_argument(
        "-r",
        "--reference",
        type=str,
        help="reference genome in fasta format",
        required=True,
    )
    required.add_argument(
        "-l",
        "--library",
        type=str,
        help="TE consensus sequences in fasta format",
        required=True,
    )

    # optional
    optional.add_argument(
        "--aligner",
        type=str,
        help="choose method for read alignment, please provide 'ngmlr' or 'minimap2' (default = 'ngmlr')",
        required=False,
    )
    optional.add_argument(
        "--assembler",
        type=str,
        help="Choose the method to be used for local contig assembly step, please provide 'wtdbg2' or 'flye' (default = 'wtdbg2')",
        required=False,
    )
    optional.add_argument(
        "--polisher",
        type=str,
        help="Choose the method to be used for local contig polishing step, please provide 'wtdbg2' or 'flye' (default = 'wtdbg2')",
        required=False,
    )
    optional.add_argument(
        "-x",
        "--presets",
        type=str,
        help="parameter presets for different sequencing technologies, please provide 'pacbio' or 'ont' (default = 'pacbio')",
        required=False,
    )
    optional.add_argument(
        "-p",
        "--polish_iterations",
        type=int,
        help="iterations of contig polishing (default = 1)",
        required=False,
    )
    optional.add_argument(
        "-o",
        "--out",
        type=str,
        help="directory to output data (default = '.')",
        required=False,
    )
    optional.add_argument(
        "-t",
        "--thread",
        type=int,
        help="max cpu threads to use (default = '1')",
        required=False,
    )
    optional.add_argument(
        "-g",
        "--gap",
        type=int,
        help="max gap size for flanking sequence alignment (default = '20')",
        required=False,
    )
    optional.add_argument(
        "-v",
        "--overlap",
        type=int,
        help="max overlap size for flanking sequence alignment (default = '20')",
        required=False,
    )
    optional.add_argument(
        "--flank_len",
        type=int,
        help="flanking sequence length (default = '500')",
        required=False,
    )
    optional.add_argument(
        "--af_flank_interval",
        type=int,
        help="5' and 3'flanking sequence interval size used for allele frequency estimation (default = '100')",
        required=False,
    )
    optional.add_argument(
        "--af_flank_offset",
        type=int,
        help="5' and 3' flanking sequence offset size used for allele frequency estimation (default = '200')",
        required=False,
    )
    optional.add_argument(
        "--af_te_interval",
        type=int,
        help="5' and 3' te sequence interval size used for allele frequency estimation (default: '50')",
        required=False,
    )
    optional.add_argument(
        "--af_te_offset",
        type=int,
        help="5' and 3' te sequence offset size used for allele frequency estimation (default: '50')",
        required=False,
    )
    optional.add_argument(
        "--different_contig_name",
        action="store_true",
        help="If provided then TELR does not require the contig name to match before and after annotation liftover (default: require contig name to be the same before and after liftover)",
        required=False,
    )
    optional.add_argument(
        "--minimap2_family",
        action="store_true",
        help="If provided then minimap2 will be used to annotate TE families in the assembled contigs (default: use repeatmasker for contig TE annotation)",
        required=False,
    )
    optional.add_argument(
        "-k",
        "--keep_files",
        action="store_true",
        help="If provided then all intermediate files will be kept (default: remove intermediate files)",
        required=False,
    )
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # check if optional arguments are valid
    if args.aligner is None:
        args.aligner = "ngmlr"
    elif args.aligner not in ["ngmlr", "minimap2"]:
        print("Please provide a valid alignment method (ngmlr/minimap2), exiting...")
        sys.exit(1)

    if args.assembler is None:
        args.assembler = "wtdbg2"
    elif args.assembler not in ["wtdbg2", "flye"]:
        print("Please provide a valid assembly method (wtdbg2/flye), exiting...")
        sys.exit(1)

    if args.polisher is None:
        args.polisher = "wtdbg2"
    elif args.polisher not in ["wtdbg2", "flye"]:
        print("Please provide a valid polish method (wtdbg2/flye), exiting...")
        sys.exit(1)

    if args.presets is None:
        args.presets = "pacbio"
    elif args.presets not in ["pacbio", "ont"]:
        print("Please provide a valid preset option (pacbio/ont), exiting...")
        sys.exit(1)

    if args.polish_iterations is None:
        args.polish_iterations = 1
    elif args.polish_iterations < 1:
        print("Please provide a valid number of iterations for polishing, exiting...")

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    mkdir(verbose(False), args.out)

    if args.thread is None:
        args.thread = 1

    if args.flank_len is None:
        args.flank_len = 500

    if args.af_flank_interval is None:
        args.af_flank_interval = 100
    else:
        if args.af_flank_interval <= 0:
            print(
                "Please provide a valid flanking sequence interval size (positive integer) for allele frequency estimation, exiting..."
            )
            sys.exit(1)

    if args.af_flank_offset is None:
        args.af_flank_offset = 200
    else:
        if args.af_flank_offset < 0:
            print(
                "Please provide a valid flanking sequence offset size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.af_te_interval is None:
        args.af_te_interval = 50
    else:
        if args.af_te_interval <= 0:
            print(
                "Please provide a valid TE interval size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.af_te_offset is None:
        args.af_te_offset = 50
    else:
        if args.af_te_offset < 0:
            print(
                "Please provide a valid TE offset size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.gap is None:
        args.gap = 20

    if args.overlap is None:
        args.overlap = 20

    return vars(args)

class verbose:
    def __init__(self, verbose=False):
        self.verbose = verbose
    
    def print(self, string):
        if(self.verbose):
            print(string)

if __name__ == "__main__":
    main()
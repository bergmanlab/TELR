import os
import sys
import random
import json
import subprocess
telr_dir = f"{__file__.split('/evaluation')[0]}/telr"
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, memory_format, abs_path, mkdir
from interpret_config import config_from_file

def main():
    config = parse_args()
    config.update(get_file_paths())

    if not "resume" in config:
        #config.update(handle_file_paths(config))
        config.update(setup_run(config))
    else:
        run_id = config["resume"]
        tmp_dir = f"{config['out']}/stelr_eval_run_{run_id}"
        config_file = f"{tmp_dir}/config.json"
        with open(config_file,"r") as config_from_file:
            config = json.load(config_from_file)
        config["run_id"] = run_id
        config["resume"] = run_id
        config["tmp_dir"] = tmp_dir
        config["config_file"] = config_file

    snakefile = os.path.split(__file__)[0] + "/evaluation.smk"
    run_workflow(snakefile, config)


def get_file_paths():
    eval_src = abs_path(__file__)
    eval_src_dir = os.path.split(eval_src)[0]
    telr_dir = eval_src_dir.split("evaluation")[0] + "telr"
    env_dir = eval_src_dir.split("src")[0] + "envs"
    
    path_dict = {
        "telr":f"{telr_dir}/stelr.py",
        "telr_dir":telr_dir,
        "telr_conda":f"{env_dir}/stelr_sniffles1.yaml"
    }

    return path_dict

def parse_args():
    params = {("out",):abs_path()}
    args = {
        "-t":[("Resources","Threads"),int],
        "--threads":[("Resources","Threads"),int],
        "--mem":[("Resources","Estimated memory"),memory_format],
        "--memory":[("Resources","Estimated memory"),memory_format],
        "-o":[("out",),abs_path],
        "--resume":[("resume",),int]
    }
    unparsed = []
    for index in range(1,len(sys.argv)):
        prev = sys.argv[index-1]
        this = sys.argv[index]
        if prev in args:
            params[args[prev][0]] = args[prev][1](this)
        elif this in args:
            pass
        else:
            unparsed.append(this)

    config = config_from_file(unparsed[0])
    #print(params)
    for param in params:
        setdict(config,param,params[param])
    
    return config

def setup_run(config):
    run_id = random.randint(1000000,9999999) #generate a random run ID
    tmp_dir = f"{config['out']}/stelr_eval_run_{run_id}"
    mkdir(tmp_dir)
    
    source = config["read source"]
    #print(json.dumps(config, indent=4))
    if source == "reads from file":
        reads_name =  config["options for reads from file"]["Long read sequencing data"]
    elif source == "simulated reads":
        cov = config["options for simulated reads"]["coverage"]
        simulation_type = config["options for simulated reads"]["Simulation type"]
        if simulation_type == "proportion":
            proportion = config["options for simulated reads"]["Simulation type params"]["proportion"]
            reads_name =  f"{cov}x_{proportion[0]}-{proportion[1]}.fastq"
        else:
            genotype = config["options for simulated reads"]["Simulation type params"]["genotype"]
            reads_name =  f"{cov}x_{simulation_type}_{genotype}.fastq"

    config["output"] = [
        f"{reads_name}.telr.contig"
    ]

    config_path = f"{tmp_dir}/config.json" # path to config file
    with open(config_path, "w") as conf:
        json.dump(config, conf) #write config file as json    
    
    config["run_id"] = run_id
    config["tmp_dir"] = tmp_dir
    config["config_file"] = config_path

    setdict(config,("Resources","Estimated memory"),f'mem_mb={getdict(config,("Resources","Estimated memory"))}')

    return config

def run_workflow(snakefile, config):
    print(f"\nSTELR_evaluation run ID {config['run_id']}\n")
    command = [
        "snakemake",
        "-s", snakefile,
        "--configfile", config["config_file"],
        "--use-conda"
    ]
    optional_args = {("Resources","Threads"):"--cores",("Resources","Estimated memory"):"--resources"}
    for arg in optional_args:
        if getdict(config,arg):
            command += [optional_args[arg], str(getdict(config,arg))]
    try:
        if not "resume" in config:
            subprocess.call(command, cwd=config["tmp_dir"])
        else:
            subprocess.call(command + ["--unlock"], cwd=config["tmp_dir"])
            subprocess.call(command + ["--rerun-incomplete", "--rerun-triggers","mtime"], cwd=config["tmp_dir"])
        #for output_file in config["output"]:
        #    os.rename(f"{config['tmp_dir']}/{output_file}",f"{config['out']}/{output_file}")
        if not config["keep_files"]:
            rmtree(config['tmp_dir'])
    except Exception as e:
        print(e)
        print("STELR_evaluation failed!")
        sys.exit(1)
    print(f"STELR_evaluation run {config['run_id']} finished.")

if __name__ == '__main__':
    main()
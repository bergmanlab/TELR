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
    
    path_dict = {
        "telr":f"{telr_dir}/stelr.py",
        "telr_dir":telr_dir
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

def format_input(value, required_format, input_for_arg = "unknown"):
    #print([value, required_format, input_for_arg])
    if required_format == "file":
        try:
            path = abs_path(value)
            open(path, "r").close()
            return path
        except:
            if not os.path.isfile(value):
                print(f"Error: {input_for_arg} -- file '{value}' could not be found.", file=sys.stderr)
            else:
                print(f"Error: {input_for_arg} -- file '{value}' could not be read.", file=sys.stderr)
            sys.exit(1)
    elif required_format == "ploidy":
        try:
            return {
                "diploid":"diploid",
                "2":"diploid",
                "tetraploid":"tetraploid",
                "4":"tetraploid"
            }[value]
        except:
            print(f"Invalid ploidy type. Given ploidy type was '{value}'. Valid ploidy types are diploid or tetraploid.", file=sys.stderr)
            sys.exit(1)            
    elif required_format == "memory":
        if "mem_mb=" in value:
            return value
        else:
            try:
                return f"mem_mb={int(value)}"
            except:
                memory_units = {
                    ("tb","TB"):1000000,
                    ("gb","GB"):1000,
                    ("mb","MB"):1,
                    ("kb","KB"):0.001
                    }
                for unit in memory_units:
                    for abbreviation in unit:
                        if abbreviation in value:
                            try:
                                value = value.replace(abbreviation, "")
                                return f"mem_mb={int(value)*memory_units[unit]}"
                            except:
                                print("Error: input format for memory not recognized. Please enter the estimated memory in mb.")
                                sys.exit(1)
    else:
        format_function = {"integer":int,"string":str,"decimal value":float,"genbank accession number":str}[required_format]
        if type(value) is format_function:
            return value
        else:
            try:
                return format_function(value)
            except:
                print(f"Error: input {input_for_arg} formatted incorrectly. Expected format is {required_format}.", file=sys.stderr)
                sys.exit(1)

def get_allele_proportion(params):
    #determine the read depth required for each reference
    #TODO: it would be entirely possible to allow for >2 input genomes if that's a functionality that would be useful, eg for pooled sequence sims
    cov = params["cov"]
    if "proportion" in params:
        return params["proportion"], 1 - params["proportion"]
    else:
        genotype = params["genotype"]
        allele_proportion = {
            "heterozygous":[0.5, 0.5],
            "homozygous":[0, 1],
            "simplex":[0.75, 0.25],
            "duplex":[0.5, 0.5],
            "triplex":[0.25, 0.75],
            "quadruplex":[0, 1],
            "control":[1, 0]
        }
        if genotype in allele_proportion:
            return allele_proportion[genotype]
        else:
            genotype_type = ""
            valid_inputs = "heterozygous, homozygous, simplex, duplex, triplex, quadruplex"
            if "ploidy" in params:
                ploidy = params["ploidy"]
                valid_input_for_ploidy = {
                    "diploid":"heterozygous, homozygous",
                    "tetraploid":"simplex, duplex, triplex, quadruplex"
                }
                valid_inputs = valid_input_for_ploidy[ploidy]
                genotype_type = f" for {ploidy} genome"                        
            print(f"Invalid genotype{genotype_type}. Given genotype was '{genotype}'. Valid input types are {valid_inputs}, or control.", file=sys.stderr)
            sys.exit(1)

def get_simulation(params):
    cov = params["cov"]
    allele_proportion = get_allele_proportion(params)
    if "proportion" in params:
        return f"{cov}x_{allele_proportion[0]}-{allele_proportion[1]}"
    else:
        genotype = params["genotype"]
        if "ploidy" in params:
            ploidy = params['ploidy']
        else:
            ploidy = {
                "heterozygous":"diploid","homozygous":"diploid",
                "simplex":"tetraploid","duplex":"tetraploid","triplex":"tetraploid","quadruplex":"tetraploid",
                "control":"ref"
            }[genotype]
        return f"{cov}x_{ploidy}_{genotype}"

def get_simulation_config(params):
    config = {
        "simulation":get_simulation(params),
        "sub_simulations":{},
        "accessions":{}
    }
    cov = params["cov"]
    allele_proportion = get_allele_proportion(params)
    references = [param for param in params if "reference" in param]
    if "reference" in params: config["reference"] = format_input(params["reference"], "file", "reference")
    if 0 in allele_proportion and len(references) < len(allele_proportion):
        allele_proportion = [proportion for proportion in allele_proportion if proportion != 0]
    for reference in references:
        if not "reference" in config:
            config["reference"] = params[reference]["path"]
        if len(allele_proportion) > 0:
            proportion = allele_proportion.pop(0)
            if not proportion == 0:
                sub_cov = proportion*cov
                if sub_cov == int(sub_cov):
                    sub_cov = int(sub_cov)
                ref_name = params[reference]["name"]
                config["sub_simulations"][f"{sub_cov}x_{ref_name}"] = params[reference]["path"]
                if "accession" in params[reference]:
                    config["accessions"][ref_name] = params[reference]["accession"]
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
    print(f"STELR_evaluation run ID {config['run_id']}")
    command = [
        "snakemake",
        "-s", snakefile,
        "--configfile", config["config_file"]
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
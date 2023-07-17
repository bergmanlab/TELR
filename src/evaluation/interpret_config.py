import sys
telr_dir = "/home/hkg58926/github/TELR/src/telr"
sys.path.insert(0,telr_dir)
import json
import subprocess
from STELR_utility import getdict, setdict, check_exist, memory_format, abs_path
import traceback

def config_from_file(config_file):
    config = {}
    with open(config_file, "r") as input:
        try:
            config = json.load(input)
        except:
            pass 
    if not config:
        config = empty_default_config()
        info = [
            ["SIMULATION READS",[["[","]","selection",("read source",)]],21],
            ["] Mapping Reference:",[
                ["[","]",bool,("options for simulated reads","Reference Genomes","Mapping Reference given")],
                ["Mapping Reference:","\n",str,("options for simulated reads","Reference Genomes", "Mapping Reference","name")]
                ],2],
            ["] By path:",[
                ["[","]","selection",("options for simulated reads","Reference Genomes", "Mapping Reference","get by")],
                ["By path:","\n","file_path",("options for simulated reads","Reference Genomes", "Mapping Reference","file path")]
                ],0],
            ["] By accession:",[
                ["[","]","selection",("options for simulated reads","Reference Genomes", "Mapping Reference","get by")],
                ["By accession:","\n",str,("options for simulated reads","Reference Genomes", "Mapping Reference","accession")]
                ],0],
            ["] Community Reference:",[
                ["[","]",bool,("options for simulated reads","Reference Genomes","Community Reference given")],
                ["Community Reference:","\n",str,("options for simulated reads","Reference Genomes", "Community Reference","name")]
                ],2],
            ["] By path:",[
                ["[","]","selection",("options for simulated reads","Reference Genomes", "Community Reference","get by")],
                ["By path:","\n","file_path",("options for simulated reads","Reference Genomes", "Community Reference","file path")]
                ],0],
            ["] By accession:",[
                ["[","]","selection",("options for simulated reads","Reference Genomes", "Community Reference","get by")],
                ["By accession:","\n",str,("options for simulated reads","Reference Genomes", "Community Reference","accession")]
                ],0],
            ["] pbsim2",[["[","]","selection",("options for simulated reads", "Simulator")]],2],
            ["] By name:",[
                ["[","]","selection",("options for simulated reads","pbsim2","get by")],
                ["By name:","\n",str,("options for simulated reads","pbsim2","model name")],
                ],0],
            ["] From file:",[
                ["[","]","selection",("options for simulated reads","pbsim2","get by")],
                ["From file:","\n","file_path",("options for simulated reads","pbsim2","model file")],
                ],0],
            ["] Diploid:",[["[","]","genotype",("options for simulated reads","Simulation type")]],3],
            ["] Homozygous reference",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] Homozygous alternate",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] Heterozygous",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] Tetraploid:",[["[","]","genotype",("options for simulated reads","Simulation type")]],4],
            ["] Simplex",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] Duplex",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] Triplex",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] Quadruplex",[["[","]","genotype",("options for simulated reads","Simulation type params","genotype")]]],
            ["] By ratio:",[
                ["[","]","genotype",("options for simulated reads","Simulation type params","genotype")],
                ["By ratio:","\n",str,("options for simulated reads","Simulation type params","proportion")]
                ],0],
            ["Seed:",[["Seed:","\n",int,("options for simulated reads","seed")]]],
            ["Coverage:",[["Coverage:","\n",int,("options for simulated reads","coverage")]]],
            ["READS FROM FILE",[["[","]","selection",("read source",)]],4],
            ["Mapping Reference:",[
                ["Mapping Reference:","\n",str,("options for reads from file","Mapping Reference","name")]
                ],],
            ["] By path:",[
                ["[","]","selection",("options for reads from file","Mapping Reference","get by")],
                ["By path:","\n","file_path",("options for reads from file","Mapping Reference","file path")]
            ],0],
            ["] By accession:",[
                ["[","]","selection",("options for reads from file","Mapping Reference","get by")],
                ["By accession:","\n",str,("options for reads from file","Mapping Reference","accession")]
            ],0],
            ["Long read sequencing data:",[["Long read sequencing data:","\n","file_path",("options for reads from file","Long read sequencing data")]]],
            ["] TELR version 1.X",[["[","]","selection",("TELR command",)]]],
            ["] Specific conda environment:",[
                ["[","]","selection",("TELR parameters","Specific conda environment")],
                ["] Specific conda environment:","\n",str,("TELR parameters","Specific conda environment")]
                ],0],
            ["TE Library:",[["TE Library:","\n","file_path",("TELR parameters","TE library")]]],
            ["] Different mapping reference:",[
                ["[","]",bool,("TELR parameters","Different reference")],
                ["Different mapping reference:","\n",str,("TELR parameters","options for different reference", "name")]
                ]],
            ["] By path:",[
                ["[","]","selection",("TELR parameters","options for different reference","get by")],
                ["] By path:","\n","file_path",("TELR parameters","options for different reference","file path")]
                ],0],
            ["] By accession:",[
                ["[","]","selection",("TELR parameters","options for different reference","get by")],
                ["] By path:","\n",str,("TELR parameters","options for different reference","accession")]
                ],0],
            ["] NGMLR",[["[","]","selection",("TELR parameters","Aligner")]]],
            ["] MINIMAP2",[["[","]","selection",("TELR parameters","Aligner")]]],
            ["] wtdbg2",[["[","]","selection",("TELR parameters","Assembler")]]],
            ["] flye",[["[","]","selection",("TELR parameters","Assembler")]]],
            ["] wtdbg2",[["[","]","selection",("TELR parameters","Polisher")]]],
            ["] flye",[["[","]","selection",("TELR parameters","Polisher")]]],
            ["Polishing iterations:",[["iterations:","\n",int,("TELR parameters","Polishing iterations")]]],
            ["Max gap size:",[["Max gap size:","\n",int,("TELR parameters","Flanking sequence alignment parameters","Max gap size")]]],
            ["Max overlap size:",[["Max overlap size:","\n",int,("TELR parameters","Flanking sequence alignment parameters","Max overlap size")]]],
            ["Flanking sequence length:",[["Flanking sequence length:","\n",int,("TELR parameters","Flanking sequence alignment parameters","Flanking sequence length")]]],
            ["Flank sequence interval size:",[["Flank sequence interval size:","\n",int,("TELR parameters","Allele frequency estimation parameters","Flank sequence interval size")]]],
            ["Flank sequence offset:",[["Flank sequence offset:","\n",int,("TELR parameters","Allele frequency estimation parameters","Flank sequence offset")]]],
            ["TE sequence interval size:",[["TE sequence interval size:","\n",int,("TELR parameters","Allele frequency estimation parameters","TE sequence interval size")]]],
            ["TE sequence offset:",[["TE sequence offset:","\n",int,("TELR parameters","Allele frequency estimation parameters","TE sequence offset")]]],
            ["] Use minimap2 to annotate TE family names instead of RepeatMasker",[["[","]",bool,("TELR parameters","Additional options","Use minimap2 to annotate TE family names")]]],
            ["] Keep intermediate files",[["[","]",bool,("TELR parameters","Additional options","Keep intermediate files")]]],
            ["Mapping reference annotation:",[["Mapping reference annotation:","\n","file_path",("Evaluation requirements","Mapping reference annotation")]]],
            ["Community reference annotation:",[["Community reference annotation:","\n","file_path",("Evaluation requirements","Community reference annotation")]]],
            ["Regular recombination region:",[["Regular recombination region:","\n","file_path",("Evaluation requirements","Regular recombination region")]]],
            ["Estimated memory:",[["Estimated memory:","\n",memory_format,("Resources","Estimated memory")]]],
            ["Threads:",[["Threads:","\n",int,("Resources","Threads")]]]
        ]
        line_number = 0
        last_unstalled_line = 0
        with open(config_file, "r") as input:
            for line in input:
                line_number += 1
                if info[0][0] in line:
                    last_unstalled_line = line_number
                    info_lines, dependent_lines = info[0][1], info.pop(0)[-1]
                    for info_line in info_lines:
                        unselected = False
                        data = line.split(info_line[0])[-1].split(info_line[1])[0].strip()
                        input_format, location_in_config = info_line[2:4]
                        if type(input_format) is str:
                            if input_format in ["selection","genotype"]: unselected = not bool(data)
                            data = format_input(config, line, info_line, data)
                        else: data = input_format(data)
                        if data:
                            setdict(config, location_in_config, data)
                            if not input_format == "selection":
                                config.update(check_input(config, line_number, location_in_config))
                        if type(dependent_lines) is int and unselected:
                            info = info[dependent_lines:]
                            break
        if info:
            print(f"Error reading config: read config stalled after line {last_unstalled_line}")
            sys.exit(1)
    return config

def format_input(config, line, info_line, data):
    data_format, data_map = info_line[-2:]
    def selection(data):
        if not data: 
            return getdict(config, data_map)
        else:
            selections = {
                "path":"file path",
                "accession":"accession",
                "wtdbg2":"wtdbg2",
                "flye":"flye",
                "name":"model name",#keep an eye on this as more features are added
                "model file":"model file",
                "pbsim2":"pbsim2",
                "reads from file":"reads from file",
                "simulation reads":"simulated reads",
                "from file":"model file",
                "ngmlr":"ngmlr",
                "minimap2":"minimap2",
                "telr version 1.x":"telr",
                "specific conda environment":True
            }
            for selection in selections:
                if selection in line.lower():
                    return selections[selection]
    def genotype(data):
        if not data: 
            return None
        else:
            genotypes = {
                "diploid":"diploid",
                "tetraploid":"tetraploid",
                "homozygous ref":"homozygous ref",
                "homozygous alt":"homozygous alt",
                "heterozygous":"heterozygous",
                "simplex":"simplex",
                "duplex":"duplex",
                "triplex":"triplex",
                "quadruplex":"quadruplex",
                "by ratio":"proportion"
            }
            for genotype in genotypes:
                if genotype in line.lower():
                    return genotypes[genotype]
    return {
        "selection":selection,
        "file_path":abs_path,
        "genotype":genotype
    }[data_format](data)

def check_input(config, line_number, data_map):
    def is_optional():
        parent_enabled = True
        for item in data_map:
            if "options for" in item:
                if "different reference" in item:
                    parent_enabled = getdict(config,("TELR parameters","Different reference"))
                else:
                    parent_enabled = config["read source"] in item
                break
        parent = data_map[:-1]
        return "get by" in getdict(config, parent) and parent_enabled
    def need_to_check():
        need = True
        if is_optional():
            get_by = list(data_map[:-1]) + ["get by"]
            if getdict(config, get_by) != data_map[-1]: 
                need = False
        return need
    def check_file(data):
        if need_to_check():
            if not check_exist(data):
                description = {True:"selected",False:"required"}[is_optional()]
                print(f"Input error in config line {line_number}: {description} file path could not be found.")
                sys.exit(1)
    def check_accession(data):
        if need_to_check():
            if not data[:4] == "GCA":
                data = data.upper()
                if data[:4] == "GCA":
                    setdict(config, data_map, data)
                else:
                    print(f"Input error in config line {line_number}: improper format for genbank accession number. Accession number should begin with \"GCA\".")
                    sys.exit(1)
            command = ["datasets", "summary", "genome", "accession", data]
            try:                    
                output = subprocess.run(command, capture_output=True, text=True).stdout
                genbank_info = json.loads(output)
                if genbank_info['total_count'] != 1:
                    sys.exit(1)
            except:
                print(f"Input error in config line {line_number}: could not find genome for accession #{value}", file=sys.stderr)
                sys.exit(1)
    def proportion(data):
        if getdict(config, ["options for simulated reads","Simulation type params","genotype"]) == "proportion":
            try:
                data = data.split()
                data = [chunk for chunk in data if ":" in chunk][0].split(":")
                data = [float(chunk) for chunk in data]
                if data[0] + data[1] > 0.99 and data[0] + data[1] < 1.01:
                    setdict(config, data_map, data)
                else: raise Exception(f"Proportion total is {data[0] + data[1]}")
            except Exception as e:
                print(e)
                print(f"Input error in config line {line_number}: incorrect format for proportion. Please give input in format reference:alternate, and make sure the proportions add up to 1.")
                sys.exit(1)
    def option_from_list(options):
        if not data in options:
            raise_error = True
            for option in options:
                try:
                    if data.lower() == option.lower():
                        setdict(config, data_map, option)
                        raise_error = False
                except: pass
            if raise_error:
                if type(options) in [list, tuple]:
                    if type(options[0]) in [int, float] and type(options[-1]) in [int, float]:
                        if data >= options[0] and data <= options[-1]:
                            raise_error = False
            if raise_error:
                print(f"Input error in config line {line_number}: input {data} does not match criteria given by {options}.")
                sys.exit(1)
    def genotype(data):
        try:
            ploidy = getdict(config, ("options for simulated reads","Simulation type"))
            genotypes = {
                "diploid":["homozygous ref","homozygous alt","heterozygous"],
                "tetraploid":["simplex","duplex","triplex","quadruplex"]
            }[ploidy]
            if not data in genotypes:
                if data.lower() in genotypes:
                    data = data.lower()
                    setdict(config, data_map, data)
                else: sys.exit(1)
            proportion = {"homozygous ref":[1,0],"homozygous alt":[0,1],"heterozygous":[0.5,0.5],
                "simplex":[0.75,0.25],"duplex":[0.5,0.5],"triplex":[0.25,0.75],"quadruplex":[0,1]}[data]
            setdict(config, ["options for simulated reads","Simulation type params","proportion"], proportion)
        except:
            print(f"Input error in config line {line_number}: {data} not an accepted genotype for given ploidy.")
            sys.exit(1)
    check_dict = {
        "read souce":"simulated reads or file",#implement this
        "options for simulated reads":{
            "Reference Genomes":{
                "Reference Genome":{
                    "file path":check_file,
                    "accession":check_accession
                },
                "Alternate Genome":{
                    "file path":check_file,
                    "accession":check_accession
                }
            },
            "Simulator":{"pbsim2"},
            "pbsim2":{
                "model name":{"P4C2", "P5C3", "P6C4", "R103", "R94", "R95"},
                "model file":check_file
            },
            "Simulation type":{"diploid","tetraploid","proportion"},
            "Simulation type params":{
                "genotype":genotype,
                "proportion":proportion
            },
            "seed":[0,float('inf')],
            "coverage":[0,float('inf')]
        },
        "options for reads from file":{
            "Reference Genome":{
                "file path":check_file,
                "accession":check_accession
            },
            "Long read sequencing data":check_file
        },
        "TELR parameters":{
            "Specific conda environment":check_file,
            "TE library":check_file,
            "options for different reference":{
                "file path":check_file,
                "accession":check_accession
            },
            "Aligner":["ngmlr","minimap2"],
            "Assembler":["wtdbg2","flye"],
            "Polisher":["wtdbg2","flye"],
            "Polishing iterations":[1,float('inf')],
            "Flanking sequence alignment parameters":{
                "Max gap size":[0,float('inf')],
                "Max overlap size":[0,float('inf')],
                "Flanking sequence length":[0,float('inf')]
            },
            "Allele frequency estimation parameters":{
                "Flank sequence interval size":[0,float('inf')],
                "Flank sequence offset":[0,float('inf')],
                "TE sequence interval size":[0,float('inf')],
                "TE sequence offset":[0,float('inf')]
            }
        },
        "Resources":{
            "Threads":[1,float('inf')]
        }
    }
    try:
        method = getdict(check_dict, data_map)
    except:
        return config
    data = getdict(config, data_map)
    try:
        if callable(method):
            method(data)
        else:
            option_from_list(method)
    except Exception as e:
        traceback.print_exc()
        print(f"Input error in config line {line_number}: unable to verify valid input format for '{data}'.")
        sys.exit(1)
    return config

            
def empty_default_config():
    return {
        "read source":"simulated reads or reads from file",
        "options for simulated reads":{
            "Reference Genomes":{
                "Mapping Reference given":True,
                "Mapping Reference":{
                    "name":"name of genome",
                    "get by":"file path or accession",
                    "file path":"path to file",
                    "accession":"genbank accession number"
                },
                "Community Reference given":True,
                "Community Reference":{
                    "name":"name of genome",
                    "get by":"file path or accession",
                    "file path":"path to file",
                    "accession":"genbank accession number"
                }
            },
            "Simulator":"pbsim2",
            "pbsim2":{
                "model name":"P6C4",
                "model file":"path to file",
                "get by":"model file"
            },
            "Simulation type":"diploid, tetraploid, or proportion",
            "Simulation type params":{
                "genotype":"homozygous ref, homozygous alt, heterozygous, simplex, duplex, triplex, quadruplex",
                "proportion":[0.5,0.5]
            },
            "seed":100,
            "coverage":50
        },
        "options for reads from file":{
            "Mapping Reference":{
                "name":"name of genome",
                "get by":"file path or accession",
                "file path":"path to file",
                "accession":"genbank accession number"
            },
            "Long read sequencing data":"path to file"
        },
        "TELR command":f"python3 {abs_path(__file__).split('evaluation')[0]}telr/stelr.py",
        "TELR parameters":{
            "Specific conda environment":False,
            "TE library":"path to file",
            "Different mapping reference":False,
            "options for different reference":{
                "name":"name of genome",
                "get by":"file path or accession",
                "file path":"path to file",
                "accession":"genbank accession number"
            },
            "Aligner":"NGMLR",
            "Assembler":"wtdbg2",
            "Polisher":"flye",
            "Polishing iterations":1,
            "Flanking sequence alignment parameters":{
                "Max gap size":20,
                "Max overlap size":20,
                "Flanking sequence length":500
            },
            "Allele frequency estimation parameters":{
                "Flank sequence interval size":100,
                "Flank sequence offset":200,
                "TE sequence interval size":50,
                "TE sequence offset":50
            },
            "Additional options":{
                "Use minimap2 to annotate TE family names":False,
                "Keep intermediate files":False
            }
        },
        "Evaluation requirements":{
            "Mapping reference annotatation":"path to file",
            "Community reference annotation":"path to file",
            "Regular recombination region":"path to file"
        },
        "Resources":{
            "Estimated memory":"100gb",
            "Threads":1
        }
    }


if __name__ == "__main__":
    config = config_from_file(sys.argv[1])
    print(json.dumps(config, indent=4))
import sys
import json

config = {
    "read source":"simulated reads or file",
    "options for simulated reads":{
        "Reference Genomes":{
            "Reference genome given":True,
            "Reference Genome":{
                "name":"name of genome",
                "get by":"file path or accession",
                "file path":"path to file",
                "accession":"genbank accession number"
            },
            "Alternate genome given":True,
            "Alternate Genome":{
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
        "Reference Genome":{
            "name":"name of genome",
            "get by":"file path or accession",
            "file path":"path to file",
            "accession":"genbank accession number"
        },
        "Long read sequencing data":"path to file"
    },
    "TELR parameters":{
        "TE library":"path to file",
        "Different reference":False,
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
    "Resources":{
        "Estimated memory":"100gb",
        "Threads":1
    }
}

with open(sys.argv[1], "w") as output:
    json.dump(config, output)
import os
telr_dir = f"{__file__.split("evaluation")[0]}/telr"
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, abs_path


rule all:
    input:
        config["output"]


def if_get_genome(wildcards):
    genome = os.path.basename(wildcards.reference)
    genomes_in_config = [
        ("options for simulated reads","Reference Genomes","Reference Genome"),
        ("options for simulated reads","Reference Genomes","Alternate Genome"),
        ("options for reads from file","Reference Genome"),
        ("TELR parameters","options for different reference"),
    ]
    for genome_from in genomes_in_config:
        if getdict(config,genome_from+("name",)) == genome:
            get_by = getdict(config,genome_from+("get by",))
            if get_by == "accession":
                return getdict(config,genome_from+("accession",))
            else: return ""#left off here
    if genome in config["accessions"]:
        return config["accessions"][genome]
    else:
        return ""
rule download_genome:
    params:
        accession = if_get_genome
    output:
        "{reference}.fna"
    shell:
        """
        mkdir ncbi_temp_{params.accession}
        cd ncbi_temp_{params.accession}
        datasets download genome accession {params.accession}
        unzip ncbi_dataset.zip
        cd ..
        mv ncbi_temp_{params.accession}/ncbi_dataset/data/{params.accession}/*_genomic.fna {output}
        rm -r ncbi_temp_{params.accession}
        """

def if_get_model(wildcards):
    if "model_name" in config:
        return config["model_name"]
    else:
        return ""
rule download_model:
    params:
        model = if_get_model
    output:
        "{model}.model"
    shell:
        """
        wget -L https://raw.githubusercontent.com/yukiteruono/pbsim2/master/data/{wildcards.model}.model -O "{output}"
        """


rule simulate_reads:
    input:
        reference = lambda wildcards: config["sub_simulations"][f"{wildcards.cov}x_{wildcards.ref_name}"],
        model = config["model"]
    output:
        out = "{cov}x_{ref_name}_pbsim.fq",
        log = "{cov}x_{ref_name}_pbsim.report"
    params:
        prefix = "'{cov}x_{ref_name}_pbsim'"
    resources: 
        mem_mb = 100000
    shell:
        """
        mkdir {params.prefix}
        cd {params.prefix}
        pbsim --depth {wildcards.cov} --prefix {params.prefix} --id-prefix '{wildcards.ref_name}' --hmm_model '{input.model}' '{input.reference}' 2> '../{output.log}'
        cd ..
        cat {params.prefix}/*.fastq > {output.out}
        rm -r {params.prefix}
        """

rule combine_simulation:
    input:
        lambda wildcards: [f"{sub_sim}_pbsim.fq" for sub_sim in config["sub_simulations"]]
    output:
        "{cov}x_{proportion_genotype}.fastq"
    shell:
        """
        cat {input} > {output}
        """

rule run_telr:
    input:
        reference = config["reference"],
        reads = lambda wildcards: f"{config['simulation']}.fastq",
        library = config["te_library"],
    output:
        "{sample_name}.telr.{output}"
    params:
        polish_iterations = config["polish_iterations"],
        assembler = config["assembler"],
        polisher = config["polisher"]
    threads: config["thread"]
    shell:
        """
        python3 {config[telr]} -i {input.reads} -r {input.reference} -l {input.library} -t {threads} -k -p {params.polish_iterations} --assembler {params.assembler} --polisher {params.polisher}
        """


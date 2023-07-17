import os
import subprocess
import json
from STELR_utility import get_contig_length

rule all:
    input: 
        config["output"]

def input_reads_if_in_bam_format(wildcards):
    #   this rule returns the bam format input if the input is given in bam format, or an empty list if not.
    #returning an empty list essentially gives snakemake permission to accept the input file in fasta 
    #format; otherwise it tries to check the next file back in the workflow, and causes a cyclic dependency
    #error if the input and output for the rule bam_input are the same.
    if ".bam" in config["reads"]: return config["reads"]
    else: return []
rule bam_input: #if input is given in bam format, convert it to fasta format.
    #todo -- check if this bam format is ACTUALLY aligned.
    input:
        input_reads_if_in_bam_format
    output:
        config["fasta_reads"]
    shell:
        "python3 {config[STELR_alignment]} bam2fasta '{input}' '{output}'"

"""
1st stage: identify TE insertion candidate loci
"""

'''1st stage
"Read alignment (NGMLR)"
^(or minimap2)
'''
#only run if reads are supplied in fasta format
def find_alignment(wildcards):
    bam_input = f"input/reads.bam"
    if(os.path.isfile(bam_input)): return bam_input
    else: 
        return f"{config['aligner']}_alignment.sam"
rule sort_index_bam:
    input:
        find_alignment#gives name of input file (this depends on whether the user supplied input was in bam format or it was aligned later)
    output:
        "reads_sort.bam"
    threads: config["thread"]
        #samtools
    shell:
        "python3 {config[STELR_alignment]} sort_index_bam '{input}' '{output}' '{threads}'"
    
rule align_with_ngmlr:#TODO: add timers to these alignments
    input:
        reads = config["fasta_reads"],
        reference = config["reference"]
    output:
        "ngmlr_alignment.sam"
    params:
        presets = config["presets"],
        label = lambda wildcards: {"ont":"ont","pacbio":"pb"}[config["presets"]]
    threads: config["thread"]
    shell:
        """
        ngmlr -r {input.reference} -q {input.reads} -x {params.presets} -t {threads} --rg-id reads --rg-sm reads --rg-lb {params.label} --no-progress | python3 {config[fix_ngmlr]} > {output}
        """

rule align_with_minimap2:
    input:
        reads = config["fasta_reads"],
        reference = config["reference"]
    output:
        "minimap2_alignment.sam"
    params:
        presets = lambda wildcards: {"ont":"map-ont","pacbio":"map-pb"}[config["presets"]]
    threads: config["thread"]
    shell:
        """
        minimap2 --cs --MD -Y -L -ax {params.presets} {input.reference} {input.reads} > {output}
        """


'''1st stage
SV calling (Sniffles)
'''
'''
rule detect_sv:
    input:
        bam = "reads_sort.bam",
        reference = config["reference"]
    output:
        "sv-reads_{sv_detector}.vcf"
    params:
        sample_name = "reads",
    threads: config["thread"]
    shell:
        "python3 {config[STELR_sv]} detect_sv '{input.bam}' '{input.reference}' '{output}' '{params.sample_name}' '{threads}'"
'''

rule sv_detection_sniffles1:
    input:
        "reads_sort.bam"
    output:
        "sv-reads_Sniffles1.vcf"
    threads: config["thread"]
    shell:
        "sniffles -n -1 --threads {threads} -m {input} -v {output}"

rule sv_detection_sniffles2:
    input:
        "reads_sort.bam"
    output:
        "sv-reads_Sniffles2.vcf"
    threads: config["thread"]
    shell:
        """
        sniffles --output-rnames -t 10 -i {input} -v {output}
        """
def sv_detector(wildcards):
    if "sv_detector" in config:
        return config["sv_detector"]
    else:
        return "Sniffles1"
def parse_vcf_input(wildcards):
    return f"sv-reads_{sv_detector(wildcards)}.vcf"
rule parse_vcf:
    input:
        sv_detector
    output:
        "reads.vcf_parsed.tsv.tmp"
    params:
        lambda wildcards: {
            "Sniffles1":'%CHROM\\t%POS\\t%END\\t%SVLEN\\t%RE\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n',
            "Sniffles2":'%CHROM\\t%POS\\t%END\\t%SVLEN\\t%AF\\t%ID\\t%ALT\\t%RNAMES\\t%FILTER\\t[ %GT]\\t[ %DR]\\t[ %DV]\n'
        }[sv_detector(wildcards)]
    shell:
        'bcftools query -i \'SVTYPE="INS" & ALT!="<INS>"\' -f "{params}" "{input}" > "{output}"'

rule swap_vcf_coordinate:
    input:
        "reads.vcf_parsed.tsv.tmp"
    output:
        "reads.vcf_parsed.tsv.swap"
    params:
        lambda wildcards: sv_detector(wildcards)
    shell:
        "python3 {config[STELR_sv]} swap_coordinate '{input}' '{output}' {params}"

rule rm_vcf_redundancy:
    input:
        "reads.vcf_parsed.tsv.swap"
    output:
        "reads.vcf_parsed.tsv"
    shell:
        "python3 {config[STELR_sv]} rm_vcf_redundancy '{input}' '{output}'"

rule write_ins_seqs:
    input:
        "reads.vcf_parsed.tsv"
    output:
        "reads.vcf_ins.fasta"
    shell:
        "python3 {config[STELR_sv]} write_ins_seqs '{input}' '{output}'"

'''1st stage
Filter for TE insertion candidate (RepeatMasker)
'''

rule sv_repeatmask:
    input:
        ins_seqs = "reads.vcf_ins.fasta",
        library = config["library"]
    output:
        "vcf_ins_repeatmask/{ins_seqs}.out.gff"
    params:
        repeatmasker_dir = "vcf_ins_repeatmask",
        thread = config["thread"]
    shell:
        "python3 {config[STELR_sv]} repeatmask '{params.repeatmasker_dir}' '{input.ins_seqs}' '{input.library}' '{params.thread}'"

rule sv_RM_sort:
    input:
        "vcf_ins_repeatmask/{ins_seqs}.out.gff"
    output:
        "vcf_ins_repeatmask/{ins_seqs}.out.sort.gff"
        #bedtools
    shell:
        "bedtools sort -i '{input}' > '{output}'"

rule sv_RM_merge:
    input:
        "vcf_ins_repeatmask/{ins_seqs}.out.sort.gff"
    output:
        "vcf_ins_repeatmask/{ins_seqs}.out.merge.bed"
        #bedtools
    shell:
        "bedtools merge -i '{input}' > '{output}'"

rule sv_TE_extract:
    input:
        parsed_vcf = "reads.vcf_parsed.tsv",
        ins_seqs = "reads.vcf_ins.fasta",
        ins_rm_merge = "vcf_ins_repeatmask/reads.vcf_ins.fasta.out.merge.bed"
    output:
        ins_filtered = "reads.vcf.filtered.tmp.tsv",
        loci_eval = "reads.loci_eval.tsv"
    shell:
        "python3 {config[STELR_sv]} te_extract '{input.parsed_vcf}' '{input.ins_seqs}' '{input.ins_rm_merge}' '{output.ins_filtered}' '{output.loci_eval}'"

rule seq_merge:
    input:
        "reads.vcf.filtered.tmp.tsv"
    output: 
        "reads.vcf.merged.tmp.tsv"
    params:
        window = 20
        #bedtools
    shell:
        'bedtools merge -o collapse -c 2,3,4,5,6,7,8,9,10,11,12,13,14 -delim ";" -d "{params.window}" -i "{input}" > "{output}"'

checkpoint merge_parsed_vcf:##### can we thread back to here? (probably not easily)
    input:
        "reads.vcf.merged.tmp.tsv"
    output:
        "reads.vcf_filtered.tsv"
    shell:
        "python3 {config[STELR_sv]} merge_vcf '{input}' '{output}'"

"""
2nd stage: assembly and polish local TE contig
"""

'''2nd stage
Local contig assembly and polishing (wtdbg2/flye + minimap2)
'''

rule initialize_contig_dir:
    input:
        "reads.vcf_filtered.tsv"
    output:
        "contigs/{contig}/00_vcf_parsed.tsv"
    shell:
        "python3 {config[STELR_assembly]} make_contig_dir '{input}' '{wildcards.contig}' '{output}'"

rule get_read_ids: # get a list of all the read IDs from the parsed vcf file
    input:
        "contigs/{contig}/00_vcf_parsed.tsv"
    output:
        "contigs/{contig}/00_reads.id"
    shell:
        "python3 {config[STELR_assembly]} write_read_IDs '{input}' '{wildcards.contig}' '{output}'"

rule unique_IDlist: # get a list of unique IDs from the readlist
    input:
        "contigs/{contig}/00_{readlist}.id"
    output:
        "contigs/{contig}/00_{readlist}.id.unique"
    shell:
        "cat '{input}' | sort | uniq > '{output}'"

rule filter_readlist: # use seqtk to get the fasta reads from the input reads file
    input:
        "contigs/{contig}/00_{readlist}.id.unique"
    output:
        "contigs/{contig}/00_{readlist}.fa"
    shell:
        "seqtk subseq '{config[fasta_reads]}' '{input}' | seqtk seq -a > '{output}'"

rule run_assembly:
    input:
        "contigs/{contig}/00_reads.fa"
    output:
        "contigs/{contig}/01_initial_assembly.fa"
    threads: 1
    shell:
        """
        python3 {config[STELR_assembly]} run_{config[assembler]}_assembly '{input}' '{wildcards.contig}' '{threads}' '{config[presets]}' '{output}' || true
        touch '{output}'
        """

rule run_polishing:
    input:
        reads = "contigs/{contig}/00_reads.fa",
        initial_assembly = "contigs/{contig}/01_initial_assembly.fa"
    output:
        "contigs/{contig}/02_polished_assembly.fa"
    threads: 1
    shell:
        """
        python3 {config[STELR_assembly]} run_{config[polisher]}_polishing '{input.initial_assembly}' '{output}' '{input.reads}' '{wildcards.contig}' '{threads}' '{config[polish_iterations]}' '{config[presets]}' || true
        touch '{output}'
        """

rule get_parsed_contigs:
    input:
        "contigs/{contig}/02_polished_assembly.fa"
    output:
        "contigs/{contig}/03_contig1.fa"
    shell:
        """
        python3 {config[STELR_assembly]} parse_assembled_contig '{input}' '{wildcards.contig}' '{output}' || true
        touch '{output}'
        """

def merged_contigs_input(wildcards):
    return [f"contigs/{contig}/03_contig1.fa" for contig in all_contigs(**wildcards)]
def all_contigs(wildcards):
    vcf_parsed_file = checkpoints.merge_parsed_vcf.get(**wildcards).output[0]
    with open(vcf_parsed_file) as vcf_parsed:
        contigs = []
        for line in vcf_parsed:
            contigs.append("_".join(line.split("\t")[:3]))
    return contigs
rule merged_contigs:
    input:
        merged_contigs_input
    output:
        "reads.contigs.fa"
    shell:
        """
        cat '{input}' > '{output}'
        """

"""
3rd stage: annotate TE and predict location in reference genome
"""

'''3rd stage
Contig TE annotation (minimap2 + RepeatMasker)
'''

rule get_vcf_seq:
    input:
        contig = "contigs/{contig}/00_reads.fa",
        vcf_parsed = "contigs/{contig}/00_vcf_parsed.tsv"
    output:
        "contigs/{contig}/04_vcf_seq.fa"
    shell:
        """
        python3 {config[STELR_te]} get_vcf_seq '{wildcards.contig}' '{input.vcf_parsed}' '{output}' || true
        touch '{output}'
        """

rule map_contig:
    input:
        subject = "contigs/{contig}/03_contig1.fa",
        query = "contigs/{contig}/04_vcf_seq.fa"
    output:
        "contigs/{contig}/05_vcf_mm2.paf"
    params:
        presets = lambda wildcards: {"pacbio":"map-pb","ont":"map-ont"}[config["presets"]]
    threads: 1
    shell:
        """
        minimap2 -cx '{params.presets}' --secondary=no -v 0 -t '{threads}' '{input.subject}' '{input.query}' > '{output}'
        """

rule te_contig_map:
    input:
        minimap_initial = "contigs/{contig}/05_vcf_mm2.paf",
        subject = "contigs/{contig}/03_contig1.fa",
        library = config["library"]
    output:
        "contigs/{contig}/06_te_mm2.paf"
    params:
        presets = lambda wildcards: {"pacbio":"map-pb","ont":"map-ont"}[config["presets"]]
    threads: 1
    shell:
        """
        if [ -s '{input.minimap_initial}' ]; then
            minimap2 -cx '{params.presets}' '{input.subject}' '{input.library}' -v 0 -t '{threads}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule minimap2bed:
    input:
        "{minimap_output}.paf"
    output:
        "{minimap_output}.bed"
    shell:
        "python3 {config[STELR_utility]} minimap2bed '{input}' '{output}'"

rule vcf_alignment_filter_intersect:
    input:
        vcf_seq_mm2 = "contigs/{contig}/05_vcf_mm2.bed",
        te_mm2 = "contigs/{contig}/06_te_mm2.bed"
    output:
        "contigs/{contig}/07_te2contig_filter.tsv"
    shell:
        """
        if [ -s '{input.te_mm2}' ]; then
            bedtools intersect -a '{input.te_mm2}' -b '{input.vcf_seq_mm2}' -wao > '{output}'
        else
            touch '{output}'
        fi
        """

rule vcf_alignment_filter:
    input:
        "contigs/{contig}/07_te2contig_filter.tsv"
    output:
        "contigs/{contig}/08_te2contig_filtered.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_te]} vcf_alignment_filter '{input}' '{output}'
        fi
        touch '{output}'
        """

rule te_annotation_sort:
    input:
        "contigs/{contig}/08_te2contig_filtered.bed"
    output:
        "contigs/{contig}/08_te2contig_sorted.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            bedtools sort -i {input} > {output}
        fi
        touch '{output}'
        """

rule te_annotation_merge:
    input:
        "contigs/{contig}/08_te2contig_sorted.bed"
    output:
        "contigs/{contig}/09_te2contig_merged.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            bedtools merge -d 10000 -c 4,6 -o distinct,distinct -delim "|" -i '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

checkpoint annotate_contig:
    input:
        "contigs/{contig}/09_te2contig_merged.bed"
    output:
        "contigs/{contig}/tes/annotation.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_te]} annotate_contig '{input}' '{output}'
        else
            touch '{output}'
        fi
        """

rule make_te_dirs:
    input:
        "contigs/{contig}/tes/annotation.bed"
    output:
        "contigs/{contig}/tes/{te}/00_annotation.bed"
    shell:
        "python3 {config[STELR_te]} make_te_dir '{input}' '{output}'"


## use RM to annotate config

rule rm_te_fasta:
    input:
        bed_file = "contigs/{contig}/tes/annotation.bed",
        sequence = "contigs/{contig}/03_contig1.fa"
    output:
        "contigs/{contig}/rm_01_annotated_tes.fasta"
    shell:
        """
        if [ -s '{input.bed_file}' ]; then
            bedtools getfasta -fi '{input.sequence}' -bed '{input.bed_file}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule rm_annotate:
    input:
        te_fasta = "contigs/{contig}/rm_01_annotated_tes.fasta",
        te_library = config["library"]
    output:
        "contigs/{contig}/rm_01_annotated_tes.fasta.out.gff"
    params:
        rm_dir = lambda wildcards: f"contigs/{wildcards.contig}'"
    threads: 1
    shell:
        """
        if [ -s '{input.te_fasta}' ]; then
            RepeatMasker -dir '{params.rm_dir}' -gff -s -nolow -no_is -xsmall -e ncbi -lib '{input.te_library}' -pa '{threads}' '{input.te_fasta}'
        fi
        touch '{output}'
        """

rule rm_annotation_sort:
    input:
        "contigs/{contig}/rm_01_annotated_tes.fasta.out.gff"
    output:
        "contigs/{contig}/rm_02_annotated_tes.out.sort.gff"
    shell:
        """
        if [ -s '{input.bed_file}' ]; then
            bedtools sort -i '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule rm_annotation_parse_merge:
    input:
        "contigs/{contig}/rm_02_annotated_tes.out.sort.gff"
    output:
        "contigs/{contig}/rm_03_annotated_tes_parsed.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_te]} rm_parse_merge '{input}' '{output}'
        else
            touch '{output}'
        fi
        """

rule rm_annotation_bedtools_merge:
    input:
        "contigs/{contig}/rm_03_annotated_tes_parsed.bed"
    output:
        "contigs/{contig}/rm_04_annotated_tes_merged.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            bedtools merge -c 4,6 -o distinct -delim "|" -i '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule rm_reannotate:
    input:
        repeat_masker_out = "contigs/{contig}/rm_04_annotated_tes_merged.bed",
        original_bed = "contigs/{contig}/tes/annotation.bed"
    output:
        "contigs/{contig}/rm_05_annotation.bed"
    shell:
        """
        if [ -s '{input.repeat_masker_out}' ]; then
            python3 {config[STELR_te]} rm_reannotate '{input.repeat_masker_out}' '{input.original_bed}' '{output}'
        else
            touch '{output}'
        fi
        """

# repeatmask reference genome using custom TE library
#   Not sure which step in the workflow this actually belongs to
rule ref_repeatmask:
    input:
        ref = config["reference"],
        lib = config["library"]
    output:
        "ref_repeatmask/{reference}.masked",
        "ref_repeatmask/{reference}.out.gff",
        "ref_repeatmask/{reference}.out"
    params:
        ref_rm_dir = "ref_repeatmask"
    threads: config["thread"]
    shell:
        """
        if [ ! -d '{params.ref_rm_dir}' ]; then mkdir '{params.ref_rm_dir}'
        fi
        RepeatMasker -dir '{params.ref_rm_dir}' -gff -s -nolow -no_is -e ncbi -lib '{input.lib}' -pa '{threads}' '{input.ref}'
        touch {output}
        """

rule ref_rm_process:
    input:
        gff = "ref_repeatmask/{reference}.out.gff",
        out = "ref_repeatmask/{reference}.out"
    output:
        "ref_repeatmask/{reference}.out.gff3"
    shell:
        "python3 {config[STELR_te]} parse_rm_out '{input.gff}' '{input.out}' '{output}'"
        #left off here, config[STELR_te]} repeatmask()

rule ref_te_bed:
    input:
        "ref_repeatmask/{reference}.out.gff3"
    output:
        "ref_repeatmask/{reference}.te.bed.unsorted"
    shell:
        """
        python3 {config[STELR_te]} gff3tobed '{input}' '{output}'
        touch '{output}'
        """

rule sort_ref_rm:
    input:
        "ref_repeatmask/{reference}.te.bed.unsorted"
    output:
        "ref_repeatmask/{reference}.te.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            bedtools sort -i '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

##### TELR Liftover

rule build_index:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    shell:
        """
        samtools faidx '{input}' || true
        touch {output}
        """

rule make_te_json:
    input:
        "contigs/{contig}/tes/{te}/00_annotation.bed"
    output:
        "contigs/{contig}/tes/{te}/00_annotation.json"
    shell:
        "python3 {config[STELR_liftover]} make_json '{input}' '{output}'"

rule get_genome_size:
    input:
        "{genome}.fai"
    output:
        "{genome}.size"
    shell:
        "python3 {config[STELR_liftover]} get_genome_size '{input}' '{output}'"

rule flank_bed:
    input:
        fasta = "contigs/{contig}/03_contig1.fa",
        contig_size = "contigs/{contig}/03_contig1.fa.size",
        te_dict = "contigs/{contig}/tes/{te}/00_annotation.json"
    output:
        "contigs/{contig}/tes/{te}/12_{flank}_flank.bed"
    params:
        flank_len = config["flank_len"]
    shell:
        """
        python3 {config[STELR_liftover]} flank_bed '{input.fasta}' '{input.contig_size}' '{input.te_dict}' '{params.flank_len}' '{output}'
        touch '{output}'
        """

rule flank_fasta:
    input:
        fasta = "contigs/{contig}/03_contig1.fa",
        bed = "contigs/{contig}/tes/{te}/12_{flank}_flank.bed"
    output:
        "contigs/{contig}/tes/{te}/12_{flank}_flank.fa"
    shell:
        """
        if [ -s '{input.bed}' ]; then
            bedtools getfasta -fi '{input.fasta}' -bed '{input.bed}' -fo '{output}'
        else
            touch '{output}'
        fi
        """

rule align_flank:
    input:
        flank_fa = "contigs/{contig}/tes/{te}/12_{flank}_flank.fa",
        ref_fa = config["reference"]
    output:
        "contigs/{contig}/tes/{te}/13_{flank}_flank.paf"
    params:
        preset = "asm10",
        num_secondary = 10
    shell:
        """
        if [ -s '{input.flank_fa}' ]; then
            minimap2 -cx '{params.preset}' -v 0 -N '{params.num_secondary}' '{input.ref_fa}' '{input.flank_fa}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule get_flank_alignment_info:
    input:
        "contigs/{contig}/tes/{te}/13_{flank}_flank.paf"
    output:
        "contigs/{contig}/tes/{te}/14_{flank}_flank.info"
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_liftover]} get_paf_info '{input}' '{output}'
        else
            touch '{output}'
        fi
        """

rule flank_paf_to_bed:
    input:
        "contigs/{contig}/tes/{te}/13_{flank}_flank.paf"
    output:
        "contigs/{contig}/tes/{te}/14_{flank}_flank.bed_unsorted"
    params:
        different_contig_name = config["different_contig_name"]
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_liftover]} paf_to_bed '{input}' '{output}' '{wildcards.contig}' '{params.different_contig_name}'
        else
            touch '{output}'
        fi
        """

rule sort_flank_bed:
    input:
        "contigs/{contig}/tes/{te}/14_{flank}_flank.bed_unsorted"
    output:
        "contigs/{contig}/tes/{te}/14_{flank}_flank.bed"
    shell:
        """
        if [ -s '{input}' ]; then
            bedtools sort -i '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule closest_flank_maps_to_ref:
    input:
        flank_5p = "contigs/{contig}/tes/{te}/14_5p_flank.bed",
        flank_3p = "contigs/{contig}/tes/{te}/14_3p_flank.bed"
    output:
        "contigs/{contig}/tes/{te}/15_flank_overlap.bed"
    shell:
        """
        if [ -s '{input.flank_5p}' ] && [ -s '{input.flank_3p}' ]; then
            bedtools closest -a '{input.flank_5p}' -b '{input.flank_3p}' -s -d -t all > '{output}'
        else
            touch '{output}'
        fi
        """

checkpoint json_for_report:
    input:
        overlap = "contigs/{contig}/tes/{te}/15_flank_overlap.bed",
        info_5p = "contigs/{contig}/tes/{te}/14_5p_flank.info",
        info_3p = "contigs/{contig}/tes/{te}/14_3p_flank.info"
    output:
        "contigs/{contig}/tes/{te}/15_flank_overlap.json"
    shell:
        """
        python3 {config[STELR_liftover]} bed_to_json {input.overlap} {input.info_5p} {input.info_3p} {output}
        touch {output}
        """

rule make_report:
    input:
        overlap = "contigs/{contig}/tes/{te}/15_flank_overlap.json",
        te_json = "contigs/{contig}/tes/{te}/00_annotation.json",
        ref_bed = lambda wildcards: f"ref_repeatmask/{os.path.basename(config['reference'])}.te.bed",
        reference = config['reference']
    output:
        "contigs/{contig}/tes/{te}/16_{overlap_id}_report.json"
    params: 
        flank_overlap_max = config["overlap"],
        flank_gap_max = config["gap"]
    shell:
        """
        python3 {config[STELR_liftover]} make_report {input.overlap} {wildcards.overlap_id} {input.te_json} {input.ref_bed} {input.reference} {params.flank_overlap_max} {params.flank_gap_max} {output}
        touch {output}
        """

def overlap_ids_report(wildcards):
    overlap_file = checkpoints.json_for_report.get(**wildcards).output[0]
    try:
        with open(overlap_file, "r") as overlap:
            overlap_dict = json.load(overlap)
            overlap_ids = [key for key in overlap_dict]
        return [f"contigs/{wildcards.contig}/tes/{wildcards.te}/16_{overlap_id}_report.json" for overlap_id in overlap_ids]
    except:
        return []
rule best_report:
    input:
        flanks = [
            "contigs/{contig}/tes/{te}/14_5p_flank.info",
            "contigs/{contig}/tes/{te}/14_3p_flank.info",
            "contigs/{contig}/tes/{te}/14_5p_flank.bed",
            "contigs/{contig}/tes/{te}/14_3p_flank.bed"
            ],
        te_json = "contigs/{contig}/tes/{te}/00_annotation.json",
        overlap_reports = overlap_ids_report
    output:
        "contigs/{contig}/tes/{te}/17_best_report.json"
    shell:
        "python3 {config[STELR_liftover]} choose_report {output} {input}"


def annotation_from_option(wildcards):
    return {True:f"contigs/{wildcards.contig}/10_annotation.bed",False:f"contigs/{wildcards.contig}/rm_05_rm_reannotated_tes.bed"}[config["minimap2_family"]]

'''3rd stage
Identify TE insertion breakpoint (minimap2)
'''

"""
4th stage: estimate intra-sample TE allele frequency (TAF)
"""

'''4th stage
Read extraction (samtools)
'''

rule read_context:
    input:
        vcf_parsed = "contigs/{contig}/00_vcf_parsed.tsv",
        bam = "reads_sort.bam"
    output:
        read_ids = "contigs/{contig}/00_read_context.id",
        vcf_parsed_new = "contigs/{contig}/00_parsed_vcf_with_readcount.tsv"
    params:
        window = 1000
    shell:
        "python3 {config[STELR_assembly]} read_context '{wildcards.contig}' '{input.vcf_parsed}' '{input.bam}' '{output.read_ids}' '{output.vcf_parsed_new}' '{params.window}'"

'''4th stage
Read alignment to TE contig (minimap2)
'''

rule get_reverse_complement:
    input:
        "contigs/{contig}/03_contig1.fa"
    output:
        "contigs/{contig}/10_revcomp.fa"
    shell:
        """
        if [ -s '{input}' ]; then
            python3 {config[STELR_utility]} get_rev_comp_sequence '{input}' '{output}'
        else
            touch '{output}'
        fi
        """
        
rule realignment:
    input:
        contig = "contigs/{contig}/{contig_revcomp}.fa",
        reads = "contigs/{contig}/00_read_context.fa"
    output:
        "contigs/{contig}/{contig_revcomp}_realign.sam"
    params:
        presets = lambda wildcards: {"pacbio":"map-pb","ont":"map-ont"}[config["presets"]]
    shell:
        """
        if [ -s '{input.contig}' ]; then
            minimap2 -a -x '{params.presets}' -v 0 '{input.contig}' '{input.reads}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule realignment_to_bam:
    input:
        "contigs/{contig}/{contig_revcomp}_realign.sam"
    output:
        "contigs/{contig}/{contig_revcomp}_realign.bam"
    shell:
        """
        if [ -s '{input}' ]; then
            samtools view -bS '{input}' > '{output}'
        else
            touch '{output}'
        fi
        """

rule sort_index_realignment:
    input:
        "contigs/{contig}/{contig_revcomp}_realign.bam"
    output:
        "contigs/{contig}/{contig_revcomp}_realign.sort.bam"
    threads: 1
    shell:
        """
        if [ -s '{input}' ]; then
            samtools sort -@ '{threads}' -o '{output}' '{input}'
            samtools index -@ '{threads}' '{output}'
        else
            touch '{output}'
        fi
        """

'''4th stage
Depth-based TAF estimation (SAMtools)
'''

rule estimate_te_depth:
    input:
        bam = "contigs/{contig}/{contig_revcomp}_realign.sort.bam",
        contig = "contigs/{contig}/{contig_revcomp}.fa"
    output:
        depth_5p = "contigs/{contig}/tes/{te}/{contig_revcomp}_5p_te.depth",
        depth_3p = "contigs/{contig}/tes/{te}/{contig_revcomp}_3p_te.depth"
    params:
        te_interval = config["af_te_interval"],
        te_offset = config["af_te_offset"]
    shell:
        """
        if [ -s '{input.bam}' ]; then
            python3 {config[STELR_te]} estimate_te_depth '{input.bam}' '{input.contig}' '{wildcards.te}' '{params.te_interval}' '{params.te_offset}' '{output.depth_5p}' '{output.depth_3p}'
        else
            touch '{output}'
        fi
        """

rule estimate_flank_depth:
    input:
        bam = "contigs/{contig}/{contig_revcomp}_realign.sort.bam",
        contig = "contigs/{contig}/{contig_revcomp}.fa"
    output:
        depth_5p = "contigs/{contig}/tes/{te}/{contig_revcomp}_5p_flank.depth",
        depth_3p = "contigs/{contig}/tes/{te}/{contig_revcomp}_3p_flank.depth"
    params:
        flank_len = config["af_flank_interval"],
        flank_offset = config["af_flank_offset"]
    shell:
        """
        if [ -s '{input.bam}' ]; then
            python3 {config[STELR_te]} estimate_flank_depth '{input.bam}' '{input.contig}' '{wildcards.te}' '{params.flank_len}' '{params.flank_offset}' '{output.depth_5p}' '{output.depth_3p}'
        else
            touch '{output}'
        fi
        """
    
rule estimate_coverage:
    input:
        te_5p = "contigs/{contig}/tes/{te}/{contig_revcomp}_5p_te.depth",
        te_3p = "contigs/{contig}/tes/{te}/{contig_revcomp}_3p_te.depth",
        flank_5p = "contigs/{contig}/tes/{te}/{contig_revcomp}_5p_flank.depth",
        flank_3p = "contigs/{contig}/tes/{te}/{contig_revcomp}_3p_flank.depth"
    output:
        "contigs/{contig}/tes/{te}/{contig_revcomp}.freq"
    shell:
        "python3 {config[STELR_te]} estimate_coverage '{input.te_5p}' '{input.te_3p}' '{input.flank_5p}' '{input.flank_3p}' '{output}'"

rule get_allele_frequency:
    input:
        fwd = "contigs/{contig}/tes/{te}/03_contig1.freq",
        rev = "contigs/{contig}/tes/{te}/10_revcomp.freq"
    output:
        "contigs/{contig}/tes/{te}/11_allele_frequency.json"
    shell:
        "python3 {config[STELR_te]} get_af '{input.fwd}' '{input.rev}' '{output}'"

'''
Write Output
'''

rule individual_json:
    input:
        liftover_file = "contigs/{contig}/tes/{te}/17_best_report.json",
        af_file = "contigs/{contig}/tes/{te}/11_allele_frequency.json",
        vcf_parsed = "contigs/{contig}/00_vcf_parsed.tsv",
        annotation_file = "contigs/{contig}/tes/{te}/00_annotation.bed",
        contig_file = "contigs/{contig}/03_contig1.fa"
    output:
        "contigs/{contig}/tes/{te}/18_output.json"
    shell:
        """
        python3 {config[STELR_output]} make_json_output {input} {output}
        touch {output}
        """


def all_contigs(wildcards):
    vcf_parsed_file = checkpoints.merge_parsed_vcf.get(**wildcards).output[0]
    with open(vcf_parsed_file) as vcf_parsed:
        contigs = []
        for line in vcf_parsed:
            contigs.append("_".join(line.split("\t")[:3]))
    return contigs
def all_tes(wildcards): #expects annotation file to be in contigs/{contig}/tes/
    tes = {}
    for contig in all_contigs(wildcards):
        annotation_file = checkpoints.annotate_contig.get(contig=contig).output[0]
        te_dir = annotation_file[:annotation_file.rindex("/")]
        with open(annotation_file, "r") as input:
            for line in input:
                entry = line.replace("\n","").split("\t")
                if len(entry) == 6:
                    tes[f"te_{entry[1]}_{entry[2]}"] = contig
    return tes
def get_output_jsons(wildcards):
    tes = all_tes(wildcards)
    return [f"contigs/{tes[te]}/tes/{te}/18_output.json" for te in tes]


rule final_output:
    input:
        reference = config["reference"],
        reference_index = lambda wildcards: f"{config['reference']}.fai",
        json_files = get_output_jsons
    output:
        contig_fa_outfile = "reads.telr.contig.fasta",
        te_fa_outfile = "reads.telr.te.fasta",
        bed_outfile = "reads.telr.bed",
        json_outfile = "reads.telr.json",
        expanded_json_outfile = "reads.telr.expanded.json",
        vcf_outfile = "reads.telr.vcf"
    shell:
        "python3 {config[STELR_output]} write_output {output} {input}"
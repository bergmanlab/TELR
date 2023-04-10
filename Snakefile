rule targets:
    input:
        file1 = "output.txt"

rule alignment:
    input:
        bam = reads.bam
        reads = reads.fasta
        reference = reference.fasta
    shell:
        "python3 TELR_alignment.py alignment [...]"
    output:
        f"{out}/{sample_name}.tmp.sam"

rule in_lib:
    shell:
        "python3 TELR_input.py input_library {tmp_dir} {library}"
    output:
        library.fasta

rule in_ref:
    shell:
        "python3 TELR_input.py input_reference {tmp_dir} {reference}"
    output:
        reference.fasta

rule input: # generalized?
    shell:
        "python3 TELR_input.py input {config[{type}]} {tmp_dir} {input}"
    output:
        "inputs/{type}.fasta" #make inputs folder duh

rule bam_input:
    input:
        reads.bam
        shell:
        #bam2fasta
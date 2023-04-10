rule targets:
    input:
        file1 = "output.txt"

rule alignment:
    input:
        bam = reads.bam
        reads = reads.fasta
        reference = reference.fasta
    execute:
        python3 TELR_alignment.py alignment [...]
    output:
        f"{out}/{sample_name}.tmp.sam"

rule in_lib:
    execute:
        python3 TELR_input.py input_library {tmp_dir} {library}
    output:
        library.fasta

rule in_ref:
    execute:
        python3 TELR_input.py input_reference {tmp_dir} {reference}
    output:
        reference.fasta

rule input: # generalized?
    execute:
        python3 TELR_input.py input {config[{type}]} {tmp_dir} {input}
    output:
        {type}.fasta

rule bam_input:
    input:
        reads.bam
    execute:
        #bam2fasta
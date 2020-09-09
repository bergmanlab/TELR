# Usage

## Required arguments
TELR requires long reads in FASTA/FASTQ format or read alignments in BAM format, a reference genome assembly in FASTA format (which must be the same as the one used to align the reads, if read alignments are provided), and TE consensus sequence in FASTA format. Here is a template with the names of the required and optional parameters.
```
python3 telr.py -i (--reads) <reads in fasta/fastq format or read alignments in bam format> \
        -r (--reference) <reference genome in fasta format> \
        -l (--library) <TE consensus sequences in fasta format>
```

## Optional arguments
In addition to the required program options listed above, there are some optional arguments. The full list of program options with descriptions can also be obtained by running `telr.py -h`.
- `-x (--presets) <arg>` Preset for different sequencing technologies (default = 'pacbio').
- `-p (--polish) <int>` Rounds of contig polishing using polisher from [wtdbg2](https://github.com/ruanjue/wtdbg2) (default = 1).
- `-o (--out) <arg>` Output directory (default = '.').
- `-t (--thread) <int>` Maximum cpu threads to use (default = '1').
- `-g (--gap) <int>` Maximum gap size between sequence alignments of two contig flanks (default= '20').
- `-v (--overlap) <int>` Maximum overlap size between sequence alignments of two contig flanks (default= '20').
- `-k (--keep_files)` If provided then all intermediate files will be kept (default: remove intermediate files).
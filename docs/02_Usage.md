# Usage

## Command line help page
```
usage: telr.py [-h] -i READS -r REFERENCE -l LIBRARY [-x PRESETS] [-p POLISH]
               [-o OUT] [-t THREAD] [-g GAP] [-v OVERLAP] [-k]

Script to detect TEs in long read data

required arguments:
  -i READS, --reads READS
                        reads in fasta/fastq format or read alignments in bam
                        format
  -r REFERENCE, --reference REFERENCE
                        reference genome in fasta format
  -l LIBRARY, --library LIBRARY
                        TE consensus sequences in fasta format

optional arguments:
  -h, --help            show this help message and exit
  -x PRESETS, --presets PRESETS
                        parameter presets for different sequencing
                        technologies (default = 'pacbio')
  -p POLISH, --polish POLISH
                        rounds of contig polishing (default = 1)
  -o OUT, --out OUT     directory to output data (default = '.')
  -t THREAD, --thread THREAD
                        max cpu threads to use (default = '1')
  -g GAP, --gap GAP     max gap size for flanking sequence alignment (default
                        = '20')
  -v OVERLAP, --overlap OVERLAP
                        max overlap size for flanking sequence alignment
                        (default = '20')
  -k, --keep_files      If provided then all intermediate files will be kept
                        (default: remove intermediate files)
```

## Required arguments
TELR requires long reads in FASTA/FASTQ format or read alignments in BAM format from BWA-MEM (use `-M` and `-x` parameter), Minimap2 (with Cigar & MD string) or NGMLR, a reference genome assembly in FASTA format (which must be the same as the one used to align the reads, if read alignments are provided), and TE consensus sequence in FASTA format. Here is a template with the names of the required and optional parameters.
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

<p align="center">
    <img src="https://github.com/bergmanlab/TELR/blob/master/TELR.png?raw=true" alt="TELR"/>
</p>

## Introducation
TELR (pronounced Teller) is a fast non-reference transposable element (TE) detector from long read sequencing data (PacBio or Oxford Nanopore). It uses 

## Install dependencies
You can use conda to install dependencies and create running environment for TELR.
```
conda create -n TELR_env python=3.6 -y
conda activate TELR_env
conda install -c conda-forge -y biopython
conda install -y pandas=1.0.0
conda install -y repeatmasker=4.0.7
conda install -y samtools=1.9
conda install -y bcftools=1.9
conda install -y bedtools
conda install -y ngmlr=0.2.7
conda install -y sniffles=1.0.11
conda install -y wtdbg
conda install -y seqtk
conda install -y minimap2
conda install -y svim=1.3
conda install -c bioconda -y mummer
conda install -c conda-forge gnuplot
```

## Quick Start
```
git clone git@github.com:bergmanlab/TELR.git
cd TELR
./telr.py --help
```

## Running on a test dataset
```
cd test
../telr.py -o test_output -i reads.fasta -r ref_38kb.fasta -l library.fasta
```

## Parameters
```
usage: telr.py [-h] -i READ -r REFERENCE -l LIBRARY [-x PRESETS] [-o OUT]
               [-t THREAD] [-g GAP] [-p OVERLAP]

Script to detect TEs in long read data

required arguments:
  -i READ, --read READ
      reads in fasta/fastq format
  -r REFERENCE, --reference REFERENCE
      reference genome in fasta format
  -l LIBRARY, --library LIBRARY
      TE consensus sequences in fasta format

optional arguments:
  -h, --help
      show this help message and exit
  -x PRESETS, --presets PRESETS
      parameter presets for different sequencing technologies (default = 'pacbio')
  -p, --polish
      provide rounds of contig polishing (default = 1)
  -o OUT, --out OUT
      directory to output data (default = '.')
  -t THREAD, --thread THREAD
      max cpu threads to use (default = '1')
  -g GAP, --gap GAP
      max gap size for flanking sequence alignment (default = '20')
  -v OVERLAP, --overlap OVERLAP
      max overlap size for flanking sequence alignment (default = '20')
```

## Output
The results of TELE pipeline are output to the directory <output>.
- `<sample>.final.bed`: non-reference TE insertion annotation predicted by TELR pipeline in bed format (0-based).
- `<sample>.final.fa`: TE insertion sequences from local assembly of reads supporting TE insertions.
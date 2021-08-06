<p align="center">
    <img src="https://github.com/bergmanlab/TELR/blob/master/img/TELR.png?raw=true" alt="TELR"/>
</p>

## Introduction
TELR (pronounced Teller) is a fast non-reference transposable element (TE) detector from long read sequencing data (PacBio or Oxford Nanopore). TELR uses long reads mapped to a reference genome to identify insertions using [Sniffles](https://github.com/fritzsedlazeck/Sniffles), then filters insertions by matching insertion supporting reads with user supplied TE consensus sequences. For each TE insertion candidate locus, TELR performs a local assembly of all reads supporting TE insertion, annotates the TE sequence in assembled contigs, then maps the flanks back to the reference genome. Finally, TELR generates the coordinates of the non-reference TE insertions, the estimated allele frequency and the assembled TE sequences.

The TELR pipeline consists of four main stages: (1) general SV detection and filter for TE insertion candidate, (2) local reassembly and polishing of the TE insertion, (3) identification of TE insertion coordinates, and (4) estimation of intra-sample TE insertion allele frequency.

- In stage 1, long reads are aligned to the reference genome using NGMLR (https://github.com/philres/ngmlr). The alignment output in BAM format is provided as input for Sniffles (https://github.com/fritzsedlazeck/Sniffles) to detect structural variations (SVs). TELR then filter for TE insertion candidates from SVs reported by Sniffles using following criteria: 1) The type of SV is insertion. 2) Insertion sequence is available. 3) The insertion sequences include hits from user provided TE consensus library using RepeatMasker (http://www.repeatmasker.org}).

- In stage 2, reads that support the TE insertion candidate locus based on Sniffles output are used as input for wtdbg2 (https://github.com/ruanjue/wtdbg2) to assemble local contig that covers the TE insertion for each TE insertion candidate locus. The local assemblies are then polished using minimap2 and wtdbg2.

- In stage 3, TE consensus library is aligned to the assembled TE insertion contigs using minimap2 and used to define TE-flank boundaries. TE region in each contig is annotated with family info using RepeatMasker. Sequences flanking the TE insertion are then re-aligned to the reference genome using minimap2 to determine the precise TE insertion coordinates and target site duplication (TSD).

- In stage 4, raw reads aligned to the reference genome are extracted within a 1kb interval on either side of the insertion breakpoints initially defined by Sniffles. The reads are then aligned to the assembled polished contig to identify reads that support the TE insertion and reference alleles, respectively, in following steps: 1) Reads are aligned to the forward strand of the contig, 5' flanking sequence depth (5p_flank_cov) and 5' TE depth (5p_te_cov) are calculated. 2) Reads are aligned to the reverse complement strand of the contig, 5' flanking sequence depth (3p_flank_cov) and 5' TE depth (3p_te_cov) are calculated. 3) The TE allele frequency is estimated as (5p_te_cov/5p_flank_cov + 3p_te_cov/3p_flank_cov)/2.

<p align="center">
<img src="https://github.com/bergmanlab/TELR/blob/master/img/TELR_workflow.png?raw=true"/>
</p>

The current version of TELR shows good performance on real Drosophila melanogaster data sets, including datasets with heterozygous TE insertions.

TELR is written in python3 and is designed to run on linux operating system.

## Documentation
The following sections will provide you installation instructions, usage guide, and descriptions of output files.
  - [Installation](docs/01_Installation.md)
  - [Usage](docs/02_Usage.md)
  - [Output Files](docs/03_Output_Files.md)

## Getting Help
Please use the [Github Issue page](https://github.com/bergmanlab/TELR/issues) if you have questions.
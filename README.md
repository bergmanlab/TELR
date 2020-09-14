<p align="center">
    <img src="https://github.com/bergmanlab/TELR/blob/master/TELR.png?raw=true" alt="TELR"/>
</p>

## Introducation
TELR (pronounced Teller) is a fast non-reference transposable element (TE) detector from long read sequencing data (PacBio or Oxford Nanopore).

TELR uses long reads mapped to a reference genome to identify insertions using [Sniffles](https://github.com/fritzsedlazeck/Sniffles), then filters insertions by matching insertion supporting reads with user supplied TE consensus sequences. For each TE insertion candidate locus, TELR performs a local assembly of all reads supporting TE insertion, annotates the TE sequence in assembled contigs, then maps the flanks back to the reference genome. Finally, TELR generates the coordinates of the non-reference TE insertions plus the assembled TE sequences.

The current version of TELR shows good performance on real drosophila melanogaster data sets with lots of heterozygous TE insertions. We are currently doing more testing on other species.

## Documentation
The following sections will provide you installation instructions, usage guide, and descriptions of output files.
  - [Installation](docs/01_Installation.md)
  - [Usage](docs/02_Usage.md)
  - [Output Files](docs/03_Output_Files.md)

## Getting Help
Please use the [Github Issue page](https://github.com/bergmanlab/TELR/issues) if you have questions.

## License

Copyright (c) 2020 Shunhua Han and Casey M. Bergman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

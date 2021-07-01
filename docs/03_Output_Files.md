# Output
TELR outputs non-referece TE insertion predictions in multiple format.
- `<sample>.telr.vcf`: non-reference TE insertion predictions in VCF format (1-based).
- `<sample>.telr.json`: non-reference TE insertion predictions in JSON format (0-based).
- `<sample>.telr.bed`: non-reference TE insertion predictions in BED format (0-based).
- `<sample>.telr.fa`: TE insertion sequences assembled by TELR.

## VCF file output by TELR
TELR generates a standard VCF file `<sample>.telr.vcf` in [v4.1 format](https://samtools.github.io/hts-specs/VCFv4.1.pdf) that has detailed information for each non-reference TE insertion.

Column | Description
-- | --
chromosome | The chromosome name where the TE insertion occurred
position | Starting breakpoint position of the TE insertions.
ID | The id of the TE insertions.
Ref | The sequence of the reference is always set to N.
Alt | The TE insertion sequence.
Quality | This is currently not indicated.
Filter | This is currently always set to be PASS
Info | Provides a list of information (see below)
FORMAT | Provides information about the next tag
Sample information | Depending on the way sniffles was run: Genotype estimation:Reads   supporting the reference: Reads supporting the variant.

### Info field description
Sniffles report multiple information in the Info field. The entries are delimited by ;.

INFO key | Description
-- | --
SVTYPE= | The type of the variant, currently this is always set to INS
END= | The position of the second breakpoint of the TE insertion
FAMILY= | TE families of the insertion, multiple families are separated by '\|'
STRANDS= | Strand that TE insertion occurs
SUPPORT_TYPE= | Type of support from flank to reference alignment (1: only one flank used to infer the insertion; 2: two flanks with gap; 3: two flanks with overlap)
RE= | Number of reads supporting the TE insertion
AF= | Allele frequency of the variant

## JSON file output by TELR
The VCF and JSON files are essentially equivalent, but the JSON file `<sample>.telr.json` can be easier to parse programmatically. For each non-reference TE insertion, the JSON file contain these keys:

Key | Description
-- | --
ID | The id of the TE insertions
chr | The chromosome name where the TE insertion occurred
start | Starting breakpoint position of the TE insertions
end | The position of the second breakpoint of the TE insertion
family | TE families of the insertion, multiple families are separated by '\|'
support_type | Type of support from flank to reference alignment (1: only one flank used to infer the insertion; 2: two flanks with gap; 3: two flanks with overlap)
strand | Strand that TE insertion occurs
frequency | Allele frequency of the variant
sequence | The TE insertion sequence
gt | Genotype of the variant
alt_count | Number of reads supporting the TE insertion

## BED file output by TELR
The BED file includes minimal info about non-reference TE insertions but is easier to use in subsquent bioinformatics analysis. For each non-reference TE insertion, the BED file contain these fields:

Column | Description
-- | --
Chromosome | The chromosome name where the TE insertion occurred
Start | Starting breakpoint position of the TE insertions
End | The position of the second breakpoint of the TE insertion
Family | TE families of the insertion, multiple families are separated by '\|'
Score | Allele frequency of the variant
Strand | Strand that TE insertion occurs

## FASTA file output by TELR
For each non-reference TE insertion, TELR reports insertion sequences in `<sample>.telr.fa`.

## Log file output by TELR
For each TELR run, a log file called `<sample>.log` is generated that records all the major steps in the program and error messages.
`<sample>.log`: log file of TELR run.

## Debugging info output by TELR
For each TELR run, a tabular file `<sample>.loci_eval.tsv` is generated that report why each potential insertion locus from SV detection is filtered out in TELR pipeline. The messages are very self-explanatory as listed below.
- VCF sequence not repeatmasked
- Contig assembly failed
- Sniffles VCF sequence not mapped to assembled contig
- VCF sequence doesn't overlap contig annotation
- No flanks mapped to reference
- No flanks mapped to correct chromosome
- No flanks have unique hit
- Two flanks mapped to different chromosomes or strands
- contigs without RM annotation
- Overlap/gap between contigs flanks exceeds threshold
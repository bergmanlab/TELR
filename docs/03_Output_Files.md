# Output
TELR outputs non-referece TE insertion predictions in multiple format.
- `<sample>.telr.vcf`: non-reference TE insertion predictions in VCF format (1-based).
- `<sample>.telr.json`: non-reference TE insertion predictions in JSON format (0-based).
- `<sample>.telr.expanded.json`: non-reference TE insertion predictions in JSON format (0-based) with expanded information.
- `<sample>.telr.bed`: non-reference TE insertion predictions in BED format (0-based).
- `<sample>.telr.fasta`: TE insertion sequences assembled by TELR.
- `<sample>.telr.contig.fasta`: Local contig sequences assembled by TELR that include new TE insertions.

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
SUPPORT_TYPE= | Type of support from flank to reference alignment (single_side or both_sides)
RE= | Number of reads supporting the TE insertion
AF= | Allele frequency of the variant
TSD_LEN= | Length of the TSD sequence if available
TSD_SEQ= | TSD sequence if available

## JSON file output by TELR
The VCF and JSON files are essentially equivalent, but the JSON file `<sample>.telr.json` can be easier to parse programmatically. For each non-reference TE insertion, the JSON file contain these keys:

Key | Description
-- | --
type | The type of the TE insertions (non-reference only in current version)
ID | The unique id of the TE insertions
chrom | The chromosome name where the TE insertion occurred
start | Starting breakpoint position of the TE insertions
end | The position of the second breakpoint of the TE insertion
family | TE families of the insertion, multiple families are separated by '\|'
strand | Strand that TE insertion occurs
support | Type of support from flank to reference alignment (single_side or both_sides)
tsd_length | Length of the TSD sequence if available
tsd_sequence | TSD sequence if available
te_sequence | The TE insertion sequence
genotype | Genotype of the variant
num_sv_reads | Number of reads supporting the SV allele
num_ref_reads | Number of reads supporting the reference allele
allele_frequency | Allele frequency of the variant

## Expanded JSON file output by TELR
Comapred to the basic JSON file described above, the expanded JSON file `<sample>.telr.expanded.json` includes more QC metrics that could help with filtering TE sequences for subsquent analysis. For each non-reference TE insertion, the expanded JSON file contain these additional keys compared to the basic JSON file:

Key | Description
-- | --
gap_between_flank | The size of the gap between 3' and 5' flanking sequence alignment to the reference genome (the value is negative if two alignments overlap)
te_length | The length of the new TE insertion sequence
contig_id | Unique ID for the local contig assembly
contig_length | The length of the local contig assembly
contig_te_start | Starting position of new TE insertion in the contig assembly
contig_te_end | End position of new TE insertion in the contig assembly
5p_flank_align_coord | Coordinate of 5' flanking sequence alignment to the reference genome
5p_flank_mapping_quality | Mapping quality of 5' flanking sequence alignment
5p_flank_num_residue_matches | Number of residue matches of 5' flanking sequence alignment
5p_flank_alignment_block_length | Alignment block length of 5' flanking sequence alignment
5p_flank_sequence_identity | Sequence identity of 5' flanking sequence alignment (calculated as 5p_flank_num_residue_matches/5p_flank_alignment_block_length)
3p_flank_align_coord | Coordinate of 3' flanking sequence alignment to the reference genome
3p_flank_mapping_quality | Mapping quality of 3' flanking sequence alignment
3p_flank_num_residue_matches | Number of residue matches of 3' flanking sequence alignment
3p_flank_alignment_block_length | Alignment block length of 3' flanking sequence alignment
3p_flank_sequence_identity | Sequence identity of 3' flanking sequence alignment (calculated as 3p_flank_num_residue_matches/3p_flank_alignment_block_length)

## BED file output by TELR
The BED file includes minimal info about non-reference TE insertions but is easier to use in subsquent bioinformatics analysis. For each non-reference TE insertion, the BED file contain these fields:

Column | Description
-- | --
Chromosome | The chromosome name where the TE insertion occurred
Start | Starting breakpoint position of the TE insertions
End | The position of the second breakpoint of the TE insertion
Family | TE families of the insertion, multiple families are separated by '\|'
Score | '.'
Strand | Strand that TE insertion occurs

## TE FASTA file output by TELR
For each non-reference TE insertion, TELR reports new TE insertion sequences in `<sample>.telr.te.fasta`.

## Contig FASTA file output by TELR
For each non-reference TE insertion, TELR reports assembled contig sequences in `<sample>.telr.contig.fasta`.

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
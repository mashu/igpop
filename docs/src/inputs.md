# Input files

## Alleles (`--alleles`)

A FASTA file. The record name up to the first whitespace is the allele name,
and the text before `*` is the gene: `IGHV1-69*01` belongs to gene
`IGHV1-69`. Sequences may be on either strand, may be partial, and may
contain IMGT gap dots (removed) and `N` (never reported as a variant).

## Coordinates (`--coords`)

Optional for human IGH: without `--coords`, igpop uses the shipped table
`data/coordinates.tsv` and the aliases in `data/gene_aliases.tsv`.

A tab-separated table that says where each gene lies, six columns; lines
starting with `#` and a header line starting with `gene` are ignored.

| Column | Meaning |
|---|---|
| gene | Gene name, matching the allele names. |
| build | `GRCh38`, `GRCh37` (or `38`, `hg38`, `37`, `hg19`). |
| contig | Chromosome, with or without `chr`. |
| start, end | 1-based inclusive coordinates of the gene. |
| strand | `+`, `-`, or `.` to let igpop choose: both orientations are aligned and the better-scoring one is used. |

A table may hold rows for several builds. A run uses the rows of `--build`,
and for a gene without such a row the row of another build.

Several genes may share exactly the same interval (for example `IGHV3-23`
and `IGHV3-23D`). igpop treats them as one locus and says so in the
`shares_locus_with` column.

### Files shipped in `data/`

| File | Content |
|---|---|
| `coordinates.tsv` | Human IGH loci: 123 on GRCh38 (110 lifted from the IgSNPer GRCh37 coordinates, 12 native, IGHV4-38-2 on an unplaced scaffold) and 5 rows on GRCh37, used for the genes GRCh38 lacks. |
| `lift_report.tsv` | Evidence for every lift (see below). |
| `gene_aliases.tsv` | IGHV4-4 allele-to-locus assignments from IgSNPer. |

`lift_report.tsv` records how each GRCh37 locus was placed on GRCh38: the
gene with 1 kb of flank on each side was searched in the GRCh38 IGH region
by edit distance. A locus was moved only if the placement was unique and
the gene itself was at least 95% identical. Columns: `gene`, `source`
build, `GRCh37` and `GRCh38` intervals, number of `placements`,
`flank_identity`, `gene_edit_distance`, `gene_identity`, `decision`.

## Gene aliases (`--gene-aliases`)

`allele<TAB>gene`, one per line. Use it when alleles of one named gene
belong to different loci, as with `IGHV4-4*01` (locus `IGHV4-4B`) and
`IGHV4-4*07` (locus `IGHV4-4A`).

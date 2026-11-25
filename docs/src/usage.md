# Usage and options

```bash
igpop variants --alleles alleles.fasta -o results/ighv
```

With no further options igpop uses the shipped human IGH coordinates and
gene aliases, genome build GRCh38, and the default remote datasets. See
[Quick start](quickstart.md).

A run writes nine files, described in [Output files](outputs.md):

```
results/ighv.report.md               read this first: results and checks in words
results/ighv.final.tsv               alleles that passed
results/ighv.alleles.tsv             every allele: verdict, why it failed, haplotype check
results/ighv.variants.tsv            one row per allele variant, all evidence
results/ighv.alignments.txt          alignments to read by eye, one block per gene
results/ighv.alignments.tsv          where each allele was placed
results/ighv.score_model.tsv         the curve behind the score
results/ighv.score_levels.tsv        the two groups counted per level of support
results/ighv.haplotype_calibration.tsv  what a found haplotype is worth
results/ighv.dbsnp_datasets.tsv      which dbSNP frequency datasets were used
```

## Options

### What to analyse

| Option | Default | Meaning |
|---|---|---|
| `--alleles` | required | Allele FASTA. |
| `--coords` | shipped human IGH table | Coordinates table, see [Input files](inputs.md). |
| `--build` | `38` | `38`/`GRCh38`/`hg38` or `37`/`GRCh37`/`hg19`. Genes are placed on this build if it has them. |
| `--place-unknown` | off | Place an allele whose gene is not in the coordinates table on the locus its sequence fits best, with a note, instead of leaving it `not_placed`. |
| `--single-build` | off | Do not judge genes that are missing from `--build` on another build; they are `not_placed`. |
| `--gene-aliases` | shipped table (with the shipped coordinates) | `allele<TAB>gene` table for alleles whose locus differs from the gene in their name. |
| `-o`, `--out` | `igpop` | Output prefix. |

### Score and verdict

| Option | Default | Meaning |
|---|---|---|
| `--candidate-marker` | `_S` | Alleles whose name contains this text are candidates, all others are named trusted. |
| `--trusted` | none | File of trusted allele names, one per line; replaces the marker rule. |
| `--calibration` | automatic | Use the `score_model.tsv` of an earlier run. See [Automatic calibration](calibration.md). |
| `--common` | `0.01` | Allele frequency a variant must demonstrably reach in some population to be classed `common`. |
| `--supported-score` | `2` | A variant with at least this score is `supported`. |
| `--confirm-score` | `0.1` | A variant with at most this score is `needs_confirmation`. |
| `--coincidence` | `0.05` | An allele found as an exact haplotype is `supported` if fewer than this share of arbitrary changes are found as often. |
| `--prior` | none | Share of variants like yours expected to be real (between 0 and 1); fills `probability_real`. |

### Data sources

| Option | Default | Meaning |
|---|---|---|
| `--reference` | 1000 Genomes FASTA for the build | Path or URL of a FASTA with `.fai` (and `.gzi` if bgzipped). |
| `--vcf` | 1000 Genomes 30x annotated calls (GRCh38) | Site-level VCF with the `*_unrel` INFO fields; `{contig}` is replaced by the contig name. |
| `--gnomad` | gnomAD v4.1 genomes (GRCh38) | gnomAD sites VCF; `{contig}` as above. |
| `--no-gnomad` | off | Do not use gnomAD. |
| `--no-frequencies` | off | Read no population data at all: alignments and variants only. |
| `--phased-vcf` | 1000 Genomes 30x phased panel (GRCh38) | Phased VCF with genotypes for the [haplotype check](haplotypes.md); `{contig}` as above. |
| `--pedigree` | 1000 Genomes sample table | Family, sample, father, mother, sex, population, superpopulation. |
| `--no-haplotypes` | off | Skip the haplotype check. |
| `--dbsnp` | latest dbSNP for the build | dbSNP VCF: rsIDs, COMMON flag and relayed frequencies. |
| `--no-dbsnp` | off | Do not use dbSNP: no rsIDs, no COMMON flag, no relayed frequencies. |
| `--dbsnp-exclude` | `SGDP_PRJ` | dbSNP frequency datasets, comma-separated, to leave out of the evidence. |
| `--cache-dir` | `.igpop_cache` | Where downloaded VCF regions are kept. |
| `--no-cache` | off | Neither read nor write the cache. |

### Alignment

| Option | Default | Meaning |
|---|---|---|
| `--flank` | `100` | Reference bases added on each side of a locus. |
| `--min-identity` | `0.9` | Alignments below this are reported but give no variants. |
| `--no-paralog-check` | off | Do not score each allele against every other locus (faster). |

## Downloads and the cache

The log ends with the amount downloaded, for example
`downloaded 63.53 MB in 163 requests`. Fetched VCF regions are stored in
`--cache-dir`, so a second run with the same loci downloads only the
reference window. Delete the directory to force a fresh download, for
example after a dataset release changes.

Set `IGPOP_TRACE=1` to print every request, its size and duration, and any
retry. Requests that fail with HTTP 429, 5xx or a dropped connection are
retried up to five times with exponential backoff, and at most three
connections are open at once.

## More than one build

A gene that has coordinates only on another build is judged on that build,
in the same run (see
[Datasets](datasets.md#genes-that-the-reference-genome-lacks)). Files
given with `--reference`, `--vcf`, `--gnomad`, `--dbsnp` and `--phased-vcf`
belong to `--build`; the other build uses its defaults. `--no-gnomad`,
`--no-dbsnp` and `--no-haplotypes` apply to both. `--single-build` keeps
the run to `--build`.

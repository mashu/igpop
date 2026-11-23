# Datasets and combining evidence

## What is used

| Dataset | Build | Samples | Role |
|---|---|---|---|
| 1000 Genomes 30x (NYGC), annotated site-level calls | GRCh38 | 3,202 sequenced; counts taken from the 2,504 unrelated | Frequencies per superpopulation, consortium site filters. |
| gnomAD v4.1 genomes | GRCh38 | 76,215 | Frequencies per genetic-ancestry group, gnomAD filters. |
| dbSNP (latest release) | GRCh38 / GRCh37 | — | rsIDs, the COMMON flag, and frequencies relayed from other projects (TOPMed, ALFA, gnomAD exomes, ExAC, ToMMo, Korean panels, ...), used as evidence. |
| 1000 Genomes 30x phased panel | GRCh38 | 2,598 founders | Haplotype check only. |

An rsID by itself is never evidence: dbSNP records every submitted
change. Only reported frequencies count.

The related 698 samples of 1000 Genomes are excluded from counts because
children repeat their parents' alleles; they enter only through the
Mendelian error rate.

## Which dataset a variant comes from

Every row shows it three ways:

- `evidence` lists the sources in which the variant was observed;
- `class_basis` names the dataset and population that decided the class;
- the `kgp_*`, `gnomad_*` and `dbsnp_*` columns hold each source's own
  numbers.

Nothing is merged into a single frequency.

## Why counts are not added together

Adding allele counts is the correct way to combine two samples only if they
are **independent** and drawn from the **same population**. Neither holds
here:

- **Overlap.** gnomAD genomes include the 1000 Genomes samples. Summing
  would count the same people twice.
- **Different populations.** The datasets group people differently (five
  superpopulations versus nine ancestry groups) and in different
  proportions. A pooled frequency would mostly reflect who was sampled.
- **Shared systematic error.** Both are short reads mapped to GRCh38. A
  paralog artefact appears in both, so agreement between them does not
  cancel it. More samples shrink sampling error only.

## How evidence is combined

Two results use the datasets together, and neither adds counts across
datasets.

**The score** summarises all datasets into one number per variant (its
*support*: the largest lower frequency bound over every population of
every dataset) and compares that number between trusted variants and
arbitrary base changes. Because both groups are measured through the same
datasets, overlap and shared errors affect both alike. See
[The score](score.md).

**The class** asks a yes/no question: is there a population in which the
variant demonstrably reaches `--common`? Every population of every dataset
is tested with an exact binomial test, and the best p-value is multiplied
by the number of populations tested. Overlap between datasets makes this
correction stricter than necessary, never too lenient.

## What the second dataset adds

- **Replication.** A variant observed in both (`evidence` lists two names)
  was called by two pipelines in mostly different people.
- **Less cell-line bias.** Most gnomAD genomes are not from B-cell lines,
  so its frequencies are not distorted by V(D)J deletions in the way 1000
  Genomes frequencies near the D–J region are.
- **Depth for rare variants.** About 150,000 chromosomes instead of about
  5,000: a variant absent from gnomAD is bounded far more tightly.

Disagreement is informative too. When the 1000 Genomes and gnomAD
frequencies of a variant differ strongly, compare the population columns:
either the variant is concentrated in a population one dataset samples
poorly, or one of the two sites is flagged.

## Frequencies relayed by dbSNP

Besides rsIDs, a dbSNP record carries a `FREQ` field with allele
frequencies submitted by other projects: ALFA (`dbGaP_PopFreq`), TOPMed,
gnomAD exomes, ExAC, the Japanese ToMMo panel, Korean panels, ESP, and
smaller national cohorts. Several are far larger than 1000 Genomes, and
they add people that 1000 Genomes and gnomAD genomes cover poorly. IgSNPer
based its verdict on six of them (`dbGaP_PopFreq`, `ExAC`, `TOPMED`,
`TOMMO`, `KOREAN`, `GoESP`), taking the first one in that order with a
frequency above 0.0001.

igpop uses the same datasets as evidence, in the same way as the
populations of 1000 Genomes and gnomAD:

- **No priority order.** Every usable dataset is one more population in
  which the variant is tested; the best one decides, with the correction
  for the number of populations tested. This stays valid when datasets
  overlap (ExAC is part of gnomAD exomes): overlap makes the correction
  conservative, not wrong.
- **Allele-specific.** A frequency counts only if it is for exactly the
  allele's base change.
- **Dataset size matters.** dbSNP gives a frequency and no counts. igpop
  takes the size of a dataset as one over the smallest non-zero frequency
  it reports in the regions read (one copy among all its chromosomes), and
  the count of an allele as frequency times size. The size estimate cannot
  exceed the real size. The same exact test and confidence bounds then
  apply as for 1000 Genomes: 3 copies in 300 chromosomes do not pass,
  3,000 in 300,000 do.
- **An approximation, stated as such.** The number of chromosomes
  genotyped at one particular site is not published; where a dataset
  covers a site poorly, the reconstructed count is too generous. This is
  the price of using datasets that publish frequencies only. To leave any
  of them out, name them in `--dbsnp-exclude`.
- **The COMMON flag is honoured.** A variant that dbSNP flags as COMMON is
  classed `common` with `class_basis = dbSNP:COMMON`, but only for the most
  frequent alternate allele of that site.

`class_basis` shows which population decided (`dbSNP:TOPMED`,
`gnomAD_v4.1:afr`, ...), and `dbsnp_freq` lists every reported frequency.

Which datasets are used is decided from the data, and written to
`<prefix>.dbsnp_datasets.tsv` together with the estimated size:

1. A dataset that re-reports a source igpop read directly in the same run
   is left out (`1000Genomes`, `1000Genomes_30X`, `GnomAD_genomes`). On a
   build where igpop has no primary source they are kept.
2. A dataset must be large enough that seeing nothing would bound the
   frequency below `--common`: at least `3 / common` chromosomes.
3. `SGDP_PRJ` is excluded by default (`--dbsnp-exclude`): almost half of
   its values in IGH are exactly 0.5, i.e. one genome per site.

In the IGHV example the used datasets are gnomAD exomes, TOPMed, ExAC,
ToMMo, ALFA, ESP, Korea4K, Estonian, ALSPAC, TwinsUK, KOREAN, GoNL,
Northern Sweden, MGP, Vietnamese and FINRISK.

## The phased panel

The 1000 Genomes 30x phased panel supplies individual chromosomes for the
[haplotype check](haplotypes.md). It does not enter the score or the
verdict.

## Other builds

1000 Genomes 30x and gnomAD v4.1 do not exist on GRCh37, and neither does
the phased panel. With `--build 37` the evidence is what dbSNP relays for
that build (including its copies of 1000 Genomes and gnomAD), scored with
the built-in calibration. There is no haplotype check. Both primary sources
are plain options (`--vcf`, `--gnomad`), so another callset with the same
INFO fields can be supplied.

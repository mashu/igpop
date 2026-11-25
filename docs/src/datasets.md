# Datasets and combining evidence

## What is used

| Dataset | Build | Samples | Role |
|---|---|---|---|
| 1000 Genomes 30x (NYGC), annotated site-level calls | GRCh38 | 3,202 sequenced; counts taken from the 2,504 unrelated | Frequencies per superpopulation, consortium site filters. |
| gnomAD v4.1 genomes | GRCh38 | 76,215 | Frequencies per genetic-ancestry group, gnomAD filters. |
| dbSNP (latest release) | GRCh38 / GRCh37 | — | rsIDs, the COMMON flag, and frequencies relayed from other projects (TOPMed, ALFA, gnomAD exomes, ExAC, ToMMo, Korean panels, ...), used as evidence. |
| ALFA (NCBI), release 4 | GRCh38 | differs by position, up to several hundred thousand | Counts in nine populations: European, African American, African Others, East Asian, South Asian, Other Asian, Latin American 1 and 2, Other. |
| GenomeAsia pilot | GRCh38 | 1,739 | Counts in seven groups: Northeast, Southeast and South Asian, Oceanian, African, American, West Eurasian. |
| ToMMo 61KJPN | GRCh38 | about 61,000 | Counts in Japanese people. |
| Westlake BioBank for Chinese (WBBC) | GRCh38 | 4,480 | Counts in four regions of China: North, Central, South, Lingnan. |
| 1000 Genomes 30x phased panel | GRCh38 | 2,598 founders | Haplotype check only. |

The last four are read for GRCh38 only and can be switched off together
(`--no-more-datasets`, or *Read more populations* in the web app). They
add populations the first two cover thinly or as one block: Japanese,
regions of China, Southeast Asia, and the populations of ALFA one by one
(dbSNP passes on only its total). ALFA and ToMMo are also among the
studies dbSNP relays; when they are read directly, their dbSNP copies are
left out, so that the same people are not consulted twice. A dataset that
cannot be reached is reported and left out; the run continues.

Not included: datasets whose licence does not allow use by a public tool
(for example SG10K_Health of Singapore, the only large panel with a Malay
group), and datasets that publish frequencies without counts.

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

## Why not dbSNP alone

dbSNP already passes on frequencies from 1000 Genomes and gnomAD, so why
read those two directly? Because dbSNP gives one overall frequency per
dataset and nothing else. Reading the files themselves adds:

| | From dbSNP | Read directly |
|---|---|---|
| Populations | One frequency for the whole dataset. | Counts for the five continental groups of 1000 Genomes and the nine ancestry groups of gnomAD. |
| Counts | A frequency without the number of chromosomes; igpop estimates the size of each study. | Copies and chromosomes examined, which the tests need. |
| A variant with no record | Nobody submitted it; nothing follows. | It was looked for and not seen: 0 of about 5,000 or 150,000 chromosomes. |
| Site quality | None. | Filters, missing genotypes, Hardy-Weinberg and Mendelian errors (1000 Genomes); filter and region flags (gnomAD). |
| Whole sequences | None. | The phased genotypes of 1000 Genomes, for the [haplotype check](haplotypes.md). |

Populations matter most. Many real IGHV variants are common in one
population and rare overall. rs370240779 in IGHV1-18 is carried by 478 of
40,476 African/African American chromosomes in gnomAD (1.2%: common, and
supported), but by only 516 of 149,388 overall (0.35%). From the overall
figure alone it is `inconclusive`.

In the other direction, dbSNP adds rsIDs and studies igpop does not read
itself: TOPMed, the gnomAD exomes, ToMMo, Korean and several European
cohorts. Where igpop reads 1000 Genomes or gnomAD directly, it leaves out
dbSNP's copy of that dataset, so that the same people are not counted
twice.

On the 519 IGHV alleles of the example that differ from the reference:

| Sources read | Alleles supported |
|---|---:|
| 1000 Genomes, gnomAD and dbSNP | 439 |
| without gnomAD | 423 |
| without dbSNP | 427 |
| 1000 Genomes only | 404 |

No single source is essential, and each contributes. The web app reads
gnomAD from the project's copy at Amazon, because the server at Google
does not allow requests from web pages; with gnomAD switched off in its
options, its results correspond to the second line.

## If a source cannot be read

Every request is retried (five times with growing pauses on the command
line, three times in the browser). After that:

| Source | The run | Signalled by |
|---|---|---|
| dbSNP | continues without rsIDs, the COMMON flag and relayed frequencies | a warning above the results and in the report |
| Phased panel | continues without the haplotype check (`hap_status = not_checked`) | a warning above the results and in the report |
| 1000 Genomes frequencies | stops | an error message |
| gnomAD | stops (in the web app, switch gnomAD off under Options and run again) | an error message |
| Reference sequence | stops: alleles cannot be placed without it | an error message |

A run that continued without a source says so only in its warnings: the
verdict of a single allele does not show it. Read the warnings before the
alleles, and run again when a source was missing.

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
- **The COMMON flag is honoured, with care.** dbSNP sets the flag on a
  record, and a record can list several alternate alleles. igpop gives it
  to one allele only: the one with the highest frequency in any study of
  the record that is large enough to show a frequency of 1%. A study of a
  few people, in which every allele it saw appears at one half, does not
  decide. The flagged allele is classed `common` with
  `class_basis = dbSNP:COMMON`, unless the data igpop reads directly
  (1000 Genomes, gnomAD) show the allele to be rare in every population:
  first-hand counts outrank a second-hand flag.

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

## Genes that the reference genome lacks

Some IGHV genes lie on stretches of DNA that only some people carry. The
GRCh38 reference represents a chromosome without IGHV1-8, IGHV3-9, IGHV4-31
and the IGHV4-4*07-like gene (IGHV4-4A in the shipped table): their
sequences are nowhere in GRCh38, so their variants cannot be carried over
from another build either. GRCh37 has these genes. igpop therefore judges
each gene on the build that has it, in one run:

| | GRCh38 (all other genes) | GRCh37 (genes missing from GRCh38) |
|---|---|---|
| 1000 Genomes | 30x call set, 2,504 unrelated samples | phase 3, the same 2,504 people at low coverage |
| gnomAD genomes | v4.1, 76,215 genomes, 9 groups | v2.1.1, 15,708 genomes, 6 groups |
| dbSNP | GRCh38 file | GRCh37 file |
| Haplotype check | phased 30x panel, 5,196 founder chromosomes | phase 3 phased genotypes, 5,008 chromosomes |
| Site quality | reported | not published for phase 3 |

The score is calibrated on the GRCh38 genes only and applied to both. The
GRCh37 data are smaller, so the same variant gets at most the same support
there: the verdicts for these genes err towards `inconclusive`.

Read them with more caution all the same. People who lack the gene
contribute no reads, and reads of a similar gene can be placed on it by
mistake. The output marks every such allele: column `build`, and a note
"gene absent from GRCh38: judged on GRCh37 ...".

`--single-build` switches this off; the genes are then `not_placed`.

IGHV4-38-2 is a third case. GRCh38 has it, but only on
`chr14_KI270726v1_random`, a piece of sequence that is not part of a
chromosome. The 1000 Genomes call set reports 227 variant sites in the 296
bases of the gene there and gives no counts per population; gnomAD and the
phased panel do not cover it. igpop compares the alleles with that
sequence (an identical allele is `matches_reference`) and reads no
population data: an allele that differs is `not_assessed`.

With `--build 37` the whole run uses the GRCh37 sources of the table.
Both primary sources are plain options (`--vcf`, `--gnomad`), so another
call set with the same INFO fields can be supplied for the main build.

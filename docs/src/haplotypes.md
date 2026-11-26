# The haplotype check

The frequency evidence judges every difference of an allele on its own. Two
differences can each be frequent and still never occur on the same
chromosome, as in a chimera of two real alleles. And a single rare
difference cannot be told from an error by its frequency alone. The
haplotype check asks the direct question: **do chromosomes with exactly
this sequence exist, and in whom?**

## What is done

igpop reads the phased 1000 Genomes 30x panel for the gene (genotypes of
3,202 people, each chromosome separately) and counts chromosomes of
**founders**: the 2,598 people without a parent in the panel, 5,196
chromosomes. Children of the trios are left out because their chromosomes
are copies of their parents'.

For each allele it counts

- **`hap_carrying`**: chromosomes that carry all of the allele's variants;
- **`hap_exact`**: of those, chromosomes with no other variant anywhere
  inside the stretch the allele covers, i.e. the same sequence as the
  allele; and the number of people carrying them (`hap_exact_donors`);

in total and per superpopulation (AFR, AMR, EAS, EUR, SAS).

## `hap_status`

| Status | Meaning | How to read it |
|---|---|---|
| `found_exact` | At least one chromosome has exactly the allele's sequence. | Support for the allele as a whole; see below for how many chromosomes count. |
| `found_with_extra` | The combination exists, but every carrier chromosome has further variants inside the allele. | The variants are real together; the allele as written may be incomplete, or the extra variants are panel artefacts. Compare with the alignment. |
| `not_found` | Every variant is in the panel, but no chromosome carries them all. | Evidence against this combination among 5,196 chromosomes. With a truly rare allele this is expected. |
| `incomplete` | Some variants are not sites of the panel (`hap_untested_str`). `hap_carrying` refers to the remaining ones; there is no `hap_exact`. | No answer for the full combination. |
| `untestable` | None of the variants is a site of the panel. | No answer. |
| `not_checked` | Allele not placed, or no phased panel for this build. | No answer. |

For an allele identical to the reference (`matches_reference`),
`hap_exact` counts chromosomes with no variant inside the allele.

## How much is a "found" worth? The calibration

A sequencing or PCR error turns a real allele into a sequence one base
away. Sometimes that sequence exists too, by coincidence. igpop measures
how often: it takes every founder chromosome of the panel, changes one base
at a time, at every position of every gene of the run, and looks the
changed sequence up among the founder chromosomes.

In the IGHV example:

| The changed sequence is found on at least | Share of arbitrary changes |
|---:|---:|
| 1 chromosome | 3.8% |
| 2 chromosomes | 1.7% |
| 3 chromosomes | 0.9% |
| 5 chromosomes | 0.5% |
| 10 chromosomes | 0.3% |
| 50 chromosomes | 0.1% |

This share is the column `hap_coincidence` of an allele: how often an
arbitrary change does as well as this allele did.

**An allele is supported by its haplotype** when

- its exact sequence is found (`found_exact`),
- on at least two founder chromosomes (a single one can be a genotyping
  error), and
- fewer than 5% of arbitrary changes are found that often
  (`--coincidence 0.05`).

In the IGHV example that means two chromosomes or more. The table is
recomputed in every run and written to `<prefix>.haplotype_calibration.tsv`.

Such an allele gets the verdict `supported`, with `supported_by` =
`haplotype` (or `frequency+haplotype`), even when its single differences
are too rare to pass by frequency. Of the 566 IGHV alleles, 71 are
supported by haplotype only and 91 by both.

Alleles with two or more differences that pass this test serve a second
purpose: they join the named trusted alleles as the yardstick of the
[score](calibration.md), because they are real alleles found by the data,
rare ones included.

## When the check speaks against an allele

Frequencies are counted one difference at a time. They can show that every
difference of an allele exists in people; they cannot show that the
differences occur together. A chimera of two real alleles, made in a PCR
reaction, consists only of real variants.

So when the panel has every difference of an allele and still no
chromosome with exactly its sequence (`found_with_extra` or `not_found`),
an allele that would be supported by frequency is reported as
`inconclusive`, with the filter `sequence_not_in_panel`. Among 324
artificial chimeras of two database alleles of the same gene, 56 could be
looked up; 22 of them were on no chromosome, and none of those is reported
as supported.

## What the haplotype check cannot do

- **It often has no answer.** The phased panel contains only sites that
  pass the consortium's filters (missing genotypes below 5%,
  Hardy-Weinberg, Mendelian errors at most 5%) and leaves out variants
  seen in a single chromosome. In IGHV many real, frequent variants sit at
  sites that fail these filters. In the example, 240 of 566 alleles
  get `found_exact`, 4 `found_with_extra`, none `not_found`, and 322 cannot
  be fully tested. `incomplete` and `untestable` say nothing about the
  allele.
- **For one difference it is not new evidence.** For an allele with a
  single difference, "exact haplotype found" and "the variant has this
  frequency in 1000 Genomes" are largely the same fact. The check adds the
  most for alleles with several differences.
- **The comparison is with single-base errors.** Other artefacts are not
  arbitrary: reads of a duplicated gene assigned to the wrong gene produce
  the same pseudo-variant in repertoire data and in short-read population
  data. Look at `notes` (ambiguous or shared locus) for such alleles.
- **Same people as the 1000 Genomes counts.** The panel is not a second
  cohort; it answers a different question on the same genomes.
- `not_found` is bounded: nothing in 5,196 chromosomes still allows a
  frequency up to about 0.07% (`hap_exact_AF_upper95`).

## Population specificity

| Column | Meaning |
|---|---|
| `hap_exact_AFR` ... `hap_exact_SAS` | Carrier chromosomes / founder chromosomes in each superpopulation. |
| `hap_exact_max_pop` | Superpopulation with the highest lower frequency bound. |
| `hap_exact_pop_heterogeneity_p` | Test of "same frequency in all five superpopulations". A small value means the haplotype is population-specific. With fewer than about ten carriers the test is approximate. |

Example: `IGHV1-24*01_S1381` is found on 72 chromosomes of 64 people, 70 of
the chromosomes African (70/1380 AFR, 1/706 AMR, 0 EAS, 1/1052 EUR, 0 SAS;
heterogeneity p = 1e-36).

Population specificity is reported, not used for the verdict. It shows
that the signal is biological and not a technical artefact, which would
appear in all populations alike; but old, real variants are spread evenly
too.

## Caveats

- Phasing is statistical. Inside one gene (a few hundred bases) switch
  errors are rare, but rare variants are phased less reliably than common
  ones.
- "No other variant inside the allele" covers the allele's own stretch
  only, not the flanks. A truncated allele is compared over its length.
- Structural variants of the panel that overlap the gene count as "other
  variants".
- The panel exists for GRCh38 only.

## Options

| Option | Meaning |
|---|---|
| `--coincidence` | Share of arbitrary changes that may do as well as an allele supported by haplotype (default 0.05). Lower is stricter. |
| `--no-haplotypes` | Skip the check (saves about 30 MB of download for a full IGHV set). |
| `--phased-vcf` | Another phased VCF with genotypes; `{contig}` is replaced. |
| `--pedigree` | Sample table for it: family, sample, father, mother, sex, population, superpopulation. |

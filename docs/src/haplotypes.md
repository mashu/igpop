# The haplotype check

The verdict judges every difference of an allele on its own. Two
differences can each be frequent and still never occur on the same
chromosome. The haplotype check asks the stricter question: **is the whole
combination carried together on one chromosome, and in whom?**

It is extra evidence. It does not change the verdict; it is reported next
to it (`hap_*` columns of `alleles.tsv`, a column in the report lists) so
that you can weigh it yourself.

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
  allele;

in total and per superpopulation (AFR, AMR, EAS, EUR, SAS).

## `hap_status`

| Status | Meaning | How to read it |
|---|---|---|
| `found_exact` | At least one chromosome has exactly the allele's sequence. | Support for the allele as a whole. Look at `hap_exact`: dozens of chromosomes are strong evidence; one or two can also arise from an imputed genotype. |
| `found_with_extra` | The combination exists, but every carrier chromosome has further variants inside the allele. | The variants are real together; the allele as written may be incomplete, or the extra variants are panel artefacts. Compare with the alignment. |
| `not_found` | Every variant is in the panel, but no chromosome carries them all. | Evidence against this combination among 5,196 chromosomes. With a truly rare allele this is expected. |
| `incomplete` | Some variants are not sites of the panel (`hap_untested_str`). `hap_carrying` refers to the remaining ones; there is no `hap_exact`. | No answer for the full combination. |
| `untestable` | None of the variants is a site of the panel. | No answer. |
| `not_checked` | Allele not placed, or no phased panel for this build. | No answer. |

For an allele identical to the reference (`matches_reference`),
`hap_exact` counts chromosomes with no variant inside the allele.

## Why so many alleles are `incomplete` or `untestable`

The phased panel is deliberately strict. It contains only sites that pass
the consortium's filters (missing genotypes below 5%, Hardy-Weinberg,
Mendelian errors at most 5%) and leaves out variants seen in a single
chromosome. In IGHV many real, frequent variants sit at sites that fail
these filters, because duplicated genes and gene deletions disturb
short-read genotyping. In the IGHV example, 209 of 513 placed alleles get
`found_exact`, 3 `found_with_extra`, none `not_found`, and 301 cannot be
fully tested.

So:

- `found_exact` with a good count is positive evidence you can lean on.
- `incomplete` and `untestable` say nothing about the allele. They say the
  panel has no opinion.
- `not_found` is the only negative outcome, and it is bounded: nothing in
  5,196 chromosomes still allows a frequency up to about 0.07%
  (`hap_exact_AF_upper95`).

## Population specificity

| Column | Meaning |
|---|---|
| `hap_exact_AFR` ... `hap_exact_SAS` | Carrier chromosomes / founder chromosomes in each superpopulation. |
| `hap_exact_max_pop` | Superpopulation with the highest guaranteed frequency. |
| `hap_exact_pop_heterogeneity_p` | Test of "same frequency in all five superpopulations". A small value means the haplotype is population-specific. |

Example: `IGHV1-24*01_S1381` is found on 72 chromosomes, 70 of them African
(70/1380 AFR, 1/706 AMR, 0 EAS, 1/1052 EUR, 0 SAS; heterogeneity
p = 1e-36). A candidate seen in a donor of matching ancestry fits; the same
candidate "discovered" only in donors of other ancestry would be worth a
second look.

A population-specific haplotype is normal biology, not a warning.

## How independent is this evidence?

Partly. The phased panel comes from the same 3,202 genomes as the 1000
Genomes counts used for the verdict, so it is not a second cohort. What is
new is the question: per-variant frequencies cannot tell whether variants
travel together; phased chromosomes can. gnomAD and the datasets relayed by
dbSNP publish no individual genotypes, so the check cannot be repeated
there.

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
| `--no-haplotypes` | Skip the check (saves about 30 MB of download for a full IGHV set). |
| `--phased-vcf` | Another phased VCF with genotypes; `{contig}` is replaced. |
| `--pedigree` | Sample table for it: family, sample, father, mother, sex, population, superpopulation. |

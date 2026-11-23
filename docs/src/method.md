# Method

## 1. Placing an allele

For each gene, igpop reads the reference sequence of the given interval plus
`--flank` bases on each side. The allele is aligned to that window with a
semi-global affine-gap alignment (Gotoh): the allele is aligned end to end,
the window may overhang on both sides. Scores are +2 match, −3 mismatch,
−5 gap open, −2 gap extend; `N` scores 0. No read mapper is involved.

Two checks guard the placement:

- **Identity.** Below `--min-identity` the allele is reported as
  `low_identity` and yields no variants.
- **Competing loci.** The allele is also scored against every other locus in
  the table. If another locus scores at least as well, the status is
  `ambiguous_locus`: the sequence alone does not say which gene it belongs
  to, so its variants are reported but should be read with that in mind.

## 2. Variants

Every difference becomes a `(position, REF, ALT)` record in VCF convention.
Substitutions are single bases. Insertions and deletions are left-aligned
and trimmed to their shortest form (Tan, Abecasis & Kang 2015), and VCF
records are split per ALT and trimmed the same way, so that the same event
has the same representation on both sides.

Reference positions the allele does not cover are unknown; igpop does not
assume the reference base there.

## 3. Frequencies and bounds

For each population, the allele count `AC` out of `AN` chromosomes gives
exact (Clopper-Pearson) one-sided 95% bounds on the allele frequency.

The bounds make datasets of different sizes comparable. Three copies among
200 chromosomes and 300 among 20,000 are both 1.5%, but only the second
supports the claim "at least 1%".

## 4. Score and class

The **score** compares the variant's support with that of arbitrary base
changes and of trusted variants; it is described in [The score](score.md).

The **class** is a formal yes/no statement. In every population where the
allele was seen at least twice, an exact binomial test asks whether the
data are compatible with a frequency below `--common`. The smallest
p-value is multiplied by the number of populations tested (Bonferroni). At
0.05 or below the variant is `common`. Variants seen but not passing are
`uncommon`; variants never seen are `absent`.

Testing the best-supported population rather than a pooled frequency
follows gnomAD's "filtering allele frequency" idea. A pooled frequency
depends on how many samples each population contributed, which is a
property of the study design, not of the variant.

## 5. Population enrichment

`kgp_pop_heterogeneity_p` is the p-value of a likelihood-ratio (G) test that
the frequency is equal across the five superpopulations (4 degrees of
freedom). It answers "is this variant enriched somewhere?" separately from
"is it common?".

## 6. Site quality

Whether a frequency can be trusted is a separate question from how large it
is, and igpop reports it in separate columns. See
[Site quality filters](site-qc.md).

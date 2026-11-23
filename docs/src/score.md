# The score

Every single-base variant gets a **score**: how many times more often its
level of population support is found among variants of trusted alleles than
among arbitrary base changes.

- a score of **100** means "this much support is 100 times more typical of
  real variants than of arbitrary changes";
- a score of **1** means the evidence does not point either way;
- a score of **0.01** means "this is what an arbitrary change looks like":
  the evidence speaks about 100 to 1 against.

## Why a count of uncommon variants is not enough

Counting has two problems.

**1. "Seen in a database" is weak evidence, and how weak depends on the
database.** Large databases contain a huge number of rare variants. In the
IGHV genes, about one in four of *all possible* single-base changes has
been seen at least once in gnomAD. An invented base change therefore has a
good chance of being "known". A fixed frequency cut-off does not fix this:
any cut-off is too strict for a small dataset and too lenient for a large
one.

**2. All "uncommon" variants are counted the same, but they are not the
same.** These five alleles each have exactly one uncommon variant
(real results, KI IGHV set):

| Allele | The uncommon variant | Seen in 1000 Genomes | Seen in gnomAD | Score |
|---|---|---|---|---|
| IGHV5-51\*01_S3648 | no rsID | 0 of 5,008 | 0 of 152,378 | 0.006 |
| IGHV1-18\*01_S2706 | rs368426916 | 0 of 5,008 | 1 of 150,172 | 0.017 |
| IGHV1-18\*03 | rs373407210 | 2 of 4,936 | 39 of 150,064 | 0.084 |
| IGHV1-18\*01_S1946 | rs370240779 | 5 of 4,938 | 264 of 149,388 | 1.2 |
| IGHV3-13\*01_S2991 | rs535166996 | 24 of 4,918 | 83 of 149,588 | 17.5 |

By count they are identical: one uncommon variant each. By score they span
a 3,000-fold range, from "looks like an error" to "well supported, just short
of the formal 1% test". The count and the score agree in direction (alleles
with more uncommon variants have lower scores), which is why both are
reported; the score adds *how* uncommon.

## How the score is made

Think of a diagnostic test. A positive result means little if healthy
people often test positive too. What matters is how much more often
patients test positive than healthy people. The score is that ratio.

1. **Measure the support of a variant.** In every population of 1000
   Genomes and gnomAD, and in every dataset relayed by dbSNP, igpop computes the lowest frequency the
   data are compatible with (the lower 95% confidence bound). The largest
   of these is the variant's *support*. A bound is used instead of the raw
   frequency so that 3 copies in 200 chromosomes do not count like 300 in
   20,000.
2. **See what arbitrary changes look like.** igpop takes every possible
   single-base change in the gene intervals of the run (three per position)
   and measures its support the same way. This is the "healthy" group: most
   of these changes are not real.
3. **See what real variants look like.** It does the same for the variants
   carried by trusted alleles. This is the "patient" group.
4. **Divide.** Support is sorted into seven levels. For each level, the
   score is the share of trusted variants at that level divided by the
   share of arbitrary changes at that level.

For the KI IGHV set this gives (`<prefix>.score_calibration.tsv`):

| Support (lowest compatible frequency) | Arbitrary changes | Trusted variants | Score |
|---|---|---|---|
| none | 49.1% | 0% | 0.006 |
| below 0.001% | 14.3% | 0% | 0.017 |
| 0.001–0.01% | 20.4% | 0% | 0.017 |
| 0.01–0.1% | 10.2% | 0.6% | 0.084 |
| 0.1–1% | 4.2% | 4.7% | 1.2 |
| 1–10% | 0.9% | 15.2% | 17.5 |
| 10% or more | 0.9% | 79.5% | 83 |

Read the last line as: 80% of trusted variants, but only 0.9% of arbitrary
changes, have this much support, so a variant at this level is about 80
times more likely to belong to the first group. (A level with no trusted
variant still gets a small score above zero, because half a count is added
everywhere.)

## What makes this work across datasets

- **Dataset size is handled automatically.** A big dataset sees many
  arbitrary changes at low frequency, so low-frequency sightings in it earn
  a low score. Nobody has to choose a cut-off per dataset.
- **Overlap between datasets does no harm.** gnomAD contains the 1000
  Genomes samples, and all short-read datasets share mapping errors.
  Multiplying evidence across datasets, as a naive Bayes classifier does,
  would count the same evidence several times. Here the datasets are first
  summarised into one number, and both groups are measured through the same
  datasets, so shared quirks affect both sides of the ratio equally.
- **No weights are set by hand.** A source in which arbitrary changes look
  as good as real ones produces a ratio near 1 and contributes nothing.

## Details

- **Trusted alleles.** By default every placed allele whose name does not
  contain `_S` (`--candidate-marker`); or give a list with `--trusted`.
  At least 30 distinct trusted single-base variants are needed.
- **Small allele sets.** D and J sets have too few variants to calibrate;
  igpop then uses its built-in human IGHV table, or one you pass with
  `--calibration`. See [Automatic calibration](calibration.md).
- **From score to verdict.** At least 10: `supported`. At most 0.1:
  `needs_confirmation`. In between: `inconclusive`. Ten to one is the
  usual mark for strong evidence in either direction.
- **From score to probability.** Probability = prior x score /
  (prior x score + 1 - prior), where the prior is the share of such
  variants you expect to be real (`--prior`).
- **Smoothing.** Half a count is added to every level so that an empty
  level gives neither zero nor infinity. More support can never mean a
  lower score; neighbouring levels that would break this are pooled
  (`pooled_group` in the table).
- **Per allele.** `min_score` is the lowest score among the allele's
  variants: an allele is only as well supported as its weakest variant.
  Scores of several variants are not multiplied, because variants on the
  same chromosome are inherited together and are not independent evidence.
- **Not scored.** Insertions and deletions (there is no "every possible
  indel" to compare with).

## What the score does not say

- It measures **population support**, not truth. A real allele private to
  one family scores low. A mapping artefact present in every dataset can
  score high; the [site quality](site-qc.md) columns address that.
- It is **relative to your trusted set**. If the trusted alleles are
  unusual, so is the scale.
- It is coarse by design: seven levels, each backed by enough variants to
  estimate a ratio.

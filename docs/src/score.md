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
| IGHV5-51\*01_S3648 | no rsID | 0 of 5,008 | 0 of 152,378 | 0.005 |
| IGHV1-18\*01_S2706 | rs368426916 | 0 of 5,008 | 1 of 150,172 | 0.013 |
| IGHV1-18\*03 | rs373407210 | 2 of 4,936 | 39 of 150,064 | 0.34 |
| IGHV1-18\*01_S1946 | rs370240779 | 5 of 4,938 | 264 of 149,388 | 2.0 |
| IGHV3-13\*01_S2991 | rs535166996 | 24 of 4,918 | 83 of 149,588 | 22 |

By count they are identical: one uncommon variant each. By score they span
a 4,000-fold range, from "looks like an error" to "well supported, just short
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
| none | 49.1% | 0% | 0.005 |
| below 0.001% | 14.3% | 0% | 0.013 |
| 0.001–0.01% | 20.4% | 0% | 0.013 |
| 0.01–0.1% | 10.2% | 3.3% | 0.34 |
| 0.1–1% | 4.2% | 8.5% | 2.0 |
| 1–10% | 0.9% | 19.2% | 22 |
| 10% or more | 0.9% | 69.0% | 72 |

Read the last line as: 69% of the variants of real alleles, but only 0.9%
of arbitrary changes, have this much support, so a variant at this level
is about 70 times more likely to belong to the first group. (A level with no trusted
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

## What is fixed, what is measured, and what is a convention

Three things are easy to confuse.

- **The levels are fixed.** Support is sorted into seven levels: never
  seen, and six steps of ten from below 0.001% to 10% or more. They are
  the same in every run.
- **The score of each level is measured.** A run counts how many variants
  of its real alleles and how many of all possible single-base changes
  fall into each level, and divides the two shares. In the IGHV example:

  | Level | Real variants (of 213) | Possible changes (of 38,316) | Score | Roughly between |
  |---|---:|---:|---:|---|
  | never seen | 0 | 18,812 | 0.005 | 0 and 0.08 |
  | below 0.01% (two levels pooled) | 0 | 13,297 | 0.013 | 0 and 0.2 |
  | 0.01 to 0.1% | 7 | 3,899 | 0.34 | 0.17 and 0.69 |
  | 0.1 to 1% | 18 | 1,614 | 2.0 | 1.3 and 3.1 |
  | 1 to 10% | 41 | 332 | 22 | 16 and 30 |
  | 10% or more | 147 | 362 | 72 | 63 and 83 |

  The last column is a 95% interval. A level with few real variants has a
  score that is known only roughly.
- **The limits 10 and 0.1 are conventions.** A ratio of ten or more is the
  usual reading of "strong evidence"; nothing in the data singles it out.
  They can be changed (`--supported-score`, `--confirm-score`). In this
  example no level scores between 2.0 and 22, so any limit from 3 to 20
  gives the same verdicts.

So a score of 2 is not a rejection. It says: at this level of support
there are twice as many real variants as arbitrary changes, proportionally.
That speaks mildly for the difference and is reported as `inconclusive`,
which means "the population data do not settle it".

## From a score to a probability

How likely a difference is to be real depends on the score and on what
was known before looking at population data: how reliable the method is
that produced the allele, in how many people it was seen.

```text
odds after = odds before x score
```

| Believed before | Score 2 | Score 5 | Score 22 |
|---|---:|---:|---:|
| 50% | 67% | 83% | 96% |
| 80% | 89% | 95% | 99% |
| 95% | 97% | 99% | 99.8% |

An allele found in several people by a careful method starts high, and a
score of 2 leaves it high. The verdict does not use such a prior, because
igpop cannot know it; `--prior` adds the probability as a column, and the
web app has a slider for it in the evidence of each difference.

## Seeing what a difference was compared with

In the web app, "populations and reasons" on a difference shows the
comparison itself: every real variant of the run as a dot at its support
(point at it for its rsID and gene), all possible changes as bars, the
seven levels with their counts and scores, and the difference marked among
them.

# The score

Every single-base variant gets a **score**: how many times more typical its
population support is of the variants of trusted alleles than of arbitrary
base changes.

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
| IGHV1-18\*01_S2706 | rs368426916 | 0 of 5,008 | 1 of 150,172 | 0.005 |
| IGHV1-18\*03 | rs373407210 | 2 of 4,936 | 39 of 150,064 | 0.48 |
| IGHV1-18\*01_S1946 | rs370240779 | 5 of 4,938 | 264 of 149,388 | 3.8 |
| IGHV3-13\*01_S2991 | rs535166996 | 24 of 4,918 | 83 of 149,588 | 7.2 |

By count they are identical: one uncommon variant each. By score they span
a more than 1,000-fold range, from "looks like an error" to "fairly well
supported, just short of the formal 1% test". The count and the score agree in direction (alleles
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
4. **Compare the two groups at every support.** The score of a support is
   how much more densely the real variants lie there than the arbitrary
   changes. igpop estimates this as a smooth curve, not in steps:

   - A change that **no dataset has seen** is a case of its own. Its score
     is the share of real variants never seen divided by the share of
     arbitrary changes never seen.
   - For a change that **was seen**, a logistic regression on
     x = log10(support) estimates, at every x, the odds that a seen change
     is one of the real variants rather than one of the arbitrary changes.
     Multiplied by the ratio of the two group sizes, these odds are the
     ratio of the two densities, which is the score:
     log(score) = b0 + b1 x + b2 x².

For the KI IGHV set (213 real variants, 38,316 possible changes, of which
18,812 were never seen) the curve is
log(score) = 4.96 + 1.11 x − 0.22 x², and a change never seen scores 0.005.
Some points on it:

| Support (lowest compatible frequency) | Score | 95% interval | Arbitrary changes supported at least as well |
|---|---:|---:|---:|
| never seen | 0.005 | up to 0.08 | 100% |
| 0.01% | 0.05 | 0.02 to 0.13 | 16% |
| 0.03% | 0.19 | 0.10 to 0.35 | 15% |
| 0.1% | 0.70 | 0.47 to 1.1 | 6.0% |
| 0.3% | 2.1 | 1.6 to 2.9 | 3.5% |
| 1% | 6.4 | 4.8 to 8.6 | 1.8% |
| 3% | 16 | 12 to 21 | 1.5% |
| 10% | 38 | 31 to 47 | 0.9% |
| 30% | 75 | 63 to 90 | 0.6% |

Read the line for 3% as: a variant shown to reach 3% in some population is
16 times more typical of real variants than of arbitrary changes, and only
1.5% of all possible changes are supported that well. The score reaches 10
at a support of 1.7% and falls to 0.1 at 0.018%.

The curve itself is written to `<prefix>.score_model.tsv`; the two groups
counted in seven levels of support, for reading, to
`<prefix>.score_levels.tsv`.

### Why a curve and not levels

Earlier versions sorted support into seven fixed levels (steps of ten) and
gave every variant of a level the same score. That is easy to read, but a
level lumps together variants that differ ten-fold in support: with levels,
a variant at 0.9% scored 2 and one at 1.1% scored 22, although the data
hardly distinguish them. On the curve they score 5.9 and 6.9. The levels
are still shown, as a ruler; they no longer decide anything.

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
  igpop then uses its built-in human IGHV curve, or one you pass with
  `--calibration`. See [Automatic calibration](calibration.md).
- **From score to verdict.** At least 10: `supported`. At most 0.1:
  `needs_confirmation`. In between: `inconclusive`. Ten to one is the
  usual mark for strong evidence in either direction.
- **From score to probability.** Probability = prior x score /
  (prior x score + 1 - prior), where the prior is the share of such
  variants you expect to be real (`--prior`).
- **Shape of the curve.** More support can never mean a lower score. The
  curved term (x²) is kept only if it improves the fit (Akaike's
  criterion) and the curve still rises over the whole range; otherwise a
  straight line in x is used. A seen change never scores below one that
  was never seen. For the never-seen score half a count is added to both
  groups, so that an empty group gives neither zero nor infinity.
- **How well the score is known.** `score_low` and `score_high` are a 95%
  interval from the uncertainty of the fitted curve. It is narrow where
  many real variants lie and wide where few do (below about 0.05%).
- **How often chance does as well.** `chance` is the share of all possible
  single-base changes in the genes of the run whose support is at least
  that of the variant.
- **Per allele.** `min_score` is the lowest score among the allele's
  variants: an allele is only as well supported as its weakest variant.
  Scores of several variants are not multiplied, because variants on the
  same chromosome are inherited together and are not independent evidence.
- **Insertions and deletions borrow the curve.** The two groups are
  single-base changes, because every possible one can be listed. An
  insertion or deletion gets the score the curve gives its support, as a
  lower bound (`score_basis` is `borrowed`). The reason: of 988 artificial
  single-base insertions and deletions, 2.8% were found in any dataset and
  0.1% had a support of 0.1% or more, against 54% and 7% of artificial
  base changes. An arbitrary insertion or deletion is far less likely to
  be "known", so the same support is at least as strong evidence. A lower
  bound can speak for a variant but not against it: where the curve is at
  1 or below, the insertion or deletion stays without a score and its
  verdict rests on its class (`absent` is `needs_confirmation`, otherwise
  `inconclusive`).
- **Matched as the exact change.** Every difference is looked up by
  position, reference bases and alternate bases together, after both the
  allele's difference and the records of the population files are brought
  to the same left-aligned, shortest form. A base change recorded at the
  position of a deletion is a different variant with its own rsID and
  counts, and of a record that lists several alternate alleles only the
  one that is the same change is used.

## What the score does not say

- It measures **population support**, not truth. A real allele private to
  one family scores low. A mapping artefact present in every dataset can
  score high; the [site quality](site-qc.md) columns address that.
- It is **relative to your trusted set**. If the trusted alleles are
  unusual, so is the scale.
- Below the support of the rarest real variant (0.018% in the IGHV
  example) the curve is an extrapolation; the interval shows it.

## What is measured and what is a convention

- **The curve is measured.** It comes from the real variants and the
  possible changes of the run (or, for small runs, from the built-in IGHV
  curve). Nothing about it is set by hand.
- **The limits 10 and 0.1 are conventions.** A ratio of ten or more is the
  usual reading of "strong evidence" for likelihood ratios; nothing in the
  data singles it out. They can be changed (`--supported-score`,
  `--confirm-score`). Because the score is smooth, a limit is a line drawn
  across a slope, not a step: a variant scoring 9 is nearly as well
  supported as one scoring 11, and the table shows the number so that you
  can read it that way. In the IGHV example a limit of 5 instead of 10
  supports 409 instead of 393 alleles, and a limit of 3 supports 422; the
  share of artificial substitutions that pass rises from about 1.5% to
  about 1.9% and 2.6%. The web app has a slider for the limit that shows
  these consequences while you move it.

So a score of 3 or 5 is not a rejection. It says that such support is
three or five times more typical of real variants than of arbitrary
changes. That speaks for the difference without being strong evidence, and
is reported as `inconclusive`, which means "the population data do not
settle it".

## From a score to a probability

How likely a difference is to be real depends on the score and on what
was known before looking at population data: how reliable the method is
that produced the allele, in how many people it was seen.

```text
odds after = odds before x score
```

| Believed before | Score 0.1 | Score 2 | Score 5 | Score 20 |
|---|---:|---:|---:|---:|
| 50% | 9% | 67% | 83% | 95% |
| 80% | 29% | 89% | 95% | 99% |
| 95% | 66% | 97% | 99% | 99.7% |

An allele found in several people by a careful method starts high, and a
score of 2 leaves it high. The verdict does not use such a prior, because
igpop cannot know it. It is one number for a whole run, a judgement of the
method that produced the alleles: `--prior` adds the probability as a
column, and the web app has one optional setting for it above the list of
alleles (off unless you turn it on).

## Seeing what a difference was compared with

In the web app, "show" in the Details column of a difference shows the
comparison itself on one axis of support: the score curve with its
interval and the two limits, every real variant of the run as a dot (point
at it for its rsID and gene), all possible changes as bars, and the
difference marked on the curve.

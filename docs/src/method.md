# Method

igpop answers one question about a germline allele: **do the differences
between this allele and the reference genome exist in people?** It does not
discover alleles and it does not read sequencing reads. It takes an allele
that somebody inferred, places it on the gene the user names, and looks up
each difference in population data that were produced independently.

This page says what happens, in order, which statistical method is used at
each step and why, and how the method relates to existing ones. The pages
on [the score](score.md), [the haplotype check](haplotypes.md), [the
datasets](datasets.md) and [site quality](site-qc.md) go deeper into single
steps.

## The steps of a run

| | Step | Result |
|---|---|---|
| 1 | Place each allele on its gene | an alignment to the reference, or a reason why not |
| 2 | List the differences | one `(position, reference bases, allele bases)` per difference |
| 3 | Read the population data of those genes | per difference and per group of people: copies seen, chromosomes examined |
| 4 | Turn counts into confidence bounds | per group: the lowest and highest frequency the counts are compatible with |
| 5 | Summarise in one number | the **support** of the difference |
| 6 | Compare with real and with arbitrary changes | the **score** |
| 7 | Test for "common" | the **class**: `common`, `uncommon`, `absent`, `not_assessed` |
| 8 | Decide per difference | `supported`, `inconclusive`, `needs_confirmation`, `not_assessed` |
| 9 | Look for the whole sequence on phased chromosomes | the **haplotype check** |
| 10 | Decide per allele | the verdict of the allele and its reason |

Steps 1 and 2 use only the allele and the reference. Steps 3 to 5 use only
the population data. Steps 6 to 10 use only the numbers of the steps before
them. Site quality (how far the calls at a position can be trusted) is read
in step 3 and reported next to every verdict; it does not change one.

## 1. Placing an allele

**The user says where each gene is.** Every run has a gene list with the
coordinates of each gene (the built-in list for human IGHV, IGHD and IGHJ,
or the user's own). An allele named `IGHV1-69*04` is compared with the
interval of IGHV1-69 and nowhere else.

A general-purpose aligner would choose the place by similarity. In the
immunoglobulin loci that is the wrong choice to leave to a program: many
genes are near-copies of each other (IGHV1-69 and IGHV1-69D, IGHV3-23 and
IGHV3-23D, IGHV3-30, IGHV3-30-3, IGHV3-30-5 and IGHV3-33), an allele of one
gene can be more similar to the reference copy of another, and some genes
are not in the reference at all. The gene an allele belongs to is knowledge
the person who inferred it has. igpop takes it as given and then checks it.

For each gene igpop reads the reference sequence of the interval plus
`--flank` bases on each side and aligns the allele to it: the whole allele,
end to end, against any stretch of the window (semi-global alignment with
affine gaps; +2 match, −3 mismatch, −5 gap open, −2 gap extend). The
alignment is exact dynamic programming, not a heuristic.

Two checks guard the placement:

- **Identity.** Below `--min-identity` the allele is `low_identity` and
  gives no differences: it does not belong at these coordinates.
- **Competing genes.** The allele is also scored against every other gene
  of the list. If another gene fits at least as well, the placement is
  marked `ambiguous_locus`. Its differences are still judged, and the
  output says that the sequence alone does not tell the two genes apart.

**Genes the reference lacks.** Four IGHV genes are absent from GRCh38 and
present in GRCh37. Their alleles are placed and judged on GRCh37, with the
older data that exist for it; the output marks them. See
[Datasets](datasets.md#genes-that-the-reference-genome-lacks).

## 2. From an alignment to differences

Every mismatch, insertion and deletion becomes one record: position,
reference bases, allele bases. Insertions and deletions are shifted as far
left as they go and trimmed to their shortest form (Tan, Abecasis and Kang
2015), and the records of the population files are split and trimmed in the
same way. The same event therefore has the same three values on both sides.

**Matching is exact.** A record of the population data counts for a
difference only if position, reference bases and allele bases are all
identical. A different base at the same position is a different variant and
is never used. An rsID by itself is not evidence either: dbSNP assigns
one to every change anybody submitted.

Reference positions the allele does not cover are unknown. igpop does not
assume the reference base there.

## 3. Reading the population data

### What is read

For the main build (GRCh38) seven collections of people are read, each
through its tabix index, for the gene intervals only:

| Dataset | People | Groups | A difference without a record |
|---|---:|---|---|
| 1000 Genomes, 30x | 2,504 unrelated | 5 continental groups | was looked for: 0 copies |
| gnomAD v4.1 genomes | 76,215 | 9 genetic-ancestry groups | was looked for: 0 copies |
| ALFA (NCBI) | varies by position | 9 populations | says nothing |
| GenomeAsia pilot | 1,739 | 7 groups | was looked for: 0 copies |
| ToMMo 61KJPN | about 61,000 | Japanese | was looked for: 0 copies |
| Westlake BioBank (WBBC) | 4,480 | 4 regions of China | was looked for: 0 copies |
| dbSNP | — | studies relayed by dbSNP | says nothing |

A further file, the phased panel of 1000 Genomes, holds individual
chromosomes and is used only for the haplotype check (step 9).

### One group of people is one stratum

A **stratum** is one group of people of one dataset: the East Asians of
1000 Genomes, the Finns of gnomAD, the Japanese of ToMMo, one study relayed
by dbSNP. For every difference, igpop collects from every stratum two
numbers: how many chromosomes carry it (`AC`) and how many were examined
(`AN`). Everything that follows is computed from this list and from nothing
else. The list itself is written to `populations.tsv`.

How a stratum gets its two numbers depends on what the dataset says:

- **The dataset has a record of exactly this change.** The counts of the
  record are used, group by group.
- **The dataset has no record, and it examined everybody at every
  position** (whole-genome call sets: 1000 Genomes, gnomAD, GenomeAsia,
  ToMMo, WBBC). The change was looked for and not found: zero copies among
  all chromosomes of each group. Absence is evidence, and a large dataset
  without the variant bounds its frequency tightly.
- **The dataset has no record, and it does not cover every position**
  (ALFA, studies relayed by dbSNP). Nothing follows; the stratum is left
  out for this difference.
- **A study relayed by dbSNP** gives a frequency and no counts. igpop
  estimates the size of the study from the smallest frequency it reports
  in the genes read (one copy among all its chromosomes) and turns the
  frequency into a count. Such counts are marked `estimated`. A study is
  used only if it is large enough that seeing nothing would show a
  frequency below `--common`.

Only whole chromosomes count. What a call set reports on an unplaced piece
of sequence is not used as evidence.

### The same people are not counted twice

- The 698 children and other relatives in 1000 Genomes are left out of the
  counts: a child repeats its parents' alleles.
- A study that dbSNP relays is dropped when igpop reads the same dataset
  directly (dbSNP's copies of 1000 Genomes, gnomAD genomes, ALFA, ToMMo).
- **Counts are never added across strata.** gnomAD contains the 1000
  Genomes samples, ExAC is part of the gnomAD exomes, and the datasets
  group people differently. A pooled frequency would count some people
  twice and would mostly reflect who happened to be sampled. igpop instead
  asks each stratum separately and takes the best one (steps 5 and 7).
  Taking a maximum does not count anybody twice.

Which strata a run used, how large each is and why any was left out is
written to `datasets.tsv`. In the IGHV example a difference on GRCh38 can
be looked up in 52 strata: 35 groups of the six datasets with counts and 17
studies relayed by dbSNP.

## 4. From counts to confidence bounds

Three copies among 200 chromosomes and 300 among 20,000 are both 1.5%. The
second is strong evidence that the frequency is above 1%; the first is not.

For every stratum igpop computes the **exact binomial (Clopper-Pearson)
bounds** at 95%, one-sided: the lowest and the highest frequency the counts
are compatible with.

| Copies / chromosomes | Observed | Lowest plausible | Highest plausible |
|---|---:|---:|---:|
| 3 / 200 | 1.5% | 0.41% | 3.8% |
| 300 / 20,000 | 1.5% | 1.36% | 1.65% |
| 0 / 5,008 | 0% | 0% | 0.06% |
| 0 / 152,000 | 0% | 0% | 0.002% |

*Why exact bounds.* The usual approximation (frequency ± 1.64 standard
errors) fails for the counts that matter most here: 0, 1, 2 or 3 copies.
The exact bounds hold for any count and any sample size, at the price of
being slightly cautious.

## 5. Support: one number per difference

The **support** of a difference is the largest *lowest plausible frequency*
over all its strata: the lower bound in the population where that bound is
highest. A difference never seen has support 0.

- It is a lower bound, so a handful of people cannot produce a high value.
- It looks at populations one by one, so a variant frequent in one
  population and rare overall is not diluted.

This is the idea behind gnomAD's *filtering allele frequency* (Whiffin and
colleagues 2017): judge a variant by the highest frequency that some
population credibly shows. gnomAD uses it to argue that a variant is too
common to cause a rare disease; igpop uses it to argue that a variant is
too common to be an error.

**Read it as a measure, not as a guarantee.** In a single population the
bound holds with 95% confidence. The support is the best of about fifty
such bounds, and the best of fifty is right less often than each one. igpop
therefore never interprets the support by itself. What a given support is
worth is measured in the next step.

## 6. The score: what a given support is worth

A database of 76,000 genomes contains an enormous number of rare variants.
In the IGHV genes about half of *all possible* single-base changes have a
record somewhere. "It is in a database" is therefore weak evidence, and how
weak depends on the databases.

The score measures it. igpop computes the support of two groups of changes
through exactly the same data and the same code:

- **real variants**: the differences of alleles accepted as real (by
  default those of the trusted alleles of the run, plus candidates whose
  whole sequence is established in step 9);
- **all possible changes**: every single-base change at every position of
  the genes of the run, three per position. Most of these do not exist.
  They stand for what an error would look like.

The **score** of a support is how many times more often real variants have
that support than possible changes do: a likelihood ratio.

| Score | Meaning |
|---:|---|
| 100 | this support is 100 times more typical of real variants |
| 1 | the data do not point either way |
| 0.01 | this is what an arbitrary change looks like |

*How it is estimated.* Among changes that were seen, a logistic regression
on the logarithm of the support gives the odds of "real" against
"possible" at every support; corrected for the sizes of the two groups,
these odds are the ratio. The curve is not allowed to fall with rising
support. Changes never seen get one score of their own: the share of real
variants never seen divided by the share of possible changes never seen.
[The score](score.md) has the details, the curve and its uncertainty.

*Why this design.*

- **No cut-off per dataset.** A large dataset sees many arbitrary changes
  at low frequency, so low frequencies in it earn a low score. This is
  measured, not set by hand, and it adapts when datasets are added.
- **Overlap and shared errors cancel.** Both groups are measured through
  the same strata. Whatever inflates the support of real variants
  (overlapping samples, the best of fifty bounds, an artefact of short
  reads) inflates the support of arbitrary changes in the same way.
- **It errs on the cautious side.** Some of the "possible changes" are real
  variants nobody has named yet. They make arbitrary changes look better
  supported than errors really are, which lowers every score.

*What the score is not.* It is not the probability that a variant is real.
That would need to know how many of the submitted alleles are wrong, which
nobody knows. It says how the evidence shifts the balance, whatever the
starting point, in the same way as the likelihood ratios used to weigh
evidence for and against a disease variant.

An insertion or deletion borrows the curve of single-base changes, and
only when that speaks for it (score above 1): random insertions and
deletions are found in population data far less often than random base
changes, so the same support is at least as strong evidence for them.

## 7. The class: is the difference demonstrably common?

Independently of the score, igpop makes one formal statement that needs no
calibration and that any geneticist will recognise:

> In at least one population, this variant has a frequency of 1% or more.

**The test.** In every stratum with at least two copies (a single copy can
be one miscalled genotype), an exact binomial test asks: *if the true
frequency were just below 1%, how likely would it be to see this many
copies or more?* A small probability (p-value) means the counts are hard
to explain unless the frequency really is 1% or more. The limit is
`--common`.

**Why the p-value is multiplied (Bonferroni).** A difference is tested in
up to 52 strata and the best result is kept. With 52 attempts, a "1 in 20"
coincidence turns up in one of them most of the time. The smallest p-value
is therefore multiplied by the number of strata of the run, and the
difference is `common` only if the product is 0.05 or less. In words: the
chance of calling a rare variant common, in *any* of the populations, stays
below 5%.

Three properties of this correction matter here:

- **It holds when strata overlap.** Bonferroni needs no assumption about
  how the tests depend on each other. With overlapping datasets it is
  stricter than necessary, never too lenient.
- **The number is fixed before any variant is looked at**: every group
  that has chromosomes in the genes read, and every relayed study in use.
  It does not depend on where a particular variant happened to be found,
  which would be choosing the correction after seeing the data.
- **It is per difference.** igpop judges each difference of each allele on
  its own evidence and does not correct across differences. Among very
  many truly rare variants, up to 5% could be called common by chance;
  the score, not this label, is what separates the doubtful from the sound.

The test and the bounds of step 4 are the same calculation read two ways:
the p-value is 0.05 or less exactly when the lowest plausible frequency is
1% or more. The correction asks for a stricter bound (99.9% instead of 95%
with 52 strata).

**dbSNP's COMMON flag** is honoured as a second-hand version of the same
statement, for the one allele of the record it belongs to, unless the data
igpop reads directly show the variant to be rare in every population.

**The other classes.** Seen, but not shown to be common: `uncommon`.
Looked for by at least one complete dataset and seen nowhere: `absent`. No
population data for its position: `not_assessed`.

**What the class is for.** It gives a plain reason that does not depend on
a calibration ("common in gnomAD Africans, p = 3e-12"), and it is a safety
net: a variant shown to be common is `supported` whatever the score curve
of the run says. With the default limits it rarely decides by itself. In
the IGHV example all 319 differences classed `common` also score 2 or
more, and a further 99 are supported by their score alone.

## 8. The verdict of a difference

| Class | Score | Verdict |
|---|---|---|
| `common` | any | `supported` |
| `uncommon` or `absent` | 2 or more | `supported` |
| `uncommon` or `absent` | between 0.1 and 2 | `inconclusive` |
| `uncommon` or `absent` | 0.1 or less | `needs_confirmation` |
| `uncommon`, no score | — | `inconclusive` |
| `absent`, no score | — | `needs_confirmation` |
| `not_assessed` | — | `not_assessed` |

The two limits are `--supported-score` and `--confirm-score`, and the web
app has a slider for the first. A limit of 2 means: this much support is at
least twice as typical of real variants. That is the strength the
guidelines for clinical variants call *supporting* evidence (a likelihood
ratio of about 2; Tavtigian and colleagues 2018). It is a deliberately low
bar for a single difference, because an allele must clear it at *every*
difference. [The score](score.md#what-is-measured-and-what-is-a-convention) shows
what each limit costs.

For a difference that is not supported, the output also says what the data
show (`never_observed`, `single_observation`, `rare_in_every_population`,
`common_not_demonstrated`): see [Why an allele did not pass](filters.md).

## 9. The haplotype check

Frequencies show that each difference exists. They do not show that the
differences occur *together*, on one chromosome, as the allele claims.

The phased panel of 1000 Genomes gives, for 2,598 unrelated people, the two
copies of each chromosome separately. igpop reconstructs the sequence of
each of these 5,196 chromosomes over the gene and counts those identical
to the allele over its whole length.

Finding a sequence there is not proof either, because a chromosome that
differs from a real allele by one common variant will also be found. The
check therefore has its own yardstick, built like the score: igpop takes
the real chromosomes of the panel, changes one base at a time, and counts
how often such an arbitrary variant of a real sequence is found on as many
chromosomes. An allele is **established as a haplotype** if its exact
sequence is on at least two chromosomes and fewer than 5% of the arbitrary
variants are found that often. See [The haplotype check](haplotypes.md).

## 10. The verdict of an allele

1. An allele identical to the reference over its length:
   `matches_reference`.
2. By frequency, **an allele is as good as its weakest difference**: its
   verdict is the lowest verdict among its differences. One doubtful base
   makes the allele doubtful, however many others are sound.
3. **An allele established as a haplotype is `supported`** whatever its
   single differences say. Its whole sequence was found in people, which
   is stronger than any statement about its parts.
4. **An allele whose differences all pass, but whose sequence a complete
   lookup finds on no phased chromosome, is `inconclusive`**
   (`sequence_not_in_panel`): each difference exists, the combination was
   not seen.

Every verdict comes with its reason in words and with the numbers behind
it, down to the counts in each population.

## Site quality: reported, not used

Whether a call set can be trusted at a position is a different question
from how frequent a variant is there. igpop reads what each dataset says
about its own sites (filters, flags for segmental duplications and
low-complexity sequence, missing genotypes, Hardy-Weinberg and Mendelian
errors) and shows it next to every difference and in `site_quality.tsv`.

It does not use it in a verdict, for a reason specific to these genes: in
IGHV most positions lie in regions that the datasets flag as duplicated,
and their general-purpose filters reject many variants that long-read
assemblies and repertoire studies confirm. Filtering on them would discard
most of the evidence; ignoring them silently would overstate it. igpop
shows both. A supported allele whose evidence rests on flagged sites
deserves a second look, above all when it has no haplotype support. See
[Site quality filters](site-qc.md).

## The statistical methods at a glance

| Where | Method | Why this one |
|---|---|---|
| Placement | Exact semi-global alignment with affine gaps | The place is given; only the best alignment within it is needed, and it should not depend on heuristics. |
| Same variant on both sides | Left-alignment and trimming | One event must have one representation, or lookups miss. |
| Frequency of one group | Exact binomial (Clopper-Pearson) bounds, one-sided 95% | Correct for counts of 0 to 3, where approximations fail. |
| Several groups | The best group, never a pooled count | Datasets overlap and sample populations unevenly. |
| Worth of a support | Likelihood ratio from a logistic regression, real variants against all possible changes | Measures how informative the databases are in these genes instead of assuming it. |
| "Common" | Exact binomial test against a frequency below 1%, in each group | The standard statement about a variant, with a known error rate. |
| Best of many groups | Bonferroni, by the number of groups of the run | Valid however the groups overlap; fixed in advance. |
| Whole sequence | Exact lookup on phased chromosomes; share of arbitrary one-base variants found as often | Finding a sequence means little unless arbitrary sequences are not found. |
| Allele from differences | The weakest difference | An allele is wrong if one base is wrong. |

## How this relates to existing methods

igpop does not replicate one existing method. Its parts are established;
their combination, applied to judging an inferred allele, is new.

- **IgSNPer**, its predecessor, looked up each difference in dbSNP and
  counted the uncommon ones, with one frequency per study and a fixed
  cut-off. igpop keeps the idea and replaces the count by evidence per
  population, with sample size taken into account.
- **gnomAD's filtering allele frequency** is the model for the support:
  the highest frequency a population credibly shows.
- **Likelihood ratios against a background of possible variants** are how
  CADD scores variants (observed against simulated changes) and how GATK's
  variant recalibration works (known sites against all calls). igpop
  applies the idea to one feature, the support, and fits it per run.
- **Evidence strength as a likelihood ratio** is the framework of the
  clinical guidelines (ACMG/AMP in its Bayesian form). igpop's limit of 2
  corresponds to their weakest level.
- **Alleles from population data** have been derived in the other
  direction: Lym1K and pmIG assembled IG alleles from the phased
  haplotypes of 1000 Genomes and kept those found on several chromosomes.
  The objections raised against them (short reads placed on the wrong
  copy of a duplicated gene; genes missing from the reference) apply to
  igpop's evidence as well, which is why site quality is shown and why the
  allele itself must come from elsewhere: igpop confirms, it does not
  discover.

What is specific to igpop: the user fixes the gene, so no aligner decides
between near-identical genes; every available population is consulted
separately and absence is counted as evidence; the worth of the evidence
is calibrated on the genes themselves in every run; and the whole sequence
is checked on phased chromosomes against its own background.

What the method cannot do is on the page [Limitations](limitations.md).

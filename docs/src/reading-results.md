# Reading the results

## The idea in three sentences

An allele is described by its differences from the reference genome. A
difference that many people carry is an established variant; a difference
nobody carries among several hundred thousand sequenced genomes behaves like
a sequencing or PCR error. igpop measures, for every difference, which of
the two it resembles, and an allele is as believable as its least
believable difference. Second, independently: if chromosomes with exactly
the allele's sequence exist in phased genomes, the allele is real whatever
its single differences look like.

## Verdicts

Every allele (in `alleles.tsv` and in the report) and every variant (in
`variants.tsv`) gets one of these:

| Verdict | What it means | What to do |
|---|---|---|
| `supported` | Every difference is an established population variant (`supported_by` = `frequency`), or the exact sequence of the allele is established on phased chromosomes (`haplotype`), or both. | Nothing. Population data back this allele. |
| `inconclusive` | At least one difference is seen in people, but so rarely that real rare variants and errors look alike. | Decide with your own evidence: seen in several unrelated individuals, in both chains of a haplotype analysis, in genomic DNA? |
| `needs_confirmation` | At least one difference is absent, or nearly absent, from all population data: it looks like an arbitrary base change. | Do not accept on repertoire data alone. Confirm independently (genomic sequencing, more individuals), then it may well be a real private allele. |
| `not_assessed` | No population data could be read for at least one difference. | See the reason and the notes; usually the gene lies outside the chromosomes of the reference, where no population data are read. |
| `matches_reference` | No difference from the reference genome. | Nothing to assess. |
| `not_placed` | The allele could not be aligned to a locus. | See the reason: gene missing from the coordinates, or the sequence is too different from its gene. |

The verdict is about **population support, not truth**. `supported` does not
prove an allele exists as one piece (each difference is checked on its own),
and `needs_confirmation` does not prove an error: a real allele private to
one family has no population support either.

## Two more columns to read with the verdict

- **`filters`** says in fixed keywords why an allele did not pass, e.g.
  `rare_in_every_population:1,low_score:1`. See
  [Why an allele did not pass](filters.md).
- **`supported_by`** says which evidence made an allele `supported`:
  `frequency`, `haplotype` or `frequency+haplotype`.
- **`hap_status`** says whether the whole combination of variants was
  found on one chromosome in phased genomes. It often has no answer;
  when it finds the allele on enough chromosomes, the allele is
  `supported`. See [The haplotype check](haplotypes.md).

## The `reason` and `notes` columns

`reason` names the variants behind the verdict:

```
1 of 2 variant(s) without population support: rs368426916(G>A:212,106184983)
```

Read a variant as `rsID(reference base > allele base : position in your
allele, position in the genome)`. A `.` instead of an rsID means dbSNP has
no entry for this change. Bases are written on the reference strand, so for
IGH genes (minus strand) they are the complement of what your FASTA shows;
the position in the allele counts along your sequence as supplied.

`notes` lists things that weaken or complicate the reading:

| Note | Meaning |
|---|---|
| fits another locus at least as well | The allele is as similar to another gene (a duplicate such as IGHV1-69 / IGHV1-69D). Its variants may belong to the other locus. |
| locus shared with ... | Two genes have identical coordinates in the table, so population data cannot tell them apart. |
| shorter than the gene | The allele is truncated. The missing part was not checked. |
| variant(s) at sites failing 1000 Genomes site QC | The site is hard to genotype in 1000 Genomes (see [Site quality filters](site-qc.md)). Frequencies there are less reliable in either direction. |

## From verdict to numbers

For each variant, `variants.tsv` holds what the verdict is made of:

| Column | In plain words |
|---|---|
| `class` | `common`: shown to reach 1% in at least one population. `uncommon`: seen, but not shown to reach 1%. `absent`: never seen. |
| `support` | The lowest plausible frequency in the best population: the lower end of the 95% confidence range, so small studies cannot inflate it. |
| `score` | A likelihood ratio: how many times more often this level of support occurs among variants of real alleles than among arbitrary base changes. Around 70 = typical of real variants; 1 = no information either way; 0.01 = typical of arbitrary changes. |
| `verdict` | The verdict of the variant by frequency: `supported` if `class` is `common` or the score is at least 2; `needs_confirmation` if the score is at most 0.1; otherwise `inconclusive`. |

A score of 2 means that such support is twice as typical of real variants
as of arbitrary changes; a score of 0.1 means ten to one for an arbitrary
change. Both limits can be changed
(`--supported-score`, `--confirm-score`).

## Why not simply count uncommon variants?

`alleles.tsv` keeps the counts (`common`, `uncommon`, `absent`), and they
agree with the verdict in most cases. They are not enough on their own
because "uncommon" lumps together two very different things: a variant
carried by 1 in 300 people, and a base change never observed in anyone. The
score separates them, and the verdict follows the score. See
[The score](score.md).

## Worked examples

**`IGHV1-18*01`: supported.** One difference, rs72695948, frequent in every
population. Score 72.

**`IGHV1-18*01_S2706`: needs_confirmation.** Two differences. One is the
frequent rs72695948 above; the other, rs368426916, has an rsID but was found in only a handful of
genomes: its lowest plausible frequency is a few in a million
(`rare_in_every_population,low_score`). An rsID is a label, not evidence:
dbSNP records every submitted change. The allele needs independent
confirmation of that one base.

**`IGHV1-2*04_S3434`: inconclusive.** A single difference, seen in people
at roughly 0.1 to 1%. Real rare variants and errors that happen to hit a
real variant both occur at this level (score 2), and its site is not in the
phased panel, so the haplotype check has no answer.

**`IGHV1-24*01_S7069`: supported by haplotype.** By frequency its one
difference is equally weak. But chromosomes with exactly this sequence
exist in 9 people of the phased panel, which fewer than 1 in 300 arbitrary
changes achieve.

**`IGHV1-24*01_S1381`: supported by frequency and haplotype.** Its exact
sequence is on 72 phased chromosomes of 64 people, 70 of them African: a
real, population-specific allele.

**`IGHV1-18*03`: a trusted allele that is inconclusive.** One of its two
differences (rs373407210) is rare in every population, and the haplotype
check cannot test the allele. Trusted alleles are judged like all others,
which is how database entries with a doubtful or very rare base come to
light.

## Probabilities

A score is a likelihood ratio, not a probability, and not odds either:
odds that a variant is real = your odds beforehand x score. igpop does not
turn it into a probability, because that needs the share of your
candidates that are real before looking at population data, which only
you can judge. [The score](score.md) has a table of how much it matters.

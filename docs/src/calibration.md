# Automatic calibration

The [score](score.md) compares a variant with two reference groups: variants
of real alleles (the yardstick), and all possible single-base changes in the
same genes.
Building that comparison is called calibration. You never have to run it:
igpop does it in every run and picks the source in this order.

| Order | Source | Used when |
|---|---|---|
| 1 | A curve you pass with `--calibration` (a `score_model.tsv`) | You ask for it. |
| 2 | **This run** | The yardstick alleles of the run carry at least 30 distinct single-base variants, and these are better supported than arbitrary changes. |
| 3 | **The built-in curve** | Otherwise. |

The report (section 5) and the last column of
`<prefix>.score_levels.tsv` (`score_source`: `run`, `built-in` or
`file`) state which one was used. The curve in use is always written to
`<prefix>.score_model.tsv`.

## Calibrating on the run itself

This is the preferred case, and what happens for a full V set. It adapts the
scale to the genes, the datasets and the `--common` level of this run.

The yardstick has two parts:

- **Named trusted alleles** (below). These are mostly common alleles.
- **Alleles confirmed as haplotypes**, found automatically: alleles with
  two or more differences whose exact sequence is established on phased
  chromosomes ([The haplotype check](haplotypes.md)). They include rare
  real alleles, which named databases under-represent; without them the
  score would learn that "real" means "frequent". Column `yardstick` of
  `alleles.tsv` says `named` or `haplotype`. One caveat: every variant in
  the phased panel has at least two copies in 1000 Genomes, so these
  alleles say nothing about the lowest supports.

Requirements:

- **Trusted alleles.** By default all alleles without `_S` in the name
  (`--candidate-marker`), or a list (`--trusted`). They should be alleles
  you would defend: database alleles confirmed in genomic DNA or in many
  individuals. They do not have to be perfect; a few doubtful ones barely
  move the scale, and they will show up in the report as trusted alleles
  that are not `supported`.
- **Enough of them.** 30 distinct single-base variants is the minimum;
  below 100 the report warns that the curve is not pinned down well.

Put candidates and trusted alleles in the same FASTA and the same run.

## The built-in curve

Small sets (D and J genes, or a handful of alleles) cannot calibrate
themselves. They are scored with a curve computed once from 566 human IGHV
alleles (147 named trusted, 45 confirmed as haplotypes) on GRCh38 with
the seven datasets of that build: 213 real variants against 38,316 possible
changes. A change never seen scores 0.005; for a change that was seen,
log(score) = 4.88 + 1.00 x − 0.26 x² with x = log10(support). The two
groups, counted per level of support:

| Support | Arbitrary changes | Variants of yardstick alleles | Score across the level |
|---|---:|---:|---:|
| none | 18,172 | 0 | 0.005 |
| below 0.001% | 6,066 | 0 | 0.005 |
| 0.001 to 0.01% | 7,551 | 0 | 0.005 to 0.04 |
| 0.01 to 0.1% | 4,152 | 7 | 0.04 to 0.64 |
| 0.1 to 1% | 1,658 | 14 | 0.64 to 6.3 |
| 1 to 10% | 347 | 45 | 6.3 to 37 |
| 10% and more | 370 | 147 | 37 to 132 |

It is appropriate for human IGH with the default datasets. For another
locus, or other datasets, the report warns; the better route is a run of
that locus large enough to calibrate itself, whose curve you then reuse:

```bash
igpop variants --alleles trbv.fasta --coords trb.tsv -o trbv
igpop variants --alleles trbj.fasta --coords trb.tsv --calibration trbv.score_model.tsv -o trbj
```

## Is my calibration good?

Open section 5 and 6 of the report and check three things.

1. **Shape.** Arbitrary changes concentrate at "none", trusted variants at
   1% and more, and the score rises from far below 1 to far above 10, with
   intervals that are narrow where it crosses the limits.
2. **Trusted alleles.** Most trusted alleles with differences are
   `supported`.
3. **No calibration warning**, or one you understand.

If the shape is flat (scores near 1 everywhere), population data do not
separate your trusted alleles from arbitrary changes, and no verdict should
be relied on. The usual causes are coordinates on the wrong genome build and
a trusted list that contains the candidates.

## Changing the verdict limits

By frequency, `supported` starts at a score of 2, `needs_confirmation` at
0.1 and below. (The haplotype route has its own limit, `--coincidence`.)
To be stricter, raise the first or raise the second:

```bash
igpop variants --alleles my.fasta --supported-score 10 --confirm-score 0.3 -o strict
```

Section 5 of the report says at which support the score crosses each limit
(0.29% and 0.021% for the defaults in the IGHV example), so you can see
what a limit means in terms of frequency.

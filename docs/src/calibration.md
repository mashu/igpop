# Automatic calibration

The [score](score.md) compares a variant with two reference groups: variants
of trusted alleles, and all possible single-base changes in the same genes.
Building that comparison is called calibration. You never have to run it:
igpop does it in every run and picks the source in this order.

| Order | Source | Used when |
|---|---|---|
| 1 | A table you pass with `--calibration` | You ask for it. |
| 2 | **This run** | The trusted alleles of the run carry at least 30 distinct single-base variants. |
| 3 | **The built-in table** | Otherwise. |

The report (section 5) and the last column of
`<prefix>.score_calibration.tsv` (`score_source`: `run`, `built-in` or
`file`) state which one was used.

## Calibrating on the run itself

This is the preferred case, and what happens for a full V set. It adapts the
scale to the genes, the datasets and the `--common` level of this run.

Requirements:

- **Trusted alleles.** By default all alleles without `_S` in the name
  (`--candidate-marker`), or a list (`--trusted`). They should be alleles
  you would defend: database alleles confirmed in genomic DNA or in many
  individuals. They do not have to be perfect; a few doubtful ones barely
  move the scale, and they will show up in the report as trusted alleles
  that are not `supported`.
- **Enough of them.** 30 distinct single-base variants is the minimum;
  below 100 the report warns that the scale is coarse.

Put candidates and trusted alleles in the same FASTA and the same run.

## The built-in table

Small sets (D and J genes, or a handful of alleles) cannot calibrate
themselves. They are scored with a table computed once from 559 human IGHV
alleles on GRCh38 with 1000 Genomes 30x, gnomAD v4.1 genomes and the
datasets relayed by dbSNP:

| Support | Arbitrary changes | Trusted variants | Score |
|---|---:|---:|---:|
| none | 18,812 | 0 | 0.006 |
| below 0.001% | 5,488 | 0 | 0.017 |
| 0.001 to 0.01% | 7,809 | 0 | 0.017 |
| 0.01 to 0.1% | 3,899 | 1 | 0.084 |
| 0.1 to 1% | 1,614 | 8 | 1.156 |
| 1 to 10% | 332 | 26 | 17.502 |
| 10% and more | 362 | 136 | 82.689 |

It is appropriate for human IGH with the default datasets. For another
locus, or other datasets, the report warns; the better route is a run of
that locus large enough to calibrate itself, whose table you then reuse:

```bash
igpop variants --alleles trbv.fasta --coords trb.tsv -o trbv
igpop variants --alleles trbj.fasta --coords trb.tsv --calibration trbv.score_calibration.tsv -o trbj
```

## Is my calibration good?

Open section 5 and 6 of the report and check three things.

1. **Shape.** Arbitrary changes concentrate at "none", trusted variants at
   1% and more, and the score rises from far below 1 to far above 10.
2. **Trusted alleles.** Most trusted alleles with differences are
   `supported`.
3. **No calibration warning**, or one you understand.

If the shape is flat (scores near 1 everywhere), population data do not
separate your trusted alleles from arbitrary changes, and no verdict should
be relied on. The usual causes are coordinates on the wrong genome build and
a trusted list that contains the candidates.

## Changing the verdict limits

`supported` starts at a score of 10, `needs_confirmation` at 0.1 and below.
To be stricter, raise the first or raise the second:

```bash
igpop variants --alleles my.fasta --supported-score 20 --confirm-score 0.3 -o strict
```

Because the score has seven levels, moving a limit only matters when it
crosses one of the seven values in the calibration table.

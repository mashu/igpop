# The run report and its checks

`<prefix>.report.md` is written by every run. It has seven short numbered
sections and a haplotype section.
Sections 1 and 2 are the result; 3 to 6 tell you whether the result can be
relied on. Nothing has to be requested: calibration and all checks are
automatic.

## 1. Result at a glance

Number of alleles per verdict, split into trusted alleles and candidates.

## 2. Alleles to look at

Three lists: `needs_confirmation`, `inconclusive` and `not_assessed`, each
with the variants responsible, the outcome of the haplotype check, and any
notes. This is the working list.

## Haplotype check

How many alleles had their whole combination of variants found on one
chromosome of the phased panel, the calibration table (how often an
arbitrary change is found as often), and from how many chromosomes an
allele counts as supported by its haplotype. Explained in
[The haplotype check](haplotypes.md).

## 3. Was every allele placed on its gene?

| Placement | Meaning | Action |
|---|---|---|
| `ok` | Aligned to the locus of its gene. | None. |
| `ambiguous_locus` | Fits another locus at least as well. | Normal for duplicated genes. Read the verdict with the note in mind. |
| `low_identity` | Less than 90% identical to its locus. | Wrong gene name, wrong coordinates, or not this gene. Look at `alignments.txt`. |
| `no_coordinates` | No coordinates for its gene on this build. | Add the gene to the coordinates table, or run the other build. |

Many `low_identity` alleles at once mean the coordinates do not match the
genome build.

## 4. Population data used

Which datasets were read, from where, and how many chromosomes each
contributes. If dbSNP could not be reached, rsIDs are missing and the report
says so; run again.

## 5. How the score was calibrated

Where the scores come from (see [Automatic calibration](calibration.md)) and
the table itself:

| Column | Meaning |
|---|---|
| population support | Level of the `support` value, from "none" to "at least 10%". |
| arbitrary changes | How many of all possible single-base changes in your genes have this support. |
| trusted variants | How many variants of trusted alleles have this support. |
| score | Ratio of the two shares. |
| verdict | What a variant at this level is called. |
| probability real if half ... / if 9 in 10 ... | The score turned into a probability under two expectations. |

A healthy table has most arbitrary changes at the top ("none") and most
trusted variants at the bottom (1% and more), with scores rising from well
below 1 to well above 10.

## 6. Checks on this run

**Trusted alleles, judged like candidates.** Trusted alleles get verdicts
too. Most of those that differ from the reference should be `supported`.
In the IGHV example, 110 of 117 are. If only a minority is, something is
wrong upstream: coordinates on the wrong build, a strand problem, or a
trusted set that is not trustworthy. Individual trusted alleles that are
not `supported` are worth a look in their own right.

**Variants found only in candidate alleles, by support.** Where the
novelty in your candidates sits. Errors pile up at "none"; candidates that
reuse established variants in new combinations sit at 1% and more; the
levels in between are the inconclusive zone.

**Warnings.** Each states what happened and what to do:

| Warning | Meaning |
|---|---|
| scores come from the built-in curve | Your run was too small to calibrate itself. Expected for D, J and small sets. |
| the built-in curve was made from IGH genes, but this run holds other genes | The scale may not transfer to another locus. Prefer a larger run of that locus. |
| the calibration rests on only N trusted variants | It worked, but is coarse. Verdicts near the limits can move. |
| only N of M trusted variants reach 1% | The trusted alleles do not look like established variation. Check build, coordinates and the trusted set before trusting any verdict. |
| alleles were not placed cleanly | See section 3. |
| no dataset with allele counts was read | No score; `supported` only through dbSNP's COMMON flag. |
| the phased panel could not be read | No haplotype check in this run. Run again. |
| dbSNP could not be read | rsIDs and COMMON flags are missing. Run again. |

The warnings are also printed on screen at the end of the run.

## 7. Files of this run

A list of the output files; details in [Output files](outputs.md).

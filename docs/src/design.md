# How the code is organised

A run is four steps. Each step is one module, takes the result of the
previous one, and has no other way of influencing the outcome.

| Step | Module | In | Out |
|---|---|---|---|
| 1. Place | `locus` | allele sequences, gene coordinates, reference | `Placed`: where each allele aligned, its `Status`, its variants |
| 2. Read evidence | `evidence` | the loci that need data | `Evidence`: 1000 Genomes sites, gnomAD sites, dbSNP records, the phased `Panel` |
| 3. Assess | `assess` | `Placed` + `Evidence` | `ScoreModel` (calibration), then `VariantResult` and `AlleleResult` |
| 4. Write | `tables`, `report` | the results | the eight output files |

`main.rs` contains the command line and the four calls, nothing else.

## One interface between data and decisions

Everything that is decided about a variant is computed from its **strata**:
the list of `(population, allele count, chromosomes)` that
`Evidence::strata` returns. A stratum is a superpopulation of 1000 Genomes,
an ancestry group of gnomAD, or a dataset relayed by dbSNP.

- class: an exact test in each stratum, best one, Bonferroni-corrected;
- support: the largest lower confidence bound over the strata;
- score: the calibration table at the level of the support;
- filters: number of copies and bounds in the strata.

Adding a data source means adding strata in `evidence.rs`. Nothing in
`assess.rs` changes, and the calibration adapts by itself because arbitrary
changes and trusted variants are measured through the same function.

## Decisions are in one place

`assess.rs` opens with the rules in words; the code below follows them in
the same order. The states are enums (`Status`, `Class`, `Verdict`,
`Filter`, `ScoreSource`, `Kind`), so a new state cannot be half-handled:
the compiler lists every place that must deal with it.

## Tables are declared, not formatted

A table in `tables.rs` is a list of columns, each a name together with the
function that fills it. The header and the rows are produced from the same
list, and a block of columns that fills the wrong number of cells stops the
program instead of writing a shifted table. To add a column, add one
`.col("name", |row| ...)` line and describe it in `outputs.md`.

## Supporting modules

| Module | Content |
|---|---|
| `align`, `variant` | Semi-global affine alignment; left-aligned, trimmed variants. |
| `freq` | Allele counts, exact (Clopper-Pearson) bounds, exact binomial test, heterogeneity test; parsers for the 1000 Genomes and gnomAD records. |
| `score` | Support levels, the likelihood-ratio table, the built-in table. |
| `hap` | Phased panel as bit sets per variant; counting carrier chromosomes of founders. |
| `io`, `tabix`, `fetch` | FASTA and coordinates; remote range reads; tabix index; region cache. |
| `render` | Alignments as text. |

## Checking a change

```bash
cargo fmt --check && cargo clippy --all-targets -- -D warnings && cargo test
```

Unit tests cover the alignment, variant normalisation, statistics, the
calibration, the haplotype counting, the verdict rules and the table
mechanism; they do not use the network. For a change that should not alter
results, run the same alleles before and after with the cache in place and
compare the output files byte by byte.

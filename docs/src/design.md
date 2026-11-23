# How the code is organised

A run is four steps. Each step is one module, takes the result of the
previous one, and has no other way of influencing the outcome.

| Step | Module | In | Out |
|---|---|---|---|
| 1. Place | `locus` | allele sequences, gene coordinates, reference | `Placed`: where each allele aligned, its `Status`, its variants |
| 2. Read evidence | `evidence` | the loci that need data | `Evidence`: 1000 Genomes sites, gnomAD sites, dbSNP records, the phased `Panel` |
| 3. Assess | `assess` | `Placed` + `Evidence` | `ScoreModel` (calibration), then `VariantResult` and `AlleleResult` |
| 4. Write | `tables`, `report` | the results | the eight output files |

The four calls are in `igpop::run` (`lib.rs`), which takes its inputs as
text and returns its outputs as text. Two programs use it:

| Program | Where | Adds |
|---|---|---|
| Command line | `src/main.rs` (feature `cli`) | options, reading and writing files, blocking HTTP (`ureq`), threads (`rayon`) |
| Web app | `web/` (crate `igpop-web`, WebAssembly) and `web/static/` (the page) | HTTP through the browser, progress messages to the page |

The library itself has no network or thread dependency: the transport is
registered through `net::set_transport`, progress goes through
`log::set_sink`, and `par::map` runs in parallel only with the `cli`
feature.

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
| `io`, `tabix`, `fetch`, `net` | FASTA and coordinates; tabix index; region cache; the pluggable transport. |
| `log`, `par` | Progress messages; parallel map. |
| `render` | Alignments as text. |

## Where the time goes

Comparing every allele with every locus (to notice an allele that fits
another gene as well) is most of the computation. Two things keep it fast
without changing any result:

- `align::score_only` computes a row of the alignment in passes over plain
  slices, which the compiler turns into vector instructions; only the gap
  running along the row is carried from cell to cell.
- `locus::place_all` compares the first allele of each gene with every
  locus. A later allele of the gene differs from it by a few edits, and
  each edit can change a score by a bounded amount, so loci that cannot
  reach the best competitor found so far are skipped. The bound is exact;
  a unit test and whole-run comparisons check that results are identical.

The command line aligns on all cores and downloads over three connections.
The web app has one core and one download at a time; there the downloads
dominate.

## Checking a change

```bash
cargo fmt --check && cargo clippy --all-targets -- -D warnings && cargo test
```

Unit tests cover the alignment, variant normalisation, statistics, the
calibration, the haplotype counting, the verdict rules and the table
mechanism; they do not use the network. For a change that should not alter
results, run the same alleles before and after with the cache in place and
compare the output files byte by byte.

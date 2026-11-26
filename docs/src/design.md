# How the code is organised

## A run is four steps

Each step is one module, takes the result of the previous one, and has no
other way of influencing the outcome.

| Step | Module | In | Out |
|---|---|---|---|
| 1. Place | `locus` | allele sequences, gene coordinates, reference | `Placed`: where each allele aligned, its `Status`, its variants |
| 2. Read evidence | `evidence` | the loci that need data, the datasets of the build | `Evidence`: per dataset the records read, dbSNP, the phased `Panel` |
| 3. Assess | `assess` | `Placed` + `Evidence` | `ScoreModel` (calibration), then `VariantResult` and `AlleleResult` |
| 4. Write | `tables`, `report` | the results | the output files |

`run` makes the four calls. It takes its inputs as text and returns its
outputs as text; `plan` works out beforehand what a run will read.

## Layers

A module uses only modules of the layers above it; there are no circular
dependencies.

| Layer | Modules | What they know |
|---|---|---|
| Arithmetic | `stats`, `score`, `align` | Counts and bounds, the score curve, alignment. Nothing about files, genes or datasets. |
| Vocabulary | `variant`, `config` | What a variant is; the limits and choices of a run. |
| Transport | `net`, `io`, `tabix`, `fetch` | FASTA and coordinate files; reading a region of an indexed file, from disk or over HTTP. |
| Datasets | `population` (`layout`, `dbsnp`), `hap`, `catalogue` | What a dataset is, how each file format is read, which datasets exist for which genome build. |
| The four steps | `locus`, `evidence`, `assess` (`haplotype`, `calibration`, `variant`, `allele`), `tables`, `report` | The method. |
| A run | `plan`, `run` | What a run will read; the order of the steps. |
| Programs | `src/main.rs`, `web/` | Options, files, threads, the browser. |

Two programs use the library:

| Program | Where | Adds |
|---|---|---|
| Command line | `src/main.rs` (feature `cli`) | options, reading and writing files, blocking HTTP (`ureq`), threads (`rayon`) |
| Web app | `web/` (crate `igpop-web`, WebAssembly) and `web/static/` (the page) | HTTP through the browser, progress messages to the page |

The library itself has no network or thread dependency: the transport is
registered through `net::set_transport`, progress goes through
`log::set_sink`, and `par::map` runs in parallel only with the `cli`
feature.

A run can also be taken in two halves (`Session::start`, then
`Session::finish`), with `Session::downloads` in between naming the files
and regions the second half will read. The web app uses this to download
them ahead, three at a time per server (`web/src/ahead.rs`): the reader of
a region (`tabix::Reader`) does no I/O itself, it says which bytes it needs
and is fed them, so it asks for exactly the same bytes when downloading
ahead as during the analysis, which then finds them already there.

## One description of a dataset

Every collection of people is a `population::Dataset`: a name, its groups
of people (`Group`: code, title, continental group), where its file is,
how the file is laid out (`Layout`), and three facts that decide how it
is used:

- `complete`: it examined everybody at every position, so a variant
  without a record was looked for and not seen;
- `replaces`: the studies relayed by dbSNP that are copies of it and must
  be dropped when it is read;
- `required`: whether a run can go on without it.

`catalogue` lists the datasets of each genome build. Nothing else in the
program names a dataset: `evidence`, `assess`, the tables, the report and
the page work with whatever the catalogue holds.

**To add a dataset**, add one entry to `catalogue.rs`. If its file has
columns like an existing one (per-group `AC_x`/`AN_x` INFO fields, sample
columns, or frequencies), choose that `Layout` and name the fields;
otherwise add a variant to `Layout` in `population/layout.rs` with its
parser and a unit test. The new groups then appear by themselves in the
support, the test for `common` and its correction, the calibration, the
three evidence tables, the report and the page. Run the example before
and after: the score curve will move, because arbitrary changes are
measured through the new dataset too.

## One interface between data and decisions

Everything that is decided about a variant is computed from its **strata**:
the list of `(stratum, copies, chromosomes)` that `Evidence::strata`
returns. A stratum is a group of a dataset or a study relayed by dbSNP.

- support: the largest lower confidence bound over the strata
  (`evidence::support`);
- score: the calibration curve at the support;
- class: an exact test in each stratum, the best one, corrected for
  `Evidence::strata_in_use()`;
- filters: number of copies and bounds in the strata.

The calibration measures arbitrary changes and real variants through the
same `Evidence::strata`, so anything that changes the evidence changes
both sides of the ratio.

## Decisions are in one place

`assess/mod.rs` opens with the rules in words, and its four files follow
them in order: `haplotype`, `calibration`, `variant`, `allele`. The states
are enums (`Status`, `Class`, `Verdict`, `Filter`, `ScoreSource`, `Kind`,
`Ancestry`), so a new state cannot be half-handled: the compiler lists
every place that must deal with it.

## Tables are declared, not formatted

A table in `tables/` is a list of columns, each a name together with the
function that fills it. The header and the rows are produced from the same
list, and a block of columns that fills the wrong number of cells stops the
program instead of writing a shifted table. To add a column, add one
`.col("name", |row| ...)` line and describe it in `outputs.md`.

Evidence is written long, one row per variant and group
(`populations.tsv`), not as a column per dataset: a new dataset adds rows,
never columns.

## The page

`web/static/index.html` loads `js/main.js`; everything else is an ES
module imported from there. Modules import only from layers above their
own, and a check for circular imports is part of reviewing a change.

| Layer | Modules |
|---|---|
| Tools | `dom`, `tsv`, `words` (every term and help text shown to the user) |
| State | `state` (the run on display, and a small event channel), `engine` (messages to and from the worker) |
| Data | `populations` (the three evidence tables), `score` (the curve), `names`, `mapping` |
| Parts of the results | `variants`, `haplotype`, `alignment`, `limit`, `pages` |
| Pages | `files`, `check`, `name-editor`, `gene-list`, `progress`, `results` |
| Start | `main` |
| Apart | `worker` (runs the WebAssembly engine), `store` (downloads kept in the browser) |

The styles are split the same way, in `css/`: `tokens` (colours and
sizes, light and dark), `base` (what every page shares), `frame`, then one
file per page or part (`start`, `check`, `progress`, `results`, `allele`,
`alignment`, `pages`). A colour that carries meaning is set once: a
verdict class (`v-supported`) sets `--v`, a continental group
(`a-east_asian`) sets `--ancestry`, and marks use those.

The page computes nothing about a verdict except when the user moves the
limit, and then with the curve the engine wrote (`score_model.tsv`).

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
The web app has one core; there the downloads dominate.

## Checking a change

```bash
cargo fmt --check && cargo clippy --all-targets -- -D warnings && cargo test
(cd web && cargo check)
for f in web/static/js/*.js; do node --check "$f"; done
```

Unit tests cover the alignment, variant normalisation, statistics (among
them that the test for `common` and the lower bound agree), the dataset
layouts, the calibration, the haplotype counting, the verdict rules and
the table mechanism; they do not use the network. For a change that should
not alter results, run the same alleles before and after with the cache in
place and compare the output files.

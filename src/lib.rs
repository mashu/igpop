//! igpop: population evidence for germline alleles at given loci.
//!
//! One run ([`run`]) is four steps:
//!
//! 1. `locus`:    place every allele on the locus of its gene, call variants;
//! 2. `evidence`: read population data for those loci (only the bytes needed);
//! 3. `assess`:   class, score, verdict and filters per variant and allele;
//! 4. `tables`, `report`: the result tables and the report, as text.
//!
//! The modules form layers, each using only the ones before it:
//!
//! - numbers and sequences: `stats` (counts, exact bounds and tests),
//!   `score` (the calibration curve), `align`, `variant` (normalised
//!   variants);
//! - files: `net`, `io`, `tabix`, `fetch` (remote files, read by region),
//!   `log`, `par`;
//! - the data: `population` (datasets, how their records are laid out,
//!   dbSNP), `hap` (the phased panel), `catalogue` (what is read by
//!   default for each genome build);
//! - the steps: `locus`, `evidence`, `assess`, `tables`, `report`, `render`;
//! - a run: `config` (what was asked), `plan` (what will be read, known
//!   before any download), `run` (the steps in order).
//!
//! The library reads its inputs as text and returns its outputs as text,
//! so the same code serves the command-line program (`main.rs`) and the
//! web app (`web/`).

pub mod align;
pub mod assess;
pub mod catalogue;
pub mod config;
pub mod evidence;
pub mod fetch;
pub mod hap;
pub mod io;
pub mod locus;
pub mod log;
pub mod net;
pub mod par;
pub mod plan;
pub mod population;
pub mod render;
pub mod report;
pub mod run;
pub mod score;
pub mod stats;
pub mod tabix;
pub mod tables;
pub mod variant;

pub use config::{Choice, Config, Limits, Text};
pub use plan::{preview, AllelePreview, GenePreview, PlannedSource, Preview, BY_SEQUENCE};
pub use run::{run, Output, Session};

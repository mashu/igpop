//! What a run is asked to do: its inputs as text, its limits, and which
//! data it may read.

use anyhow::{bail, Result};
use std::path::PathBuf;

/// The limits a verdict depends on.
#[derive(Clone)]
pub struct Limits {
    /// Frequency a `common` variant must demonstrably reach.
    pub common: f64,
    /// Score from which a variant is `supported`.
    pub supported: f64,
    /// Score up to which a variant is `needs_confirmation`.
    pub confirm: f64,
    /// A haplotype supports an allele if fewer than this share of
    /// arbitrary changes are found as often.
    pub coincidence: f64,
}

/// Where one data source is read from.
#[derive(Clone, Default)]
pub enum Choice {
    /// The default for the build, if there is one.
    #[default]
    Default,
    /// This path or URL.
    Custom(String),
    /// Not used.
    Off,
}

impl Choice {
    /// `main`: for the build the run was asked for. A file given by the
    /// user belongs to that build; any other build uses its default.
    pub(crate) fn resolve(&self, default: Option<&str>, main: bool) -> Option<String> {
        match self {
            Choice::Custom(spec) if main => Some(spec.clone()),
            Choice::Default | Choice::Custom(_) => default.map(str::to_string),
            Choice::Off => None,
        }
    }
}

/// A named text: the content of an input file and what to call it in
/// messages.
#[derive(Clone)]
pub struct Text {
    pub name: String,
    pub content: String,
}

/// Everything a run needs. `Config::new` gives the defaults.
pub struct Config {
    /// Allele FASTA.
    pub alleles: Text,
    /// Coordinates table; `None` = the shipped human IGH table.
    pub coords: Option<Text>,
    /// Genome build, in any spelling `io::normalise_build` accepts.
    pub build: String,
    /// Judge a gene that has coordinates only on another build with the
    /// data of that build.
    pub other_builds: bool,
    /// `allele<TAB>gene` table; `None` = the shipped one, if the shipped
    /// coordinates are used.
    pub aliases: Option<Text>,
    /// Place an allele whose gene is not in the coordinates table on the
    /// locus its sequence fits best, instead of leaving it unplaced. The
    /// alias table asks the same for single alleles or genes with the
    /// gene `?` (`BY_SEQUENCE`).
    pub place_unknown: bool,
    /// Trusted allele names, one per line; `None` = by `candidate_marker`.
    pub trusted: Option<String>,
    /// Alleles whose name contains this are candidates, the others trusted.
    pub candidate_marker: String,
    /// A `score_model.tsv` of an earlier run.
    pub calibration: Option<Text>,
    pub limits: Limits,
    /// Reference FASTA; `None` = the default for the build.
    pub reference: Option<String>,
    pub kgp: Choice,
    pub gnomad: Choice,
    pub dbsnp: Choice,
    /// Read the further datasets of the build (ALFA, GenomeAsia, ToMMo,
    /// WBBC).
    pub more_datasets: bool,
    pub phased: Choice,
    pub pedigree: Option<String>,
    /// dbSNP datasets left out by name.
    pub dbsnp_exclude: Vec<String>,
    /// Where downloaded regions are kept; `None` = no cache.
    pub cache_dir: Option<PathBuf>,
    /// Reference bases added on each side of a locus.
    pub flank: u64,
    pub min_identity: f64,
    pub paralog_check: bool,
    /// Output prefix, used in the report to name the other files.
    pub prefix: String,
}

impl Config {
    pub fn new(alleles: Text) -> Self {
        Config {
            alleles,
            coords: None,
            build: "38".to_string(),
            other_builds: true,
            aliases: None,
            place_unknown: false,
            trusted: None,
            candidate_marker: "_S".to_string(),
            calibration: None,
            limits: Limits {
                common: 0.01,
                supported: 2.0,
                confirm: 0.1,
                coincidence: 0.05,
            },
            reference: None,
            kgp: Choice::Default,
            gnomad: Choice::Default,
            dbsnp: Choice::Default,
            more_datasets: true,
            phased: Choice::Default,
            pedigree: None,
            dbsnp_exclude: vec!["SGDP_PRJ".to_string()],
            cache_dir: None,
            flank: 100,
            min_identity: 0.9,
            paralog_check: true,
            prefix: "igpop".to_string(),
        }
    }

    pub(crate) fn check(&self) -> Result<()> {
        if self.limits.confirm >= self.limits.supported {
            bail!("the confirm score must be below the supported score");
        }
        Ok(())
    }
}

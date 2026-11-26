//! The score of a run: which calibration is used, and where it came from.

use super::Yardstick;
use crate::evidence::{self, Evidence};
use crate::locus::{Placed, Window};
use crate::score::{self, Calibration};
use crate::variant::{Kind, Lookup};
use anyhow::Result;
use std::collections::HashSet;

/// Where the scores of a run come from.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ScoreSource {
    /// Calibrated on this run.
    Run,
    /// The built-in table: too few trusted variants in this run.
    BuiltIn,
    /// A table passed with `--calibration`.
    File,
    /// No scores: no source with allele counts was read.
    None,
}

impl ScoreSource {
    pub fn as_str(self) -> &'static str {
        match self {
            ScoreSource::Run => "run",
            ScoreSource::BuiltIn => "built-in",
            ScoreSource::File => "file",
            ScoreSource::None => "none",
        }
    }
}

/// The scores in use and where they came from.
pub struct ScoreModel {
    pub calibration: Option<Calibration>,
    pub source: ScoreSource,
    /// The source in words, for the report.
    pub description: String,
}

impl ScoreModel {
    /// Tabulate this run and choose the scores: a file (its text and name)
    /// if given; otherwise
    /// this run; otherwise, with too few trusted variants, the built-in table.
    pub fn new(ev: &Evidence, windows: &[Window], yardstick: &[(&Placed, Yardstick)], file: Option<(&str, &str)>) -> Result<Self> {
        // The support of every possible change and of every real variant.
        let (mut arbitrary, mut trusted) = (Vec::new(), Vec::new());
        if ev.counted() {
            let mut positions: HashSet<(&str, u64)> = HashSet::new();
            for w in windows.iter().filter(|w| ev.covers(w)) {
                for pos in w.coord.start.max(w.start)..=w.coord.end.min(w.end()) {
                    if !positions.insert((w.coord.contig.as_str(), pos)) {
                        continue;
                    }
                    let reference = w.seq[(pos - w.start) as usize];
                    for alt in b"ACGT".iter().filter(|alt| **alt != reference && reference != b'N') {
                        let lookup = (w.coord.contig.clone(), (pos, vec![reference], vec![*alt]));
                        arbitrary.push(evidence::support(&ev.strata(&lookup)));
                    }
                }
            }
            let mut variants: HashSet<Lookup> = HashSet::new();
            for (p, _) in yardstick {
                for c in p.variants.iter().filter(|c| c.kind == Kind::Snv) {
                    let lookup = p.lookup(c);
                    if variants.insert(lookup.clone()) {
                        trusted.push(evidence::support(&ev.strata(&lookup)));
                    }
                }
            }
        }
        let (na, nt) = (arbitrary.len(), trusted.len());
        let (calibration, source, description) = match (file, ev.counted()) {
            (Some((text, name)), _) => (
                Some(score::parse(text, name)?),
                ScoreSource::File,
                format!("the file `{name}` (`--calibration`)"),
            ),
            (None, false) => (None, ScoreSource::None, String::new()),
            (None, true) => match score::calibrate(&arbitrary, &trusted) {
                Some(c) => (
                    Some(c),
                    ScoreSource::Run,
                    format!(
                        "this run: {nt} single-base variants of {} named trusted alleles and {} alleles confirmed as haplotypes, against {na} possible single-base changes in the same genes",
                        yardstick.iter().filter(|(_, y)| *y == Yardstick::Named).count(),
                        yardstick.iter().filter(|(_, y)| *y == Yardstick::Haplotype).count(),
                    ),
                ),
                None => (
                    Some(score::builtin()),
                    ScoreSource::BuiltIn,
                    format!(
                        "the built-in curve, because this run has only {nt} single-base variants of trusted alleles (at least {} are needed, and they must be better supported than arbitrary changes). \
                         It was computed from 566 human IGHV alleles on GRCh38 with 1000 Genomes 30x, gnomAD v4.1 genomes, dbSNP, ALFA, GenomeAsia, ToMMo 61KJPN and WBBC",
                        score::MIN_TRUSTED
                    ),
                ),
            },
        };
        Ok(ScoreModel {
            calibration,
            source,
            description,
        })
    }
}

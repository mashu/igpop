//! One allele: the verdicts of its variants and its haplotype check
//! combined into one verdict, with the reason in words.

use super::{assess_variant, AlleleHaplotype, Filter, HaplotypeSupport, Judge, VariantResult, Verdict};
use crate::hap;
use crate::locus::{Placed, Status};
use crate::stats;

/// Why an allele's variants calibrate the score.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Yardstick {
    /// Named as trusted and cleanly placed.
    Named,
    /// Two or more differences whose combination is established as a
    /// haplotype: real variation found by the data, rare alleles included.
    Haplotype,
    No,
}

impl Yardstick {
    pub fn as_str(self) -> &'static str {
        match self {
            Yardstick::Named => "named",
            Yardstick::Haplotype => "haplotype",
            Yardstick::No => ".",
        }
    }
}

/// Everything concluded about one allele.
/// Filter of an allele whose variants pass but whose sequence the phased
/// panel does not have.
pub const SEQUENCE_NOT_IN_PANEL: &str = "sequence_not_in_panel";

pub struct AlleleResult<'a> {
    pub placed: &'a Placed,
    /// The genome build the allele was placed and judged on.
    pub build: &'a str,
    /// Named as trusted (by `--trusted` or the candidate marker).
    pub trusted: bool,
    /// Whether, and why, its variants calibrate the score.
    pub yardstick: Yardstick,
    pub variants: Vec<VariantResult>,
    /// The verdict from the frequencies of its variants alone.
    pub by_frequency: Verdict,
    /// The verdict: `by_frequency`, raised to `supported` if the allele's
    /// exact sequence is established as a haplotype.
    pub verdict: Verdict,
    pub reason: String,
    pub notes: Vec<String>,
    pub haplotype: Option<hap::Assessment>,
    /// Share of arbitrary changes found on at least as many chromosomes
    /// as this allele's exact sequence.
    pub haplotype_coincidence: Option<f64>,
    /// Set if the haplotype check supports the allele.
    pub haplotype_support: Option<HaplotypeSupport>,
}

impl AlleleResult<'_> {
    /// Its variants calibrate the score.
    pub fn calibrates(&self) -> bool {
        self.yardstick != Yardstick::No
    }

    /// Which evidence makes the allele `supported`.
    pub fn supported_by(&self) -> &'static str {
        if self.verdict != Verdict::Supported {
            return ".";
        }
        match (self.by_frequency == Verdict::Supported, self.haplotype_support.is_some()) {
            (true, true) => "frequency+haplotype",
            (true, false) => "frequency",
            (false, true) => "haplotype",
            (false, false) => ".",
        }
    }

    pub fn ids(&self, wanted: impl Fn(&VariantResult) -> bool) -> Vec<String> {
        self.variants.iter().filter(|v| wanted(v)).map(|v| v.id.clone()).collect()
    }

    /// The lowest-scoring variant: an allele is only as well supported as
    /// its weakest variant.
    pub fn weakest(&self) -> Option<(f64, &str)> {
        let mut weakest: Option<(f64, &str)> = None;
        for v in &self.variants {
            if let Some(x) = v.score {
                if weakest.is_none_or(|(w, _)| x < w) {
                    weakest = Some((x, &v.id));
                }
            }
        }
        weakest
    }

    /// Why the allele did not pass: the filters of its variants with counts.
    pub fn filters(&self) -> String {
        if matches!(self.placed.status, Status::NoCoordinates | Status::LowIdentity) {
            return self.placed.status.as_str().to_string();
        }
        let parts: Vec<String> = Filter::ALL
            .iter()
            .filter_map(|f| {
                let k = self.variants.iter().filter(|v| v.filters.contains(f)).count();
                (k > 0).then(|| format!("{}:{k}", f.as_str()))
            })
            .collect();
        if self.contradicted() {
            return SEQUENCE_NOT_IN_PANEL.to_string();
        }
        if parts.is_empty() {
            ".".to_string()
        } else {
            parts.join(",")
        }
    }

    /// Supported by the frequencies of its variants, yet the phased panel,
    /// which has every one of them, has the exact sequence on no chromosome.
    pub fn contradicted(&self) -> bool {
        self.by_frequency == Verdict::Supported && self.verdict == Verdict::Inconclusive
    }

    pub fn hap_status(&self) -> &'static str {
        self.haplotype.as_ref().map_or("not_checked", |h| h.status(self.variants.len()))
    }
}

/// `other_build`: a build on which the allele's gene does have coordinates.
pub fn assess_allele<'a>(
    judge: &Judge,
    placed: &'a Placed,
    build: &'a str,
    trusted: bool,
    yardstick: Yardstick,
    haplotype: AlleleHaplotype,
    other_build: Option<&str>,
) -> AlleleResult<'a> {
    let variants: Vec<VariantResult> = placed.variants.iter().map(|c| assess_variant(judge, placed, c)).collect();
    let total = variants.len();
    let with = |verdict: Verdict| -> Vec<&str> { variants.iter().filter(|v| v.verdict == verdict).map(|v| v.id.as_str()).collect() };
    let weakest = |verdict: Verdict, what: &str| {
        let ids = with(verdict);
        (verdict, format!("{} of {total} variant(s) {what}: {}", ids.len(), ids.join(";")))
    };
    // By frequency, the verdict of an allele is that of its weakest variant.
    let (by_frequency, frequency_reason) = match placed.status {
        Status::NoCoordinates => (
            Verdict::NotPlaced,
            match other_build {
                Some(other) => format!("its gene has coordinates on {other} only: run again with --build {other}"),
                None => "its gene is missing from the coordinates table".to_string(),
            },
        ),
        Status::LowIdentity => (
            Verdict::NotPlaced,
            "too different from the locus of its gene (below --min-identity)".to_string(),
        ),
        _ if total == 0 => (Verdict::MatchesReference, "no difference from the reference genome".to_string()),
        _ if !with(Verdict::NeedsConfirmation).is_empty() => weakest(Verdict::NeedsConfirmation, "without population support"),
        _ if !with(Verdict::Inconclusive).is_empty() => weakest(Verdict::Inconclusive, "with weak population support"),
        _ if !with(Verdict::NotAssessed).is_empty() => weakest(Verdict::NotAssessed, "without population data"),
        _ => (Verdict::Supported, format!("all variants ({total}) have population support")),
    };
    // An allele whose exact sequence is established as a haplotype is
    // supported whatever the frequencies of its single variants.
    let undecided = matches!(
        by_frequency,
        Verdict::Inconclusive | Verdict::NeedsConfirmation | Verdict::NotAssessed
    );
    // The reverse also holds: frequencies show that each variant exists, not
    // that they occur together. A complete lookup that finds the sequence on
    // no phased chromosome leaves the combination unsupported.
    let absent = haplotype
        .assessment
        .as_ref()
        .and_then(|a| a.exact.as_ref())
        .map(|exact| stats::total(exact))
        .filter(|exact| exact.ac == 0);
    let (verdict, reason) = match (&haplotype.support, absent) {
        (_, Some(exact)) if by_frequency == Verdict::Supported => (
            Verdict::Inconclusive,
            format!(
                "every variant has population support, but none of the {} phased chromosomes has this exact sequence",
                exact.an
            ),
        ),
        (Some(h), _) if undecided => (
            Verdict::Supported,
            format!(
                "exact sequence on {} chromosomes of {} people in the phased panel; by frequency alone: {frequency_reason}",
                h.chromosomes, h.donors
            ),
        ),
        _ => (by_frequency, frequency_reason),
    };

    let mut notes: Vec<String> = Vec::new();
    if placed.status == Status::AmbiguousLocus {
        notes.push("fits another locus at least as well (see alignments.tsv)".to_string());
    }
    if let Some(p) = &placed.placement {
        if !p.shared.is_empty() {
            notes.push(format!("locus shared with {}", p.shared.join(",")));
        }
        if p.missing.0 + p.missing.1 > 0 {
            notes.push(format!(
                "shorter than the gene: {} base(s) missing at 5', {} at 3'",
                p.missing.0, p.missing.1
            ));
        }
    }
    let flagged = variants.iter().filter(|v| !v.site_qc.is_empty()).count();
    if flagged > 0 {
        notes.push(format!(
            "{flagged} difference(s) at positions where the genotypes of a dataset are less reliable (column site_quality says why)"
        ));
    }
    AlleleResult {
        placed,
        build,
        trusted,
        yardstick,
        variants,
        by_frequency,
        verdict,
        reason,
        notes,
        haplotype: haplotype.assessment,
        haplotype_coincidence: haplotype.coincidence,
        haplotype_support: haplotype.support,
    }
}

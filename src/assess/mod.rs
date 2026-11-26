//! Step 3: turn the population evidence of each variant into a class, a
//! score and a verdict, and the verdicts of its variants into the verdict
//! of an allele.
//!
//! The rules, in the order they are applied:
//!
//! 1. **Haplotype check** (`haplotype`). Is the allele's exact sequence on
//!    phased chromosomes, more often than an arbitrary change of a real
//!    chromosome would be? It runs first because it also says which
//!    candidates are real enough to calibrate the score.
//! 2. **Calibration** (`calibration`). The score curve of the run, from
//!    the variants of alleles accepted as real against every possible
//!    single-base change in the same genes.
//! 3. **Variant** (`variant`). Class: `common` if an exact binomial test
//!    shows a frequency of `Limits::common` or more in some stratum, after
//!    correction for the number of strata; else `uncommon` if seen, `absent`
//!    if looked for and not seen, `not_assessed` without data. Score: the
//!    curve at the variant's support. Verdict: `variant_verdict`.
//! 4. **Allele** (`allele`). By frequency an allele is as good as its
//!    weakest variant; an exact sequence established as a haplotype is
//!    supported whatever its single variants say, and one whose variants
//!    all pass but whose sequence a complete lookup finds on no chromosome
//!    is inconclusive.

mod allele;
mod calibration;
mod haplotype;
mod variant;

pub use allele::{assess_allele, AlleleResult, Yardstick, SEQUENCE_NOT_IN_PANEL};
pub use calibration::{ScoreModel, ScoreSource};
pub use haplotype::{AlleleHaplotype, HaplotypeCheck, HaplotypeSupport};
pub use variant::{assess_variant, VariantResult};

use crate::config::Limits;
use crate::evidence::Evidence;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Class {
    Common,
    Uncommon,
    Absent,
    NotAssessed,
}

impl Class {
    pub fn as_str(self) -> &'static str {
        match self {
            Class::Common => "common",
            Class::Uncommon => "uncommon",
            Class::Absent => "absent",
            Class::NotAssessed => "not_assessed",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    Supported,
    Inconclusive,
    NeedsConfirmation,
    NotAssessed,
    MatchesReference,
    NotPlaced,
}

impl Verdict {
    pub fn as_str(self) -> &'static str {
        match self {
            Verdict::Supported => "supported",
            Verdict::Inconclusive => "inconclusive",
            Verdict::NeedsConfirmation => "needs_confirmation",
            Verdict::NotAssessed => "not_assessed",
            Verdict::MatchesReference => "matches_reference",
            Verdict::NotPlaced => "not_placed",
        }
    }

    /// An allele with this verdict is listed in the final table.
    pub fn passes(self) -> bool {
        matches!(self, Verdict::Supported | Verdict::MatchesReference)
    }
}

/// Why a variant is not `supported`. The first five say what the data
/// show (exactly one applies), the last three how that scored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    NoPopulationData,
    NeverObserved,
    SingleObservation,
    RareInEveryPopulation,
    CommonNotDemonstrated,
    LowScore,
    IntermediateScore,
    NotScored,
}

impl Filter {
    pub const ALL: [Filter; 8] = [
        Filter::NoPopulationData,
        Filter::NeverObserved,
        Filter::SingleObservation,
        Filter::RareInEveryPopulation,
        Filter::CommonNotDemonstrated,
        Filter::LowScore,
        Filter::IntermediateScore,
        Filter::NotScored,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            Filter::NoPopulationData => "no_population_data",
            Filter::NeverObserved => "never_observed",
            Filter::SingleObservation => "single_observation",
            Filter::RareInEveryPopulation => "rare_in_every_population",
            Filter::CommonNotDemonstrated => "common_not_demonstrated",
            Filter::LowScore => "low_score",
            Filter::IntermediateScore => "intermediate_score",
            Filter::NotScored => "not_scored",
        }
    }
}

/// The verdict of one variant from its class and, where it has one, its score.
pub fn variant_verdict(class: Class, score: Option<f64>, limits: &Limits) -> Verdict {
    match (class, score) {
        (Class::Common, _) => Verdict::Supported,
        (Class::NotAssessed, _) => Verdict::NotAssessed,
        (_, Some(x)) if x >= limits.supported => Verdict::Supported,
        (_, Some(x)) if x <= limits.confirm => Verdict::NeedsConfirmation,
        (_, Some(_)) => Verdict::Inconclusive,
        (Class::Absent, None) => Verdict::NeedsConfirmation,
        (_, None) => Verdict::Inconclusive,
    }
}

/// What a variant or an allele is judged with.
pub struct Judge<'a> {
    pub ev: &'a Evidence,
    pub model: &'a ScoreModel,
    pub limits: &'a Limits,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verdict_rules() {
        let limits = Limits {
            common: 0.01,
            supported: 10.0,
            confirm: 0.1,
            coincidence: 0.05,
        };
        let v = |class, score| variant_verdict(class, score, &limits).as_str();
        assert_eq!(v(Class::Common, Some(0.04)), "supported");
        assert_eq!(v(Class::Uncommon, Some(23.0)), "supported");
        assert_eq!(v(Class::Uncommon, Some(1.5)), "inconclusive");
        assert_eq!(v(Class::Absent, Some(0.04)), "needs_confirmation");
        assert_eq!(v(Class::Absent, None), "needs_confirmation");
        assert_eq!(v(Class::Uncommon, None), "inconclusive");
        assert_eq!(v(Class::NotAssessed, None), "not_assessed");
    }
}

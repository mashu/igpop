//! One variant: its class, its score and its verdict.

use super::{variant_verdict, Class, Filter, Judge, Verdict};
use crate::evidence::{self, Stratum};
use crate::locus::{Called, Placed};
use crate::stats::{informative, PopCount};
use crate::variant::{Kind, Lookup};

/// True if every stratum large enough to tell bounds the frequency below
/// `common`, and there is at least one such stratum.
fn rare_everywhere(strata: &[(Stratum, PopCount)], common: f64) -> bool {
    let mut large = strata.iter().filter(|(_, c)| informative(c.an, common)).peekable();
    large.peek().is_some() && large.all(|(_, c)| c.upper95() < common)
}

/// Everything concluded about one variant of one allele.
pub struct VariantResult {
    pub lookup: Lookup,
    pub rsid: String,
    /// `rsid(REF>ALT:position in the allele,position in the genome)`.
    pub id: String,
    /// Sources in which the allele was observed.
    pub observed: Vec<&'static str>,
    pub support: f64,
    /// The stratum that gives the support, as `dataset:population`.
    pub support_from: Option<String>,
    pub class: Class,
    /// The stratum that made the variant common, or `dbSNP:COMMON`.
    pub basis: Option<String>,
    /// Corrected p-value of the best stratum, if one was tested.
    pub class_p: Option<f64>,
    pub score: Option<f64>,
    /// The score of an insertion or deletion is read from the curve of
    /// single-base changes.
    pub score_borrowed: bool,
    /// A 95% interval for the score: how well the calibration knows it.
    pub score_interval: Option<(f64, f64)>,
    /// Share of all possible single-base changes in the genes of the run
    /// that are supported at least as well as this variant.
    pub chance: Option<f64>,
    pub verdict: Verdict,
    /// What the data show about a variant that is not common; the first
    /// of `filters` when the variant is not supported.
    pub shows: Option<Filter>,
    pub filters: Vec<Filter>,
    /// The site criteria the variant's position fails, per dataset.
    pub site_qc: Vec<(&'static str, Vec<String>)>,
}

pub fn assess_variant(judge: &Judge, placed: &Placed, called: &Called) -> VariantResult {
    let (ev, model, limits) = (judge.ev, judge.model, judge.limits);
    let lookup = placed.lookup(called);
    let rsid = ev.dbsnp.rsids.get(&lookup).cloned().unwrap_or_else(|| ".".into());
    let (pos, reference, alternate) = (
        &called.key.0,
        String::from_utf8_lossy(&called.key.1),
        String::from_utf8_lossy(&called.key.2),
    );
    let id = format!("{rsid}({reference}>{alternate}:{},{pos})", called.allele_pos);
    let strata = ev.strata(&lookup);
    let support = evidence::support(&strata);
    // The population whose counts guarantee that frequency.
    let support_from = strata
        .iter()
        .filter(|(_, c)| c.ac > 0)
        .max_by(|a, b| a.1.lower95().total_cmp(&b.1.lower95()))
        .map(|(stratum, _)| ev.stratum_name(*stratum));

    let observed = ev.observed(&lookup);
    // Class: the best stratum among those where the allele was seen at
    // least twice (one copy can be a genotyping error), corrected for the
    // number of strata.
    let best = strata
        .iter()
        .filter(|(_, c)| c.ac >= 2)
        .map(|(name, c)| (c.tail_p(limits.common), c.lower95(), name))
        // Smallest p-value; among equals, the highest guaranteed frequency.
        .min_by(|x, y| x.0.total_cmp(&y.0).then(y.1.total_cmp(&x.1)))
        .map(|(p, _, stratum)| (p, *stratum));
    let corrected = best.map(|(p, stratum)| ((p * strata.len() as f64).min(1.0), ev.stratum_name(stratum)));
    let (class, basis) = match &corrected {
        Some((p, name)) if *p <= 0.05 => (Class::Common, Some(name.clone())),
        // dbSNP's flag is second-hand. It does not count where the data read
        // directly show the variant to be rare in every population.
        _ if ev.dbsnp.common.contains(&lookup) && !rare_everywhere(&strata, limits.common) => {
            (Class::Common, Some("dbSNP:COMMON".to_string()))
        }
        _ if !observed.is_empty() => (Class::Uncommon, None),
        _ if ev.reads(&lookup.0) => (Class::Absent, None),
        _ => (Class::NotAssessed, None),
    };

    // The curve is made from single-base changes. An insertion or deletion
    // borrows it, but only as a lower bound: of all possible single-base
    // insertions and deletions far fewer are found in population data than
    // of all possible base changes, so the same support is at least as
    // strong evidence for them. A lower bound can speak for a variant and
    // cannot speak against it: below 1 it says nothing, and the variant
    // stays without a score.
    let scored = match &model.calibration {
        Some(cal) if ev.reads(&lookup.0) && (called.kind == Kind::Snv || cal.score(support) > 1.0) => {
            Some((cal.score(support), cal.interval(support), cal.chance(support)))
        }
        _ => None,
    };
    let score = scored.map(|(x, _, _)| x);
    let verdict = variant_verdict(class, score, limits);

    // What the data show about a variant that is not common, whatever
    // its verdict: the verdict depends on a limit, this does not.
    let most = strata.iter().map(|(_, c)| c.ac).max().unwrap_or(0);
    let shows = match class {
        Class::Common => None,
        Class::NotAssessed => Some(Filter::NoPopulationData),
        Class::Absent => Some(Filter::NeverObserved),
        _ if most <= 1 => Some(Filter::SingleObservation),
        _ if rare_everywhere(&strata, limits.common) => Some(Filter::RareInEveryPopulation),
        _ => Some(Filter::CommonNotDemonstrated),
    };
    let mut filters = Vec::new();
    if verdict != Verdict::Supported {
        filters.extend(shows);
        match score {
            Some(x) if x <= limits.confirm => filters.push(Filter::LowScore),
            Some(_) => filters.push(Filter::IntermediateScore),
            None if class != Class::NotAssessed => filters.push(Filter::NotScored),
            None => {}
        }
    }
    let site_qc = ev.failures(&lookup);
    VariantResult {
        lookup,
        rsid,
        id,
        observed,
        support,
        support_from,
        class,
        basis,
        class_p: corrected.map(|(p, _)| p),
        score,
        score_borrowed: scored.is_some() && called.kind != Kind::Snv,
        score_interval: scored.map(|(_, interval, _)| interval),
        chance: scored.map(|(_, _, chance)| chance),
        verdict,
        shows,
        filters,
        site_qc,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rare_needs_a_stratum_large_enough_to_tell() {
        let s = |ac, an| (Stratum::Relayed(0), PopCount::new(ac, an));
        assert!(rare_everywhere(&[s(3, 100_000), s(0, 50)], 0.01));
        assert!(!rare_everywhere(&[s(3, 100_000), s(40, 5_000)], 0.01));
        assert!(!rare_everywhere(&[s(2, 50)], 0.01));
    }
}

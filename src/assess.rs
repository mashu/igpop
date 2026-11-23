//! Step 3: turn the population evidence of each variant into a class, a
//! score and a verdict, and the verdicts of its variants into the verdict
//! of an allele.
//!
//! The rules, in the order they are applied to a variant:
//!
//! 1. **Class** (a test). `common` if, in at least one stratum where the
//!    allele was seen twice or more, an exact binomial test rejects
//!    "frequency below `common`" at 5% after Bonferroni correction for the
//!    number of strata; or if dbSNP flags exactly this allele as COMMON.
//!    Otherwise `uncommon` if seen anywhere, `absent` if looked for and not
//!    seen, `not_assessed` if no data was consulted.
//! 2. **Score** (a likelihood ratio; single-base changes only). The level
//!    of the variant's support, looked up in the calibration table.
//! 3. **Verdict.** `supported` if common or score >= `supported`;
//!    `needs_confirmation` if score <= `confirm` (unscored: if absent);
//!    otherwise `inconclusive`.
//! 4. **Filters.** For a variant that is not supported: what the data show
//!    and how that scored, as fixed keywords.

use crate::evidence::{self, informative, Evidence, Lookup, Stratum};
use crate::freq::PopCount;
use crate::hap;
use crate::locus::{Called, Placed, Status, Window};
use crate::score::{self, Calibration, BINS};
use crate::variant::Kind;
use anyhow::Result;
use std::collections::HashSet;

pub struct Limits {
    /// Frequency a `common` variant must demonstrably reach.
    pub common: f64,
    /// Score from which a variant is `supported`.
    pub supported: f64,
    /// Score up to which a variant is `needs_confirmation`.
    pub confirm: f64,
}

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

/// True if every stratum large enough to tell bounds the frequency below
/// `common`, and there is at least one such stratum.
fn rare_everywhere(strata: &[(Stratum, PopCount)], common: f64) -> bool {
    let mut large = strata.iter().filter(|(_, c)| informative(c.an, common)).peekable();
    large.peek().is_some() && large.all(|(_, c)| c.upper95() < common)
}

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
    /// Support levels of all possible single-base changes in the loci read,
    /// and of the single-base variants of trusted alleles, in this run.
    pub arbitrary: [u64; BINS],
    pub trusted: [u64; BINS],
}

impl ScoreModel {
    /// Tabulate this run and choose the scores: a file (its text and name)
    /// if given; otherwise
    /// this run; otherwise, with too few trusted variants, the built-in table.
    pub fn new(ev: &Evidence, windows: &[Window], trusted_alleles: &[&Placed], file: Option<(&str, &str)>) -> Result<Self> {
        let (mut arbitrary, mut trusted) = ([0u64; BINS], [0u64; BINS]);
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
                        arbitrary[score::bin(evidence::support(&ev.strata(&lookup)))] += 1;
                    }
                }
            }
            let mut variants: HashSet<Lookup> = HashSet::new();
            for p in trusted_alleles {
                for c in p.variants.iter().filter(|c| c.kind == Kind::Snv) {
                    let lookup = p.lookup(c);
                    if variants.insert(lookup.clone()) {
                        trusted[score::bin(evidence::support(&ev.strata(&lookup)))] += 1;
                    }
                }
            }
        }
        let (na, nt): (u64, u64) = (arbitrary.iter().sum(), trusted.iter().sum());
        let (calibration, source, description) = match (file, ev.counted()) {
            (Some((text, name)), _) => (
                Some(score::parse(text, name)?),
                ScoreSource::File,
                format!("the file `{name}` (`--calibration`)"),
            ),
            (None, false) => (None, ScoreSource::None, String::new()),
            (None, true) => match score::calibrate(arbitrary, trusted) {
                Some(c) => (
                    Some(c),
                    ScoreSource::Run,
                    format!("this run: {nt} single-base variants of trusted alleles against {na} possible single-base changes in the same genes"),
                ),
                None => (
                    Some(score::builtin()),
                    ScoreSource::BuiltIn,
                    format!(
                        "the built-in table, because this run has only {nt} single-base variants of trusted alleles (at least {} are needed). \
                         It was computed from 559 human IGHV alleles on GRCh38 with 1000 Genomes 30x, gnomAD v4.1 genomes and dbSNP",
                        score::MIN_TRUSTED
                    ),
                ),
            },
        };
        Ok(ScoreModel {
            calibration,
            source,
            description,
            arbitrary,
            trusted,
        })
    }

    /// The counts behind the scores in use, if known.
    pub fn counts(&self) -> Option<([u64; BINS], [u64; BINS])> {
        match self.source {
            ScoreSource::Run => Some((self.arbitrary, self.trusted)),
            ScoreSource::BuiltIn => Some((score::BUILTIN_ARBITRARY, score::BUILTIN_TRUSTED)),
            ScoreSource::File | ScoreSource::None => None,
        }
    }
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
    pub class: Class,
    /// The stratum that made the variant common, or `dbSNP:COMMON`.
    pub basis: Option<String>,
    /// Corrected p-value of the best stratum, if one was tested.
    pub class_p: Option<f64>,
    pub score: Option<f64>,
    pub verdict: Verdict,
    pub filters: Vec<Filter>,
    /// 1000 Genomes site criteria the variant's site fails.
    pub site_qc: Vec<String>,
}

pub fn assess_variant(ev: &Evidence, placed: &Placed, called: &Called, model: &ScoreModel, limits: &Limits) -> VariantResult {
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

    let kgp_site = ev.kgp.get(&lookup);
    let observed: Vec<&'static str> = [
        (evidence::KGP, kgp_site.is_some_and(|s| s.total().ac > 0)),
        // AC = 0 after gnomAD's genotype filtering is not an observation.
        (evidence::GNOMAD, ev.gnomad.get(&lookup).is_some_and(|s| s.total.ac > 0)),
        ("dbSNP", !ev.dbsnp.relayed(&lookup).is_empty()),
    ]
    .iter()
    .filter(|(_, seen)| *seen)
    .map(|(name, _)| *name)
    .collect();

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
        _ if ev.dbsnp.common.contains(&lookup) => (Class::Common, Some("dbSNP:COMMON".to_string())),
        _ if !observed.is_empty() => (Class::Uncommon, None),
        _ if ev.counted() => (Class::Absent, None),
        _ => (Class::NotAssessed, None),
    };

    let score = match &model.calibration {
        Some(cal) if called.kind == Kind::Snv && ev.counted() => Some(cal.score[score::bin(support)]),
        _ => None,
    };
    let verdict = variant_verdict(class, score, limits);

    let mut filters = Vec::new();
    if verdict != Verdict::Supported {
        let most = strata.iter().map(|(_, c)| c.ac).max().unwrap_or(0);
        filters.push(match class {
            Class::NotAssessed => Filter::NoPopulationData,
            Class::Absent => Filter::NeverObserved,
            _ if most <= 1 => Filter::SingleObservation,
            _ if rare_everywhere(&strata, limits.common) => Filter::RareInEveryPopulation,
            _ => Filter::CommonNotDemonstrated,
        });
        match score {
            Some(x) if x <= limits.confirm => filters.push(Filter::LowScore),
            Some(_) => filters.push(Filter::IntermediateScore),
            None if class != Class::NotAssessed => filters.push(Filter::NotScored),
            None => {}
        }
    }
    let site_qc = kgp_site.map(|s| s.qc_failures(ev.kgp_full_total())).unwrap_or_default();
    VariantResult {
        lookup,
        rsid,
        id,
        observed,
        support,
        class,
        basis,
        class_p: corrected.map(|(p, _)| p),
        score,
        verdict,
        filters,
        site_qc,
    }
}

/// Everything concluded about one allele.
pub struct AlleleResult<'a> {
    pub placed: &'a Placed,
    /// Named as trusted (by `--trusted` or the candidate marker).
    pub trusted: bool,
    pub variants: Vec<VariantResult>,
    pub verdict: Verdict,
    pub reason: String,
    pub notes: Vec<String>,
    pub haplotype: Option<hap::Assessment>,
}

impl AlleleResult<'_> {
    /// Trusted and cleanly placed: its variants calibrate the score.
    pub fn calibrates(&self) -> bool {
        self.trusted && self.placed.status == Status::Ok
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
        if parts.is_empty() {
            ".".to_string()
        } else {
            parts.join(",")
        }
    }

    pub fn hap_status(&self) -> &'static str {
        self.haplotype.as_ref().map_or("not_checked", |h| h.status(self.variants.len()))
    }
}

/// `other_build`: a build on which the allele's gene does have coordinates.
pub fn assess_allele<'a>(
    ev: &Evidence,
    placed: &'a Placed,
    trusted: bool,
    other_build: Option<&str>,
    model: &ScoreModel,
    limits: &Limits,
) -> AlleleResult<'a> {
    let variants: Vec<VariantResult> = placed
        .variants
        .iter()
        .map(|c| assess_variant(ev, placed, c, model, limits))
        .collect();
    let total = variants.len();
    let with = |verdict: Verdict| -> Vec<&str> { variants.iter().filter(|v| v.verdict == verdict).map(|v| v.id.as_str()).collect() };
    let weakest = |verdict: Verdict, what: &str| {
        let ids = with(verdict);
        (verdict, format!("{} of {total} variant(s) {what}: {}", ids.len(), ids.join(";")))
    };
    // The verdict of an allele is that of its weakest variant.
    let (verdict, reason) = match placed.status {
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
        notes.push(format!("{flagged} variant(s) at sites failing 1000 Genomes site QC"));
    }

    // Haplotype check: the allele's variants together on one chromosome.
    let haplotype = match (&ev.panel, &placed.placement) {
        (Some(panel), Some(p)) if placed.status != Status::LowIdentity => {
            let keys: Vec<_> = placed.variants.iter().map(|c| &c.key).collect();
            Some(panel.assess(&p.contig, p.span, &keys))
        }
        _ => None,
    };
    AlleleResult {
        placed,
        trusted,
        variants,
        verdict,
        reason,
        notes,
        haplotype,
    }
}

/// Single-base variants carried only by alleles that do not calibrate the
/// score, per support level: where the novelty of the candidates sits.
pub fn candidate_only(results: &[AlleleResult]) -> [u64; BINS] {
    let known: HashSet<&Lookup> = results
        .iter()
        .filter(|r| r.calibrates())
        .flat_map(|r| r.variants.iter().map(|v| &v.lookup))
        .collect();
    let mut seen: HashSet<&Lookup> = HashSet::new();
    let mut bins = [0u64; BINS];
    for r in results.iter().filter(|r| !r.calibrates()) {
        for (v, c) in r.variants.iter().zip(&r.placed.variants) {
            if c.kind == Kind::Snv && !known.contains(&v.lookup) && seen.insert(&v.lookup) {
                bins[score::bin(v.support)] += 1;
            }
        }
    }
    bins
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

    #[test]
    fn rare_needs_a_stratum_large_enough_to_tell() {
        let s = |ac, an| (Stratum::Kgp(0), PopCount { ac, an, het: 0, hom: 0 });
        assert!(rare_everywhere(&[s(3, 100_000), s(0, 50)], 0.01));
        assert!(!rare_everywhere(&[s(3, 100_000), s(40, 5_000)], 0.01));
        assert!(!rare_everywhere(&[s(2, 50)], 0.01));
    }
}

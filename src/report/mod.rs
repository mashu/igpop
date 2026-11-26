//! The plain-language run report, `<prefix>.report.md`: what was found,
//! which alleles to look at, how the score was calibrated in this run and
//! whether the run itself looks healthy.

mod text;

use crate::assess::{AlleleResult, HaplotypeCheck, ScoreModel, ScoreSource, Verdict};
use crate::config::Limits;
use crate::evidence::Evidence;
use crate::locus::Status;
use crate::score::{self, Calibration};
use crate::variant::{Kind, Lookup};
use std::collections::HashSet;

/// Verdicts in the order they are listed.
pub const VERDICTS: [Verdict; 6] = [
    Verdict::Supported,
    Verdict::Inconclusive,
    Verdict::NeedsConfirmation,
    Verdict::NotAssessed,
    Verdict::MatchesReference,
    Verdict::NotPlaced,
];

/// The run, as its steps left it.
pub struct Input<'a> {
    pub prefix: &'a str,
    pub build: &'a str,
    /// Other builds used for genes the main build lacks: the build, its
    /// data, and the number of alleles judged on it.
    pub other_builds: Vec<(&'a str, &'a Evidence, usize)>,
    pub alleles: String,
    pub coords: String,
    pub limits: &'a Limits,
    pub results: &'a [AlleleResult<'a>],
    pub ev: &'a Evidence,
    pub model: &'a ScoreModel,
    pub haplotypes: &'a HaplotypeCheck,
    pub requests: u64,
    pub bytes: u64,
}

/// One allele, as the report shows it.
struct Row {
    allele: String,
    trusted: bool,
    verdict: &'static str,
    reason: String,
    notes: String,
    /// Haplotype check status.
    hap: &'static str,
}

/// The run reduced to what the report prints.
struct View<'a> {
    prefix: &'a str,
    build: &'a str,
    alleles: &'a str,
    coords: &'a str,
    common: f64,
    supported: f64,
    confirm: f64,
    rows: Vec<Row>,
    /// (locus status, number of alleles).
    statuses: Vec<(&'static str, usize)>,
    /// (dataset, where it was read from, what it contributed).
    datasets: Vec<(String, String, String)>,
    requests: u64,
    bytes: u64,
    calibration: Option<&'a Calibration>,
    calibration_source: &'a str,
    /// The two groups behind the scores in use, counted per support
    /// level as (possible changes, real variants).
    levels: Option<[(u64, u64); score::BINS]>,
    /// Single-base variants carried only by candidate alleles, per support level.
    candidate_only: [u64; score::BINS],
    /// Alleles supported by frequency only, by haplotype only, by both.
    routes: [usize; 3],
    /// Haplotype calibration: the limit, the chromosomes it requires, and
    /// (chromosomes, share of arbitrary changes found that often).
    coincidence_limit: f64,
    required: Option<u32>,
    coincidence: Vec<(u32, f64)>,
    warnings: Vec<String>,
}

fn count(input: &Input, status: Status) -> usize {
    input.results.iter().filter(|r| r.placed.status == status).count()
}

/// What each source contributed: (dataset, read from, contribution).
fn datasets(input: &Input) -> Vec<(String, String, String)> {
    let mut out = datasets_of(input.build, input.ev);
    for (build, ev, _) in &input.other_builds {
        out.extend(datasets_of(build, ev));
    }
    out
}

fn datasets_of(build: &str, ev: &Evidence) -> Vec<(String, String, String)> {
    let mut out: Vec<(String, String, String)> = ev
        .datasets
        .iter()
        .map(|counted| {
            let groups: Vec<&str> = counted
                .dataset
                .groups
                .iter()
                .zip(&counted.full)
                .filter(|(_, an)| **an > 0)
                .map(|(group, _)| group.title)
                .collect();
            (
                counted.dataset.title.to_string(),
                counted.dataset.spec.clone(),
                format!(
                    "allele counts in {}, up to {} chromosomes; {} variants in the regions read",
                    match groups.len() {
                        1 => format!("one group ({})", groups[0]),
                        n => format!("{n} groups ({})", groups.join(", ")),
                    },
                    counted.full_all,
                    counted.sites.len()
                ),
            )
        })
        .collect();
    if let Some(spec) = &ev.sources.dbsnp {
        let mut used: Vec<&str> = ev.dbsnp.used().map(|study| study.name.as_str()).collect();
        used.sort_unstable();
        out.push((
            "dbSNP".to_string(),
            spec.clone(),
            format!(
                "rsIDs, COMMON flag; frequencies of {} relayed studies used as evidence ({})",
                used.len(),
                used.join(", ")
            ),
        ));
    }
    if let (Some(spec), Some(panel)) = (&ev.sources.phased, &ev.panel) {
        out.push((
            "phased panel of 1000 Genomes".to_string(),
            spec.clone(),
            format!("haplotype check on {} founder chromosomes", 2 * panel.founders),
        ));
    }
    for row in &mut out {
        row.0 = format!("{build}: {}", row.0);
    }
    out
}

/// Single-base variants carried only by alleles that do not calibrate the
/// score, per support level: where the novelty of the candidates sits.
fn candidate_only(results: &[AlleleResult]) -> [u64; score::BINS] {
    let known: HashSet<&Lookup> = results
        .iter()
        .filter(|r| r.calibrates())
        .flat_map(|r| r.variants.iter().map(|v| &v.lookup))
        .collect();
    let mut seen: HashSet<&Lookup> = HashSet::new();
    let mut bins = [0u64; score::BINS];
    for r in results.iter().filter(|r| !r.calibrates()) {
        for (v, c) in r.variants.iter().zip(&r.placed.variants) {
            if c.kind == Kind::Snv && !known.contains(&v.lookup) && seen.insert(&v.lookup) {
                bins[score::bin(v.support)] += 1;
            }
        }
    }
    bins
}

/// Things about this run that limit how far its results can be relied on.
/// Each says what happened and what to do.
pub fn warnings(input: &Input) -> Vec<String> {
    let (ev, model) = (input.ev, input.model);
    let mut out = Vec::new();
    let with_variants = input.results.iter().any(|r| !r.variants.is_empty());
    if ev.sources.phased.is_some() && ev.panel.is_none() {
        out.push("the phased panel could not be read, so the haplotype check is missing (`hap_status = not_checked`). Running again usually fixes it.".to_string());
    }
    if !ev.counted() && with_variants {
        out.push(format!(
            "no dataset with allele counts was read for {}, so there is no score; a variant can only be `supported` through dbSNP's COMMON flag.",
            input.build
        ));
    }
    if model.source == ScoreSource::BuiltIn {
        out.push(
            "scores come from the built-in curve, not from your alleles. This is the intended behaviour for small sets (D, J, a handful of alleles)."
                .to_string(),
        );
        let mut foreign: Vec<&str> = input
            .results
            .iter()
            .map(|r| r.placed.gene.as_str())
            .filter(|gene| !gene.starts_with("IGH"))
            .collect();
        foreign.sort();
        foreign.dedup();
        if !foreign.is_empty() {
            out.push(format!(
                "the built-in curve was made from IGH genes, but this run holds other genes ({} ...). Variation in other loci may be distributed differently: prefer a run with at least {} trusted variants of the same locus, or pass its `score_model.tsv` with `--calibration`.",
                foreign.iter().take(3).copied().collect::<Vec<_>>().join(", "),
                score::MIN_TRUSTED
            ));
        }
    }
    if let (ScoreSource::Run, Some(calibration)) = (model.source, &model.calibration) {
        let nt = calibration.counts.trusted;
        if nt < 100 {
            out.push(format!(
                "the calibration rests on only {nt} trusted variants, so the score curve is not pinned down well (see the intervals in section 5); verdicts near the limits can change with a few more trusted alleles."
            ));
        }
        let established: u64 = calibration.levels()[score::BINS - 2..].iter().map(|(_, trusted)| trusted).sum();
        if 2 * established < nt {
            out.push(format!(
                "only {established} of {nt} trusted variants reach 1% in any population. Trusted alleles are expected to carry mostly established variants: check that the coordinates match the genome build and that the trusted set (`--trusted`, `--candidate-marker`) really is trusted."
            ));
        }
    }
    let unclean = [Status::AmbiguousLocus, Status::LowIdentity, Status::NoCoordinates].map(|s| count(input, s));
    if unclean.iter().sum::<usize>() > 0 {
        out.push(format!(
            "{} allele(s) were not placed cleanly ({} ambiguous_locus, {} low_identity, {} no_coordinates); see section 3.",
            unclean.iter().sum::<usize>(),
            unclean[0],
            unclean[1],
            unclean[2]
        ));
    }
    if ev.sources.dbsnp.is_some() && ev.dbsnp.rsids.is_empty() && with_variants {
        out.push("dbSNP could not be read: rsIDs, the COMMON flag and relayed frequencies are missing in this run. Running again usually fixes it.".to_string());
    }
    out
}

pub fn render(input: &Input) -> String {
    let view = View {
        prefix: input.prefix,
        build: input.build,
        alleles: &input.alleles,
        coords: &input.coords,
        common: input.limits.common,
        supported: input.limits.supported,
        confirm: input.limits.confirm,
        rows: input
            .results
            .iter()
            .map(|r| Row {
                allele: r.placed.allele.clone(),
                trusted: r.trusted,
                verdict: r.verdict.as_str(),
                reason: r.reason.clone(),
                notes: r.notes.join("; "),
                hap: r.hap_status(),
            })
            .collect(),
        statuses: Status::ALL.iter().map(|s| (s.as_str(), count(input, *s))).collect(),
        datasets: datasets(input),
        requests: input.requests,
        bytes: input.bytes,
        calibration: input.model.calibration.as_ref(),
        calibration_source: &input.model.description,
        levels: input.model.calibration.as_ref().map(Calibration::levels),
        candidate_only: if input.ev.counted() {
            candidate_only(input.results)
        } else {
            [0; score::BINS]
        },
        routes: ["frequency", "haplotype", "frequency+haplotype"]
            .map(|route| input.results.iter().filter(|r| r.supported_by() == route).count()),
        coincidence_limit: input.limits.coincidence,
        required: input.haplotypes.required(input.limits.coincidence),
        coincidence: match &input.haplotypes.coincidence {
            Some(c) => [1, 2, 3, 5, 10, 20, 50].map(|k| (k, c.rate(k))).to_vec(),
            None => Vec::new(),
        },
        warnings: warnings(input),
    };
    text::write_report(&view)
}

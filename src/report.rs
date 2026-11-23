//! The plain-language run report, `<prefix>.report.md`: what was found,
//! which alleles to look at, how the score was calibrated in this run and
//! whether the run itself looks healthy.

use crate::assess::{self, AlleleResult, Class, Limits, ScoreModel, ScoreSource, Verdict};
use crate::evidence::Evidence;
use crate::locus::Status;
use crate::score::{self, Calibration};
use std::fmt::Write as _;

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
    pub alleles: String,
    pub coords: String,
    pub prior: Option<f64>,
    pub limits: &'a Limits,
    pub results: &'a [AlleleResult<'a>],
    pub ev: &'a Evidence,
    pub model: &'a ScoreModel,
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
    prior: Option<f64>,
    rows: Vec<Row>,
    /// (locus status, number of alleles).
    statuses: Vec<(&'static str, usize)>,
    /// (dataset, where it was read from, what it contributed).
    datasets: Vec<(String, String, String)>,
    requests: u64,
    bytes: u64,
    calibration: Option<&'a Calibration>,
    calibration_source: &'a str,
    /// Counts behind the scores in use; `None` if they came from a file.
    counts: Option<([u64; score::BINS], [u64; score::BINS])>,
    /// Single-base variants carried only by candidate alleles, per support level.
    candidate_only: [u64; score::BINS],
    warnings: Vec<String>,
}

fn count(input: &Input, status: Status) -> usize {
    input.results.iter().filter(|r| r.placed.status == status).count()
}

/// What each source contributed: (dataset, read from, contribution).
fn datasets(input: &Input) -> Vec<(String, String, String)> {
    let ev = input.ev;
    let mut out = Vec::new();
    if let Some(spec) = &ev.sources.kgp {
        out.push((
            "1000 Genomes 30x, unrelated samples".to_string(),
            spec.clone(),
            format!("allele counts in 5 populations, {} chromosomes; site quality", ev.kgp_full_total()),
        ));
    }
    if let Some(spec) = &ev.sources.gnomad {
        out.push((
            "gnomAD genomes".to_string(),
            spec.clone(),
            format!("allele counts in 9 genetic ancestry groups, {} chromosomes", ev.gnomad_full_total),
        ));
    }
    if let Some(spec) = &ev.sources.dbsnp {
        let used: Vec<&str> = ev.dbsnp.datasets.iter().filter(|d| d.used).map(|d| d.name.as_str()).collect();
        let mut by_name = used.clone();
        by_name.sort();
        out.push((
            "dbSNP".to_string(),
            spec.clone(),
            format!(
                "rsIDs, COMMON flag; frequencies of {} submitted datasets used as evidence ({}), see dbsnp_datasets.tsv",
                used.len(),
                by_name.join(", ")
            ),
        ));
    }
    if let (Some(spec), Some(panel)) = (&ev.sources.phased, &ev.panel) {
        out.push((
            "1000 Genomes 30x phased panel".to_string(),
            spec.clone(),
            format!("haplotype check on {} founder chromosomes", 2 * panel.founders),
        ));
    }
    out
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
            "scores come from the built-in table, not from your alleles. This is the intended behaviour for small sets (D, J, a handful of alleles)."
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
                "the built-in table was made from IGH genes, but this run holds other genes ({} ...). Variation in other loci may be distributed differently: prefer a run with at least {} trusted variants of the same locus, or pass its table with `--calibration`.",
                foreign.iter().take(3).copied().collect::<Vec<_>>().join(", "),
                score::MIN_TRUSTED
            ));
        }
    }
    if model.source == ScoreSource::Run {
        let nt: u64 = model.trusted.iter().sum();
        if nt < 100 {
            out.push(format!(
                "the calibration rests on only {nt} trusted variants, so scores are coarse; verdicts near the limits can change with a few more trusted alleles."
            ));
        }
        let established: u64 = model.trusted[score::BINS - 2..].iter().sum();
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
        prior: input.prior,
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
        counts: input.model.counts(),
        candidate_only: if input.ev.counted() {
            assess::candidate_only(input.results)
        } else {
            [0; score::BINS]
        },
        warnings: warnings(input),
    };
    write_report(&view)
}

fn percent(x: f64) -> String {
    if x >= 0.995 {
        ">99%".to_string()
    } else if x < 0.005 {
        "<1%".to_string()
    } else {
        format!("{:.0}%", 100.0 * x)
    }
}

/// The verdict a variant with this score gets when it is not classed common.
fn band(score: f64, supported: f64, confirm: f64) -> &'static str {
    let limits = Limits {
        common: f64::NAN,
        supported,
        confirm,
    };
    assess::variant_verdict(Class::Uncommon, Some(score), &limits).as_str()
}

fn write_report(r: &View) -> String {
    let mut o = String::new();
    let verdicts = VERDICTS.map(|v| v.as_str());
    let count = |verdict: &str, trusted: Option<bool>| {
        r.rows
            .iter()
            .filter(|x| x.verdict == verdict && trusted.is_none_or(|t| x.trusted == t))
            .count()
    };
    let _ = writeln!(o, "# igpop report: {}\n", r.prefix);
    let _ = writeln!(
        o,
        "igpop {} | alleles: `{}` | coordinates: {} | genome build: {}\n",
        env!("CARGO_PKG_VERSION"),
        r.alleles,
        r.coords,
        r.build
    );

    let _ = writeln!(o, "## 1. Result at a glance\n");
    let _ = writeln!(o, "| verdict | alleles | trusted | candidates | meaning |");
    let _ = writeln!(o, "|---|---:|---:|---:|---|");
    let meaning = [
        "every difference from the reference genome is an established population variant",
        "at least one difference is seen in people, but too rarely to tell real from error",
        "at least one difference looks like an arbitrary base change: confirm independently",
        "no population data could be read for at least one difference",
        "no difference from the reference genome, nothing to assess",
        "the allele could not be aligned to a locus; see section 3",
    ];
    for (verdict, text) in verdicts.iter().zip(meaning) {
        let _ = writeln!(
            o,
            "| `{verdict}` | {} | {} | {} | {text} |",
            count(verdict, None),
            count(verdict, Some(true)),
            count(verdict, Some(false))
        );
    }
    let _ = writeln!(
        o,
        "\nAn allele takes the verdict of its weakest difference. \"Trusted\" alleles are the ones \
         igpop used as examples of real variation; \"candidates\" are the ones being judged.\n"
    );

    let _ = writeln!(o, "## 2. Alleles to look at\n");
    for (verdict, title) in [
        ("needs_confirmation", "Needs confirmation"),
        ("inconclusive", "Inconclusive"),
        ("not_assessed", "Not assessed"),
    ] {
        let rows: Vec<&Row> = r.rows.iter().filter(|x| x.verdict == verdict).collect();
        let _ = writeln!(o, "### {title} ({})\n", rows.len());
        if rows.is_empty() {
            let _ = writeln!(o, "None.\n");
            continue;
        }
        let _ = writeln!(o, "| allele | trusted | why | haplotype check | notes |");
        let _ = writeln!(o, "|---|---|---|---|---|");
        for x in rows {
            let _ = writeln!(
                o,
                "| {} | {} | {} | `{}` | {} |",
                x.allele,
                if x.trusted { "yes" } else { "no" },
                x.reason,
                x.hap,
                x.notes
            );
        }
        let _ = writeln!(o);
    }
    let _ = writeln!(
        o,
        "Variants are written as `rsID(reference>allele:position in the allele,position in the genome)`. \
         Every one of them has a full row of evidence in `{}.variants.tsv`.\n",
        r.prefix
    );

    let _ = writeln!(o, "## Haplotype check\n");
    if r.rows.iter().all(|x| x.hap == "not_checked") {
        let _ = writeln!(
            o,
            "Not done in this run (no phased panel for this build, or it was switched off).\n"
        );
    } else {
        let _ = writeln!(
            o,
            "Separate from the verdict: is the allele's whole combination of variants carried on one chromosome \
             in the phased 1000 Genomes panel? The verdict looks at each variant on its own.\n"
        );
        let _ = writeln!(o, "| haplotype check | alleles | trusted | candidates | meaning |");
        let _ = writeln!(o, "|---|---:|---:|---:|---|");
        for (status, text) in [
            (
                "found_exact",
                "chromosomes with exactly this sequence over the allele exist (`hap_exact` says how many)",
            ),
            (
                "found_with_extra",
                "the combination exists, but always together with further variants inside the allele",
            ),
            (
                "not_found",
                "every variant is in the panel, but no chromosome carries them all together",
            ),
            (
                "incomplete",
                "at least one variant is not a site of the panel, so the full combination cannot be looked up",
            ),
            ("untestable", "none of the variants is a site of the panel"),
            ("not_checked", "not placed, or no phased data"),
        ] {
            let k = |trusted: Option<bool>| {
                r.rows
                    .iter()
                    .filter(|x| x.hap == status && trusted.is_none_or(|t| x.trusted == t))
                    .count()
            };
            let _ = writeln!(o, "| `{status}` | {} | {} | {} | {text} |", k(None), k(Some(true)), k(Some(false)));
        }
        let _ = writeln!(
            o,
            "\n`found_exact` on more than a few chromosomes is strong support. The other outcomes are weak evidence against an allele: the panel leaves out \
             variants seen once and sites that fail quality filters, holds a few thousand chromosomes, and rests on short reads. \
             Counts and frequencies per superpopulation are in the `hap_*` columns of `{}.alleles.tsv`.\n",
            r.prefix
        );
    }

    let _ = writeln!(o, "## 3. Was every allele placed on its gene?\n");
    let _ = writeln!(o, "| placement | alleles |");
    let _ = writeln!(o, "|---|---:|");
    for (status, n) in &r.statuses {
        let _ = writeln!(o, "| `{status}` | {n} |");
    }
    let unplaced: Vec<&Row> = r.rows.iter().filter(|x| x.verdict == "not_placed").collect();
    if !unplaced.is_empty() {
        let _ = writeln!(o, "\nNot placed:\n");
        let mut reasons: Vec<&str> = unplaced.iter().map(|x| x.reason.as_str()).collect();
        reasons.sort();
        reasons.dedup();
        for reason in reasons {
            let names: Vec<&str> = unplaced.iter().filter(|x| x.reason == reason).map(|x| x.allele.as_str()).collect();
            let shown = names.iter().take(8).copied().collect::<Vec<_>>().join(", ");
            let more = if names.len() > 8 { ", ..." } else { "" };
            let _ = writeln!(o, "- {} allele(s), {reason}: {shown}{more}", names.len());
        }
    }
    let _ = writeln!(
        o,
        "\n`ok`: aligned to the locus of its gene. `ambiguous_locus`: fits another locus at least as well, \
         so its variants may belong to the other gene. `low_identity`: too different from the locus to call \
         variants. `no_coordinates`: its gene has no row in the coordinates table. \
         The alignments themselves are in `{}.alignments.txt`.\n",
        r.prefix
    );

    let _ = writeln!(o, "## 4. Population data used\n");
    if r.datasets.is_empty() {
        let _ = writeln!(o, "None.\n");
    } else {
        let _ = writeln!(o, "| dataset | read from | contributed |");
        let _ = writeln!(o, "|---|---|---|");
        for (name, spec, detail) in &r.datasets {
            let _ = writeln!(o, "| {name} | `{spec}` | {detail} |");
        }
        let _ = writeln!(
            o,
            "\nDownloaded in this run: {:.2} MB in {} requests (regions already in the cache are not downloaded again).\n",
            r.bytes as f64 / 1e6,
            r.requests
        );
    }

    let _ = writeln!(o, "## 5. How the score was calibrated\n");
    match r.calibration {
        None => {
            let _ = writeln!(
                o,
                "No score in this run: no dataset with allele counts was read. Verdicts rest on the `class` column alone.\n"
            );
        }
        Some(cal) => {
            let _ = writeln!(o, "Source: {}.\n", r.calibration_source);
            let _ = writeln!(
                o,
                "The score of a variant says how many times more often its level of population support is found \
                 among variants of trusted alleles than among arbitrary single-base changes. \
                 Above 1 the data speak for a real variant, below 1 against.\n"
            );
            let _ = writeln!(
                o,
                "| population support (frequency demonstrated in the best population) | arbitrary changes | trusted variants | score | verdict | probability real if half of such variants are real | if 9 in 10 are real |"
            );
            let _ = writeln!(o, "|---|---:|---:|---:|---|---:|---:|");
            for i in 0..score::BINS {
                let (a, t) = match &r.counts {
                    Some((a, t)) => (a[i].to_string(), t[i].to_string()),
                    None => (".".to_string(), ".".to_string()),
                };
                let _ = writeln!(
                    o,
                    "| {} | {a} | {t} | {:.3} | `{}` | {} | {} |",
                    score::LABELS[i],
                    cal.score[i],
                    band(cal.score[i], r.supported, r.confirm),
                    percent(score::posterior(0.5, cal.score[i])),
                    percent(score::posterior(0.9, cal.score[i])),
                );
            }
            let _ = writeln!(
                o,
                "\nVerdict limits: `supported` at score {} or more, `needs_confirmation` at {} or less \
                 (`--supported-score`, `--confirm-score`). A variant classed `common` (frequency of at least {}% \
                 demonstrated in one population, or dbSNP's COMMON flag for this exact allele) is `supported` whatever its score.",
                r.supported,
                r.confirm,
                100.0 * r.common
            );
            match r.prior {
                Some(prior) => {
                    let _ = writeln!(
                        o,
                        "The `probability_real` column uses your prior of {prior} (`--prior`): the share of variants like these you expect to be real before looking at population data.\n"
                    );
                }
                None => {
                    let _ = writeln!(
                        o,
                        "The two probability columns show how the same score reads under two different expectations; \
                         pass `--prior` to get a `probability_real` column for your own.\n"
                    );
                }
            }
        }
    }

    let _ = writeln!(o, "## 6. Checks on this run\n");
    let trusted_total = r.rows.iter().filter(|x| x.trusted).count();
    if trusted_total > 0 {
        let parts: Vec<String> = verdicts
            .iter()
            .filter(|v| count(v, Some(true)) > 0)
            .map(|v| format!("{} `{v}`", count(v, Some(true))))
            .collect();
        let _ = writeln!(
            o,
            "- **Trusted alleles, judged like candidates:** {} of {trusted_total}. \
             Most trusted alleles with differences should come out `supported`; if they do not, \
             check the coordinates and the choice of trusted alleles before reading the candidates.",
            parts.join(", ")
        );
    }
    let only: u64 = r.candidate_only.iter().sum();
    if only > 0 {
        let parts: Vec<String> = (0..score::BINS)
            .filter(|i| r.candidate_only[*i] > 0)
            .map(|i| format!("{} at {}", r.candidate_only[i], score::LABELS[i]))
            .collect();
        let _ = writeln!(
            o,
            "- **Single-base variants found only in candidate alleles, by population support:** {}. \
             Arbitrary errors pile up at `none`; real but rare variants sit in between.",
            parts.join(", ")
        );
    }
    if r.warnings.is_empty() {
        let _ = writeln!(o, "- **Warnings:** none.");
    } else {
        for w in &r.warnings {
            let _ = writeln!(o, "- **Warning:** {w}");
        }
    }

    let _ = writeln!(o, "\n## 7. Files of this run\n");
    for (file, text) in [
        (
            "alleles.tsv",
            "every allele, with verdict, the reasons it did not pass (`filters`) and the haplotype check",
        ),
        (
            "final.tsv",
            "only the alleles that passed (`supported` or `matches_reference`), same columns",
        ),
        (
            "variants.tsv",
            "one row per difference from the reference genome, with all population evidence",
        ),
        ("alignments.txt", "alignments stacked per gene, to check by eye"),
        (
            "alignments.tsv",
            "placement of every allele: position, strand, identity, competing loci",
        ),
        (
            "score_calibration.tsv",
            "the calibration table of section 5; reusable with --calibration",
        ),
        (
            "dbsnp_datasets.tsv",
            "frequency datasets relayed by dbSNP and whether each was used",
        ),
    ] {
        let _ = writeln!(o, "- `{}.{file}`: {text}", r.prefix);
    }
    o
}

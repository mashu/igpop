//! The plain-language run report, `<prefix>.report.md`: what was found,
//! which alleles to look at, how the score was calibrated in this run and
//! whether the run itself looks healthy.

use crate::assess::{self, AlleleResult, Class, HaplotypeCheck, Limits, ScoreModel, ScoreSource, Verdict};
use crate::evidence::Evidence;
use crate::freq::KgpRelease;
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
    /// Counts behind the scores in use; `None` if they came from a file.
    counts: Option<([u64; score::BINS], [u64; score::BINS])>,
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
    let mut out = Vec::new();
    let phase3 = ev.sources.kgp_release == KgpRelease::Phase3;
    if let Some(spec) = &ev.sources.kgp {
        out.push((
            if phase3 {
                "1000 Genomes phase 3".to_string()
            } else {
                "1000 Genomes 30x, unrelated samples".to_string()
            },
            spec.clone(),
            format!(
                "allele counts in 5 populations, {} chromosomes{}",
                ev.kgp_full_total(),
                if phase3 { "" } else { "; site quality" }
            ),
        ));
    }
    if let Some(spec) = &ev.sources.gnomad {
        out.push((
            format!("gnomAD genomes ({})", ev.sources.gnomad_name),
            spec.clone(),
            format!(
                "allele counts in {} genetic ancestry groups, {} chromosomes",
                ev.gnomad_full.iter().filter(|an| **an > 0).count(),
                ev.gnomad_full_total
            ),
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
                "rsIDs, COMMON flag; frequencies of {} submitted datasets used as evidence ({})",
                used.len(),
                by_name.join(", ")
            ),
        ));
    }
    if let (Some(spec), Some(panel)) = (&ev.sources.phased, &ev.panel) {
        out.push((
            if phase3 {
                "1000 Genomes phase 3, phased genotypes".to_string()
            } else {
                "1000 Genomes 30x phased panel".to_string()
            },
            spec.clone(),
            format!("haplotype check on {} founder chromosomes", 2 * panel.founders),
        ));
    }
    for row in &mut out {
        row.0 = format!("{build}: {}", row.0);
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
    if model.source == ScoreSource::Run {
        let nt: u64 = model.trusted.iter().sum();
        if nt < 100 {
            out.push(format!(
                "the calibration rests on only {nt} trusted variants, so the score curve is not pinned down well (see the intervals in section 5); verdicts near the limits can change with a few more trusted alleles."
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
    write_report(&view)
}

#[allow(dead_code)]
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
        coincidence: f64::NAN,
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
        "\nBy frequency, an allele takes the verdict of its weakest difference. An allele whose exact sequence \
         is established on phased chromosomes is `supported` whatever its single differences say. \
         Of the supported alleles, {} are supported by frequency only, {} by haplotype only and {} by both \
         (column `supported_by`). \"Trusted\" alleles are the ones named as examples of real variation; \
         \"candidates\" are the ones being judged.\n",
        r.routes[0], r.routes[1], r.routes[2]
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
            "\n**How much is a \"found\" worth?** igpop takes every founder chromosome of the panel, changes one base at a time, \
             and looks the changed sequence up among the chromosomes. This is what an error on a real allele would look like:\n"
        );
        let _ = writeln!(o, "| found on at least | share of arbitrary changes |");
        let _ = writeln!(o, "|---:|---:|");
        for (chromosomes, rate) in &r.coincidence {
            let _ = writeln!(o, "| {chromosomes} chromosome(s) | {:.2}% |", 100.0 * rate);
        }
        match r.required {
            Some(required) => {
                let _ = writeln!(
                    o,
                    "\nAn allele is supported by its haplotype when its exact sequence is on at least {required} founder chromosomes: \
                     fewer than {}% of arbitrary changes are found that often (`--coincidence`), and never on a single chromosome, \
                     which can be a genotyping error.",
                    100.0 * r.coincidence_limit
                );
            }
            None => {
                let _ = writeln!(
                    o,
                    "\nNo number of chromosomes reaches the limit of {}% in this run.",
                    100.0 * r.coincidence_limit
                );
            }
        }
        let _ = writeln!(
            o,
            "\nThe outcomes other than `found_exact` say little against an allele: the panel leaves out \
             variants seen once and sites that fail quality filters, holds a few thousand chromosomes, and rests on short reads. \
             For an allele with a single difference, the haplotype check and the frequency of that difference are largely the same evidence. \
             Counts per superpopulation are in the `hap_*` columns of `{}.alleles.tsv`.\n",
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

    // Alleles of genes the main build does not have on a chromosome.
    let apart = |start: &str| -> Vec<&Row> { r.rows.iter().filter(|x| x.notes.starts_with(start)).collect() };
    let (other_build, scaffold) = (apart("gene absent from"), apart("gene on "));
    if !other_build.is_empty() || !scaffold.is_empty() {
        let _ = writeln!(o, "### Genes that {} does not have on a chromosome\n", r.build);
        let genes = |rows: &[&Row]| {
            let mut genes: Vec<&str> = rows.iter().map(|x| x.allele.split('*').next().unwrap_or(&x.allele)).collect();
            genes.sort();
            genes.dedup();
            genes.join(", ")
        };
        if !other_build.is_empty() {
            let _ = writeln!(
                o,
                "- {} allele(s) of {}. These genes lie on a stretch of DNA that some people carry and others \
                 do not, and the {} reference represents a chromosome without it. Their alleles were placed on \
                 another build and judged with the population data of that build (section 4) and the score of \
                 section 5; the column `build` and the notes of each allele say which. Those data are older \
                 and smaller, and people without the gene contribute no reads, so read these verdicts with \
                 more caution than the others.",
                other_build.len(),
                genes(&other_build),
                r.build
            );
        }
        if !scaffold.is_empty() {
            let _ = writeln!(
                o,
                "- {} allele(s) of {}: the gene is in {} only on a piece of sequence that is not part of a \
                 chromosome. The alleles are compared with that sequence, but the population call sets are \
                 not reliable there and are not read: an allele that differs from it is `not_assessed`.",
                scaffold.len(),
                genes(&scaffold),
                r.build
            );
        }
        let _ = writeln!(o);
    }

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
                "The score of a variant says how many times more typical its population support is of the variants of \
                 trusted alleles than of arbitrary single-base changes. Above 1 the data speak for a real variant, below 1 against. \
                 It is a smooth curve of support, so two variants with nearly the same support get nearly the same score.\n"
            );
            let c = cal.counts;
            let _ = writeln!(
                o,
                "The curve was placed with {} trusted variants ({} of them seen in some dataset) and {} possible single-base changes ({} seen). \
                 A change that no dataset has seen scores {:.3}. For a change that was seen, with x = log10(support): \
                 log(score) = {:.3} {} {:.3} x {} {:.3} x^2.\n",
                c.trusted,
                c.trusted_seen,
                c.arbitrary,
                c.arbitrary_seen,
                cal.never,
                cal.beta[0],
                if cal.beta[1] < 0.0 { "-" } else { "+" },
                cal.beta[1].abs(),
                if cal.beta[2] < 0.0 { "-" } else { "+" },
                cal.beta[2].abs(),
            );
            let _ = writeln!(
                o,
                "| population support (frequency demonstrated in the best population) | score | 95% interval of the score | arbitrary changes supported at least as well | verdict |"
            );
            let _ = writeln!(o, "|---|---:|---:|---:|---|");
            for support in [0.0, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3] {
                let (x, (low, high)) = (cal.score(support), cal.interval(support));
                let (label, chance) = (
                    if support == 0.0 {
                        "never seen".to_string()
                    } else {
                        format!("{}%", 100.0 * support)
                    },
                    format!("{:.2}%", 100.0 * cal.chance(support)),
                );
                let _ = writeln!(
                    o,
                    "| {label} | {x:.3} | {low:.3} to {high:.3} | {chance} | `{}` |",
                    band(x, r.supported, r.confirm),
                );
            }
            let crossing = |limit: f64| match cal.support_at(limit) {
                Some(s) if s < 1e-3 => format!("a support of {:.3}%", 100.0 * s),
                Some(s) => format!("a support of {:.1}%", 100.0 * s),
                None => "no support (the curve does not reach it)".to_string(),
            };
            let _ = writeln!(
                o,
                "\nThe score reaches {} at {}, and {} at {}. The interval says how well the trusted variants pin the curve down; \
                 the last column says how often an arbitrary change does as well as the variant.\n",
                r.supported,
                crossing(r.supported),
                r.confirm,
                crossing(r.confirm),
            );
            if let Some((a, t)) = &r.counts {
                let _ = writeln!(
                    o,
                    "The two groups counted per level of support (a ruler for reading; the levels do not enter the score):\n"
                );
                let _ = writeln!(o, "| support | arbitrary changes | trusted variants |\n|---|---:|---:|");
                for i in 0..score::BINS {
                    let _ = writeln!(o, "| {} | {} | {} |", score::LABELS[i], a[i], t[i]);
                }
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
            let _ = writeln!(
                o,
                "A score is evidence, not a probability: how likely a variant is to be real also depends on what was known before looking at population data.\n"
            );
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
        ("score_model.tsv", "the score curve of section 5; reusable with --calibration"),
        ("score_levels.tsv", "the two groups counted per level of support"),
        (
            "dbsnp_datasets.tsv",
            "frequency datasets relayed by dbSNP and whether each was used",
        ),
        (
            "haplotype_calibration.tsv",
            "how often an arbitrary change is found on phased chromosomes",
        ),
    ] {
        let _ = writeln!(o, "- `{}.{file}`: {text}", r.prefix);
    }
    o
}

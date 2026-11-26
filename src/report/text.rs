//! The report as Markdown text.

use super::{Row, View, VERDICTS};
use crate::assess::{variant_verdict, Class};
use crate::config::Limits;
use crate::score;
use std::fmt::Write as _;

/// The verdict a variant with this score gets when it is not classed common.
fn band(score: f64, supported: f64, confirm: f64) -> &'static str {
    let limits = Limits {
        common: f64::NAN,
        supported,
        confirm,
        coincidence: f64::NAN,
    };
    variant_verdict(Class::Uncommon, Some(score), &limits).as_str()
}

pub(super) fn write_report(r: &View) -> String {
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
            if let Some(levels) = &r.levels {
                let _ = writeln!(
                    o,
                    "The two groups counted per level of support (a ruler for reading; the levels do not enter the score):\n"
                );
                let _ = writeln!(o, "| support | arbitrary changes | trusted variants |\n|---|---:|---:|");
                for (label, (arbitrary, trusted)) in score::LABELS.iter().zip(levels) {
                    let _ = writeln!(o, "| {label} | {arbitrary} | {trusted} |");
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

//! `alleles.tsv` and `final.tsv`: one row per allele.

use super::{dots, list, names, num, Table};
use crate::assess::{AlleleResult, Class, Yardstick};
use crate::population::KGP_GROUPS;
use crate::stats::{self, PopCount};

fn per_pop(prefix: &str) -> Vec<String> {
    KGP_GROUPS.iter().map(|g| format!("{prefix}{}", g.code)).collect()
}

fn counts(pops: &[PopCount]) -> Vec<String> {
    pops.iter().map(|c| format!("{}/{}", c.ac, c.an)).collect()
}

/// `alleles.tsv` and `final.tsv`: one row per allele.
pub fn alleles_table<'a>() -> Table<AlleleResult<'a>> {
    type R<'a> = AlleleResult<'a>;
    let of = |r: &R, class: Class| r.ids(|v| v.class == class);
    // "Uncommon" keeps IgSNPer's meaning: every variant not classed common.
    let not_common = move |r: &R| -> Vec<String> {
        [Class::Uncommon, Class::Absent, Class::NotAssessed]
            .iter()
            .flat_map(|c| of(r, *c))
            .collect()
    };
    let missing = |r: &R, end: usize| match &r.placed.placement {
        Some(p) => [p.missing.0, p.missing.1][end].to_string(),
        None => ".".to_string(),
    };
    Table::new()
        .col("allele_name", |r: &R| r.placed.allele.clone())
        .col("gene", |r: &R| r.placed.gene.clone())
        .col("build", |r: &R| r.build.to_string())
        .col("verdict", |r: &R| r.verdict.as_str().to_string())
        .col("supported_by", |r: &R| r.supported_by().to_string())
        .col("filters", |r: &R| r.filters())
        .col("reason", |r: &R| r.reason.clone())
        .col(
            "notes",
            |r: &R| if r.notes.is_empty() { ".".to_string() } else { r.notes.join("; ") },
        )
        .col("hap_status", |r: &R| r.hap_status().to_string())
        .col("locus_status", |r: &R| r.placed.status.as_str().to_string())
        .col("shares_locus_with", |r: &R| r.placed.shared())
        .col("trusted", |r: &R| {
            if r.yardstick == Yardstick::Named { "yes" } else { "no" }.to_string()
        })
        .col("yardstick", |r: &R| r.yardstick.as_str().to_string())
        .col("missing_5prime", move |r: &R| missing(r, 0))
        .col("missing_3prime", move |r: &R| missing(r, 1))
        .col("n_variants", |r: &R| r.variants.len().to_string())
        .col("min_score", |r: &R| r.weakest().map_or(".".to_string(), |(x, _)| format!("{x:.3}")))
        .col("min_score_variant", |r: &R| {
            r.weakest().map_or(".".to_string(), |(_, id)| id.to_string())
        })
        .col("common", move |r: &R| of(r, Class::Common).len().to_string())
        .col("common_str", move |r: &R| list(&of(r, Class::Common)))
        .col("uncommon", move |r: &R| not_common(r).len().to_string())
        .col("uncommon_str", move |r: &R| list(&not_common(r)))
        .col("absent", move |r: &R| of(r, Class::Absent).len().to_string())
        .col("absent_str", move |r: &R| list(&of(r, Class::Absent)))
        .col("not_assessed", move |r: &R| of(r, Class::NotAssessed).len().to_string())
        .col("site_qc_flagged", |r: &R| r.ids(|v| !v.site_qc.is_empty()).len().to_string())
        .col("site_qc_flagged_str", |r: &R| list(&r.ids(|v| !v.site_qc.is_empty())))
        .block(haplotype_names(), haplotype_cells)
}

fn haplotype_names() -> Vec<String> {
    let mut out = names(&[
        "hap_tested",
        "hap_untested_str",
        "hap_total",
        "hap_exact",
        "hap_exact_donors",
        "hap_coincidence",
        "hap_exact_AF",
        "hap_exact_AF_lower95",
        "hap_exact_AF_upper95",
    ]);
    out.extend(per_pop("hap_exact_"));
    out.extend(names(&[
        "hap_exact_max_pop",
        "hap_exact_pop_heterogeneity_p",
        "hap_carrying",
        "hap_carrying_AF",
    ]));
    out.extend(per_pop("hap_carrying_"));
    out
}

/// The haplotype check. `hap_exact_*` only when every variant was tested;
/// `hap_carrying_*` only when at least one was (or the allele has none).
fn haplotype_cells(r: &AlleleResult) -> Vec<String> {
    let Some(h) = &r.haplotype else {
        return dots(haplotype_names().len());
    };
    let total = r.variants.len();
    let untested: Vec<String> = h.untested.iter().map(|i| r.variants[*i].id.clone()).collect();
    let founders = stats::total(&h.carrying).an;
    let mut out = vec![format!("{}/{total}", total - untested.len()), list(&untested), founders.to_string()];
    match &h.exact {
        Some(exact) => {
            let t = stats::total(exact);
            out.extend([
                t.ac.to_string(),
                h.exact_donors.map_or(".".to_string(), |d| d.to_string()),
                r.haplotype_coincidence.map_or(".".to_string(), num),
                num(t.af()),
                num(t.lower95()),
                num(t.upper95()),
            ]);
            out.extend(counts(exact));
            out.push(if t.ac > 0 {
                KGP_GROUPS[stats::strongest(exact)].code.to_string()
            } else {
                ".".to_string()
            });
            out.push(num(stats::heterogeneity_p(exact)));
        }
        None => out.extend(dots(13)),
    }
    if total > 0 && untested.len() == total {
        out.extend(dots(7));
    } else {
        let carrying = stats::total(&h.carrying);
        out.extend([carrying.ac.to_string(), num(carrying.af())]);
        out.extend(counts(&h.carrying));
    }
    out
}

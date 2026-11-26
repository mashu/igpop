//! `variants.tsv`: one row per difference of an allele, with what was
//! concluded about it. The counts behind the conclusion are in
//! `populations.tsv`.

use super::{num, Table};
use crate::assess::{AlleleResult, VariantResult};
use crate::evidence::Evidence;
use crate::locus::Called;

/// One row of `variants.tsv`.
pub struct VariantRow<'a> {
    pub ev: &'a Evidence,
    pub allele: &'a AlleleResult<'a>,
    pub called: &'a Called,
    pub v: &'a VariantResult,
}

pub fn variants_table<'a>() -> Table<VariantRow<'a>> {
    type R<'a> = VariantRow<'a>;
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    Table::new()
        .col("allele_name", |r: &R| r.allele.placed.allele.clone())
        .col("gene", |r: &R| r.allele.placed.gene.clone())
        .col("build", |r: &R| r.allele.build.to_string())
        .col("rsid", |r: &R| r.v.rsid.clone())
        .col("evidence", |r: &R| {
            if r.v.observed.is_empty() {
                "-".to_string()
            } else {
                r.v.observed.join(";")
            }
        })
        .col("contig", |r: &R| r.allele.placed.contig().to_string())
        .col("genome_pos", |r: &R| r.called.key.0.to_string())
        .col("aln_pos", |r: &R| r.called.allele_pos.to_string())
        .col("aln_pos_genome_strand", |r: &R| r.called.strand_pos.to_string())
        .col("strand", |r: &R| if r.called.reverse { "-" } else { "+" }.to_string())
        .col("type", |r: &R| r.called.kind.as_str().to_string())
        .col("reference_bases", move |r: &R| text(&r.called.key.1))
        .col("query_variant", move |r: &R| text(&r.called.key.2))
        .col("verdict", |r: &R| r.v.verdict.as_str().to_string())
        .col("filters", |r: &R| {
            if r.v.filters.is_empty() {
                ".".to_string()
            } else {
                r.v.filters.iter().map(|f| f.as_str()).collect::<Vec<_>>().join(",")
            }
        })
        .col("data_show", |r: &R| r.v.shows.map_or(".".to_string(), |f| f.as_str().to_string()))
        .col("score", |r: &R| r.v.score.map_or(".".to_string(), |x| format!("{x:.3}")))
        .col("score_basis", |r: &R| {
            match (r.v.score, r.v.score_borrowed) {
                (None, _) => ".",
                (Some(_), false) => "curve",
                (Some(_), true) => "borrowed",
            }
            .to_string()
        })
        .col("score_low", |r: &R| {
            r.v.score_interval.map_or(".".to_string(), |x| format!("{:.3}", x.0))
        })
        .col("score_high", |r: &R| {
            r.v.score_interval.map_or(".".to_string(), |x| format!("{:.3}", x.1))
        })
        .col("chance", |r: &R| r.v.chance.map_or(".".to_string(), num))
        .col("support", |r: &R| {
            if r.ev.reads(&r.v.lookup.0) {
                num(r.v.support)
            } else {
                ".".to_string()
            }
        })
        .col("support_from", |r: &R| r.v.support_from.clone().unwrap_or_else(|| ".".to_string()))
        .col("class", |r: &R| r.v.class.as_str().to_string())
        .col("class_basis", |r: &R| r.v.basis.clone().unwrap_or_else(|| ".".to_string()))
        .col("class_p", |r: &R| r.v.class_p.map_or(".".to_string(), num))
        .col("site_quality", |r: &R| {
            if r.v.site_qc.is_empty() {
                ".".to_string()
            } else {
                r.v.site_qc
                    .iter()
                    .map(|(dataset, failed)| format!("{dataset}:{}", failed.join(",")))
                    .collect::<Vec<_>>()
                    .join(";")
            }
        })
        .col("dbsnp_common", |r: &R| {
            match (r.ev.sources.dbsnp.is_some(), r.ev.dbsnp.common.contains(&r.v.lookup)) {
                (false, _) => ".",
                (true, true) => "yes",
                (true, false) => "no",
            }
            .to_string()
        })
        .col("locus_status", |r: &R| r.allele.placed.status.as_str().to_string())
        .col("shares_locus_with", |r: &R| r.allele.placed.shared())
}

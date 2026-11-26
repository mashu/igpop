//! The population data behind the verdicts, in three tables that share
//! their key columns with `variants.tsv`:
//!
//! - `populations.tsv`: one row per difference and population that looked
//!   at its position, with the counts and their bounds;
//! - `site_quality.tsv`: one row per difference and dataset that has a
//!   record of it, with the dataset's quality measures;
//! - `datasets.tsv`: one row per population of the run, saying who it is.
//!
//! A difference shared by several alleles is listed once.

use super::{num, or_dot, Table};
use crate::catalogue;
use crate::evidence::{Counted, Evidence, Stratum};
use crate::population::dbsnp::Study;
use crate::population::{Group, Record};
use crate::stats::PopCount;
use crate::variant::Lookup;

/// One difference from the reference, on the build it was judged on.
pub struct Site<'a> {
    pub build: &'a str,
    pub ev: &'a Evidence,
    pub lookup: &'a Lookup,
    pub rsid: &'a str,
}

/// The columns that say which difference a row is about.
fn key<'a, T: 'a>(table: Table<T>, site: impl Fn(&T) -> &Site<'a> + Copy + 'static) -> Table<T> {
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    table
        .col("build", move |r: &T| site(r).build.to_string())
        .col("contig", move |r: &T| site(r).lookup.0.clone())
        .col("genome_pos", move |r: &T| site(r).lookup.1 .0.to_string())
        .col("reference_bases", move |r: &T| text(&site(r).lookup.1 .1))
        .col("query_variant", move |r: &T| text(&site(r).lookup.1 .2))
        .col("rsid", move |r: &T| site(r).rsid.to_string())
}

/// One row of `populations.tsv`.
pub struct PopulationRow<'a> {
    pub site: &'a Site<'a>,
    pub stratum: Stratum,
    pub count: PopCount,
}

impl<'a> PopulationRow<'a> {
    /// The rows of one difference: every population that looked at it.
    pub fn of(site: &'a Site<'a>) -> impl Iterator<Item = PopulationRow<'a>> + 'a {
        site.ev
            .examined(site.lookup)
            .into_iter()
            .map(move |(stratum, count)| PopulationRow { site, stratum, count })
    }
}

pub fn populations_table<'a>() -> Table<PopulationRow<'a>> {
    type R<'a> = PopulationRow<'a>;
    key(Table::new(), |r: &R| r.site)
        .col("dataset", |r: &R| r.site.ev.describe(r.stratum).dataset.to_string())
        .col("group", |r: &R| r.site.ev.describe(r.stratum).group)
        .col("counts", |r: &R| {
            if r.site.ev.describe(r.stratum).estimated {
                "estimated"
            } else {
                "reported"
            }
            .to_string()
        })
        .col("AC", |r: &R| r.count.ac.to_string())
        .col("AN", |r: &R| r.count.an.to_string())
        .col("AF", |r: &R| num(r.count.af()))
        .col("AF_lower95", |r: &R| num(r.count.lower95()))
        .col("AF_upper95", |r: &R| num(r.count.upper95()))
}

/// One row of `site_quality.tsv`.
pub struct QualityRow<'a> {
    pub site: &'a Site<'a>,
    pub counted: &'a Counted,
    pub record: &'a Record,
}

impl<'a> QualityRow<'a> {
    /// The rows of one difference: every dataset with a record of it.
    pub fn of(site: &'a Site<'a>) -> impl Iterator<Item = QualityRow<'a>> + 'a {
        site.ev
            .datasets
            .iter()
            .filter_map(move |counted| counted.sites.get(site.lookup).map(|record| QualityRow { site, counted, record }))
    }
}

pub fn quality_table<'a>() -> Table<QualityRow<'a>> {
    type R<'a> = QualityRow<'a>;
    let join = |items: Vec<String>| if items.is_empty() { ".".to_string() } else { items.join(",") };
    key(Table::new(), |r: &R| r.site)
        .col("dataset", |r: &R| r.counted.dataset.name.to_string())
        .col("failed", move |r: &R| match r.counted.failures(r.site.lookup) {
            failed if failed.is_empty() => "PASS".to_string(),
            failed => failed.join(","),
        })
        .col("filter", |r: &R| {
            if r.record.quality.filter.is_empty() {
                ".".to_string()
            } else {
                r.record.quality.filter.clone()
            }
        })
        .col("region_flags", move |r: &R| {
            join(r.record.quality.flags.iter().map(|f| f.to_string()).collect())
        })
        .col("AC", |r: &R| r.record.all().ac.to_string())
        .col("AN", |r: &R| r.record.all().an.to_string())
        .col("people_one_copy", |r: &R| or_dot(r.record.quality.heterozygous, |n| n.to_string()))
        .col("people_two_copies", |r: &R| or_dot(r.record.quality.homozygous, |n| n.to_string()))
        .col("hwe_p", move |r: &R| {
            join(
                r.counted
                    .dataset
                    .groups
                    .iter()
                    .zip(&r.record.quality.hwe)
                    .filter(|(_, p)| !p.is_nan())
                    .map(|(group, p)| format!("{}:{}", group.code, num(*p)))
                    .collect(),
            )
        })
        .col("excess_het_p", |r: &R| or_dot(r.record.quality.excess_het, num))
        .col("mendelian_error_rate", |r: &R| or_dot(r.record.quality.mendel, num))
        .col("vqslod", |r: &R| or_dot(r.record.quality.vqslod, |x| format!("{x:.3}")))
}

/// One row of `datasets.tsv`: a population of the run.
pub enum DatasetRow<'a> {
    Group {
        build: &'a str,
        counted: &'a Counted,
        group: &'a Group,
        chromosomes: u32,
    },
    Study {
        build: &'a str,
        spec: &'a str,
        study: &'a Study,
    },
}

impl<'a> DatasetRow<'a> {
    /// Every population the evidence of one build holds.
    pub fn of(build: &'a str, ev: &'a Evidence) -> Vec<DatasetRow<'a>> {
        let mut out: Vec<DatasetRow> = ev
            .datasets
            .iter()
            .flat_map(|counted| {
                counted
                    .dataset
                    .groups
                    .iter()
                    .zip(&counted.full)
                    .map(move |(group, chromosomes)| DatasetRow::Group {
                        build,
                        counted,
                        group,
                        chromosomes: *chromosomes,
                    })
            })
            .collect();
        let spec = ev.sources.dbsnp.as_deref().unwrap_or(".");
        out.extend(ev.dbsnp.studies.iter().map(|study| DatasetRow::Study { build, spec, study }));
        out
    }
}

pub fn datasets_table<'a>() -> Table<DatasetRow<'a>> {
    type R<'a> = DatasetRow<'a>;
    use DatasetRow::{Group, Study};
    Table::new()
        .col("build", |r: &R| match r {
            Group { build, .. } | Study { build, .. } => build.to_string(),
        })
        .col("dataset", |r: &R| match r {
            Group { counted, .. } => counted.dataset.name.to_string(),
            Study { .. } => "dbSNP".to_string(),
        })
        .col("dataset_title", |r: &R| match r {
            Group { counted, .. } => counted.dataset.title.to_string(),
            Study { .. } => catalogue::DBSNP_TITLE.to_string(),
        })
        .col("group", |r: &R| match r {
            Group { group, .. } => group.code.to_string(),
            Study { study, .. } => study.name.clone(),
        })
        .col("group_title", |r: &R| match r {
            Group { group, .. } => group.title.to_string(),
            Study { study, .. } => study.name.clone(),
        })
        .col("ancestry", |r: &R| match r {
            Group { group, .. } => group.ancestry.as_str().to_string(),
            Study { study, .. } => study.ancestry().as_str().to_string(),
        })
        .col("counts", |r: &R| match r {
            Group { .. } => "reported".to_string(),
            Study { .. } => "estimated".to_string(),
        })
        .col("chromosomes", |r: &R| match r {
            Group { chromosomes, .. } => chromosomes.to_string(),
            Study { study, .. } => study.chromosomes.to_string(),
        })
        .col("used", |r: &R| match r {
            Group { chromosomes, .. } => if *chromosomes > 0 { "yes" } else { "no" }.to_string(),
            Study { study, .. } => if study.used { "yes" } else { "no" }.to_string(),
        })
        .col("reason", |r: &R| match r {
            Group { chromosomes, .. } => if *chromosomes > 0 {
                "."
            } else {
                "no chromosomes in the regions read"
            }
            .to_string(),
            Study { study, .. } => study.reason.to_string(),
        })
        .col("every_position_examined", |r: &R| match r {
            Group { counted, .. } => if counted.dataset.complete { "yes" } else { "no" }.to_string(),
            Study { .. } => "no".to_string(),
        })
        .col("records", |r: &R| match r {
            Group { counted, .. } => counted.sites.len().to_string(),
            Study { study, .. } => study.records.to_string(),
        })
        .col("smallest_nonzero_AF", |r: &R| match r {
            Group { .. } => ".".to_string(),
            Study { study, .. } => num(study.smallest_af),
        })
        .col("people", |r: &R| match r {
            Group { counted, .. } => or_dot(counted.dataset.people, |n| n.to_string()),
            Study { .. } => ".".to_string(),
        })
        .col("read_from", |r: &R| match r {
            Group { counted, .. } => counted.dataset.spec.clone(),
            Study { spec, .. } => spec.to_string(),
        })
}

//! Step 2: read the population data for the loci of the run, and answer one
//! question about any variant: in which populations was it looked for, and
//! how often was it seen?
//!
//! A *stratum* is one group of people of one dataset: a superpopulation of
//! 1000 Genomes, a genetic-ancestry group of gnomAD, the Japanese of ToMMo,
//! a study relayed by dbSNP. Everything downstream (class, support, score)
//! is computed from the strata of a variant and nothing else. Strata are
//! never pooled: the same people can be in several datasets, and a
//! statistic that takes the best stratum does not count them twice.

use crate::fetch::{fetch_vcf, merge_regions, Cache, Region};
use crate::hap::Panel;
use crate::locus::{Placed, Status, Window};
use crate::population::dbsnp::Dbsnp;
use crate::population::{Ancestry, Dataset, Record};
use crate::stats::PopCount;
use crate::variant::Lookup;
use crate::{io, note, tabix};
use anyhow::Result;
use std::collections::HashMap;

/// One group of people of one dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stratum {
    /// A group of a dataset with counts: positions in `Evidence::datasets`
    /// and in the groups of that dataset.
    Counted { dataset: usize, group: usize },
    /// A study relayed by dbSNP, by its position among the studies in use.
    Relayed(usize),
}

/// A stratum in words.
pub struct Described {
    /// The dataset as written in the output; `dbSNP` for a relayed study.
    pub dataset: &'static str,
    pub group: String,
    pub title: String,
    pub ancestry: Ancestry,
    /// The counts are worked out from a frequency and an estimated size.
    pub estimated: bool,
}

/// What a run reads for the alleles of one genome build.
pub struct Sources {
    /// The datasets with counts, in the order they are reported.
    pub datasets: Vec<Dataset>,
    pub dbsnp: Option<String>,
    /// Studies relayed by dbSNP that are left out by name.
    pub dbsnp_exclude: Vec<String>,
    /// The phased panel and the pedigree of its samples.
    pub phased: Option<String>,
    pub pedigree: String,
    /// The frequency a `common` variant must reach (`--common`): a study
    /// relayed by dbSNP must be large enough to show it.
    pub common: f64,
}

/// True for a whole chromosome (1-22, X, Y, with or without `chr`), false
/// for an unplaced scaffold or an alternate contig. The population call
/// sets are made for whole chromosomes: what they report elsewhere is not
/// used as evidence.
pub fn whole_chromosome(contig: &str) -> bool {
    let bare = contig.strip_prefix("chr").unwrap_or(contig);
    matches!(bare, "X" | "Y") || bare.parse::<u32>().is_ok_and(|n| (1..=22).contains(&n))
}

/// The loci of the alleles that satisfy `wanted`, merged into regions.
fn loci(windows: &[Window], placed: &[Placed], wanted: impl Fn(&Placed) -> bool) -> Vec<Region> {
    let regions = windows
        .iter()
        .filter(|w| whole_chromosome(&w.coord.contig))
        .filter(|w| placed.iter().any(|p| p.gene == w.coord.gene && wanted(p)))
        .map(|w| (w.coord.contig.clone(), w.start, w.end()))
        .collect();
    merge_regions(regions, 10_000)
}

/// Frequencies are read for loci with at least one variant.
fn frequency_regions(windows: &[Window], placed: &[Placed]) -> Vec<Region> {
    loci(windows, placed, |p| !p.variants.is_empty())
}

/// The phased panel is read for loci with any placed allele.
fn haplotype_regions(windows: &[Window], placed: &[Placed]) -> Vec<Region> {
    loci(windows, placed, |p| p.placement.is_some() && p.status != Status::LowIdentity)
}

/// What the two sources that are not datasets with counts are called in
/// progress messages.
const DBSNP: &str = "dbSNP";
const PHASED: &str = "the phased chromosomes of 1000 Genomes";

/// One tabix-indexed file a run reads, and the regions it reads of it.
pub struct Download {
    /// What the file is, in words.
    pub title: String,
    /// Path or URL; `{contig}` stands for the contig name.
    pub spec: String,
    pub regions: Vec<Region>,
}

/// Everything `Evidence::load` will read region by region, so that a
/// program able to download concurrently can fetch it ahead.
pub fn downloads(sources: &Sources, windows: &[Window], placed: &[Placed]) -> Vec<Download> {
    let frequency = frequency_regions(windows, placed);
    let mut out: Vec<Download> = sources
        .datasets
        .iter()
        .map(|dataset| (dataset.title, dataset.spec.clone()))
        .chain(sources.dbsnp.clone().map(|spec| (DBSNP, spec)))
        .map(|(title, spec)| Download {
            title: title.to_string(),
            spec,
            regions: frequency.clone(),
        })
        .collect();
    if let Some(spec) = &sources.phased {
        out.push(Download {
            title: PHASED.to_string(),
            spec: spec.clone(),
            regions: haplotype_regions(windows, placed),
        });
    }
    out.retain(|d| !d.regions.is_empty());
    out
}

/// A dataset as read for the loci of a run.
pub struct Counted {
    pub dataset: Dataset,
    pub sites: HashMap<Lookup, Record>,
    /// Per group, the chromosomes of a fully called site: the largest
    /// number examined at any record of the regions read.
    pub full: Vec<u32>,
    /// The same over all people of the dataset.
    pub full_all: u32,
}

impl Counted {
    fn read(dataset: Dataset, regions: &[Region], cache: Cache) -> Result<Self> {
        let mut sites = HashMap::new();
        for (contig, line) in fetch_vcf(&dataset.spec, regions, cache)? {
            for (key, record) in dataset.layout.parse(&line) {
                sites.insert((contig.clone(), key), record);
            }
        }
        let mut full = vec![0; dataset.groups.len()];
        let mut full_all = 0;
        for record in sites.values() {
            full_all = full_all.max(record.all().an);
            for (most, count) in full.iter_mut().zip(&record.counts) {
                *most = (*most).max(count.an);
            }
        }
        Ok(Counted {
            dataset,
            sites,
            full,
            full_all,
        })
    }

    /// The site criteria the record of this variant fails.
    pub fn failures(&self, lookup: &Lookup) -> Vec<String> {
        self.sites.get(lookup).map_or(Vec::new(), |record| {
            record
                .quality
                .failures(record.all().an, self.dataset.complete.then_some(self.full_all))
        })
    }
}

pub struct Evidence {
    pub sources: Sources,
    /// The regions whose population data was read.
    regions: Vec<Region>,
    /// The datasets that were read, in the order of `sources.datasets`.
    pub datasets: Vec<Counted>,
    pub dbsnp: Dbsnp,
    pub panel: Option<Panel>,
}

impl Evidence {
    /// Read every source for the loci that need it. Frequencies: loci with
    /// at least one variant. Phased panel: loci with any placed allele.
    pub fn load(mut sources: Sources, windows: &[Window], placed: &[Placed], cache: Cache) -> Result<Self> {
        let regions = frequency_regions(windows, placed);
        if regions.is_empty() {
            sources.datasets.clear();
            sources.dbsnp = None;
        }
        let mut datasets = Vec::new();
        for dataset in std::mem::take(&mut sources.datasets) {
            note!("reading {}: {}", dataset.title, dataset.spec);
            let (required, title) = (dataset.required, dataset.title);
            match Counted::read(dataset, &regions, cache) {
                Ok(counted) => datasets.push(counted),
                Err(e) if required => return Err(e),
                Err(e) => note!("warning: {title} could not be read and is left out: {e:#}"),
            }
        }
        sources.datasets = datasets.iter().map(|counted| counted.dataset.clone()).collect();
        let mut dbsnp = Dbsnp::default();
        if let Some(spec) = &sources.dbsnp {
            note!("reading {DBSNP}: {spec}");
            // If dbSNP is unavailable the run continues without it.
            let lines = fetch_vcf(spec, &regions, cache).unwrap_or_else(|e| {
                note!("warning: dbSNP could not be read; no rsIDs, no relayed frequencies: {e:#}");
                Vec::new()
            });
            for (contig, line) in lines {
                dbsnp.add_line(&contig, &line);
            }
            // A study that reports people read directly is left out.
            let duplicated: Vec<&str> = sources
                .datasets
                .iter()
                .flat_map(|dataset| dataset.replaces.iter().copied())
                .collect();
            dbsnp.select(&sources.dbsnp_exclude, &duplicated, sources.common);
        }
        let hap_regions = haplotype_regions(windows, placed);
        if hap_regions.is_empty() {
            sources.phased = None;
        }
        let mut panel = None;
        if let Some(spec) = &sources.phased {
            note!("reading {PHASED}: {spec}");
            // Additional evidence: if the panel is unavailable the run
            // continues without the haplotype check.
            match load_panel(spec, &sources.pedigree, &hap_regions, cache) {
                Ok(loaded) => {
                    note!("  {} founders, {} haplotypes", loaded.founders, 2 * loaded.founders);
                    panel = Some(loaded);
                }
                Err(e) => note!("warning: haplotype check skipped: {e:#}"),
            }
        }
        Ok(Evidence {
            sources,
            regions,
            datasets,
            dbsnp,
            panel,
        })
    }

    /// True if at least one source with allele counts was consulted.
    pub fn counted(&self) -> bool {
        self.datasets.iter().any(|counted| counted.full_all > 0) || self.dbsnp.used().next().is_some()
    }

    /// True if the population data of this run cover the contig: a variant
    /// without a record there was looked for and is absent.
    pub fn reads(&self, contig: &str) -> bool {
        self.counted() && self.regions.iter().any(|region| region.0 == contig)
    }

    /// True if population data was read for the whole window.
    pub fn covers(&self, window: &Window) -> bool {
        self.regions
            .iter()
            .any(|(contig, from, to)| contig == &window.coord.contig && *from <= window.start && window.end() <= *to)
    }

    /// Every stratum that looked at the variant's position, with its
    /// counts: the groups of a dataset that has a record of the variant;
    /// the groups of a dataset that examines everyone everywhere and has
    /// no record, with a count of zero (the variant was looked for and not
    /// seen); and the studies dbSNP relays a frequency from.
    pub fn strata(&self, lookup: &Lookup) -> Vec<(Stratum, PopCount)> {
        let mut out = Vec::new();
        let read = self.reads(&lookup.0);
        for (dataset, counted) in self.datasets.iter().enumerate() {
            let stratum = |group| Stratum::Counted { dataset, group };
            match counted.sites.get(lookup) {
                Some(record) => out.extend(record.counts.iter().enumerate().map(|(group, count)| (stratum(group), *count))),
                None if read && counted.dataset.complete => out.extend(
                    counted
                        .full
                        .iter()
                        .enumerate()
                        .map(|(group, an)| (stratum(group), PopCount::new(0, *an))),
                ),
                None => {}
            }
        }
        out.extend(self.dbsnp.counts(lookup).map(|(study, count)| (Stratum::Relayed(study), count)));
        out.retain(|(_, count)| count.an > 0);
        out
    }

    /// The number of strata the run can test a variant in: every group
    /// with chromosomes in the regions read, and every relayed study in
    /// use. A test repeated over the strata is corrected for this number,
    /// which is fixed before any variant is looked at.
    pub fn strata_in_use(&self) -> usize {
        self.datasets
            .iter()
            .map(|counted| counted.full.iter().filter(|an| **an > 0).count())
            .sum::<usize>()
            + self.dbsnp.used().count()
    }

    /// Who a stratum is.
    pub fn describe(&self, stratum: Stratum) -> Described {
        match stratum {
            Stratum::Counted { dataset, group } => {
                let dataset = &self.datasets[dataset].dataset;
                let group = &dataset.groups[group];
                Described {
                    dataset: dataset.name,
                    group: group.code.to_string(),
                    title: group.title.to_string(),
                    ancestry: group.ancestry,
                    estimated: false,
                }
            }
            Stratum::Relayed(study) => {
                let study = self.dbsnp.used().nth(study);
                let name = study.map_or("?", |s| s.name.as_str()).to_string();
                Described {
                    dataset: "dbSNP",
                    group: name.clone(),
                    title: name,
                    ancestry: study.map_or(Ancestry::Other, |s| s.ancestry()),
                    estimated: true,
                }
            }
        }
    }

    /// A stratum as `dataset:group`.
    pub fn stratum_name(&self, stratum: Stratum) -> String {
        let described = self.describe(stratum);
        format!("{}:{}", described.dataset, described.group)
    }

    /// The sources in which the variant was seen at least once, by name.
    /// A record whose count is zero after a dataset's own genotype
    /// filtering is not an observation.
    pub fn observed(&self, lookup: &Lookup) -> Vec<&'static str> {
        let mut out: Vec<&'static str> = self
            .datasets
            .iter()
            .filter(|counted| counted.sites.get(lookup).is_some_and(|record| record.all().ac > 0))
            .map(|counted| counted.dataset.name)
            .collect();
        if !self.dbsnp.relayed(lookup).is_empty() {
            out.push("dbSNP");
        }
        out
    }

    /// The site criteria the variant's position fails, per dataset that
    /// has any to report.
    pub fn failures(&self, lookup: &Lookup) -> Vec<(&'static str, Vec<String>)> {
        self.datasets
            .iter()
            .map(|counted| (counted.dataset.name, counted.failures(lookup)))
            .filter(|(_, failed)| !failed.is_empty())
            .collect()
    }
}

/// Support of a variant: the largest lower confidence bound of its
/// frequency over all strata; 0 if it was never seen.
pub fn support(strata: &[(Stratum, PopCount)]) -> f64 {
    strata.iter().map(|(_, c)| c.lower95()).fold(0.0, f64::max)
}

fn load_panel(spec: &str, pedigree: &str, regions: &[Region], cache: Cache) -> Result<Panel> {
    let file = spec.replace("{contig}", &regions[0].0);
    let header = cache.text(&file, "samples", || Ok(vec![tabix::sample_line(&file)?]))?;
    let pedigree = cache.text(&file, "pedigree", || {
        Ok(String::from_utf8_lossy(&io::read_whole(pedigree)?)
            .lines()
            .map(str::to_string)
            .collect())
    })?;
    let mut panel = Panel::new(header.first().map(String::as_str).unwrap_or(""), &pedigree.join("\n"))?;
    for (contig, line) in fetch_vcf(spec, regions, cache)? {
        panel.add_line(&contig, &line);
    }
    Ok(panel)
}

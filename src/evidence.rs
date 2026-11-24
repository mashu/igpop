//! Step 2: read the population data for the loci of the run, and answer one
//! question about any variant: in which populations was it seen, how often?
//!
//! A *stratum* is one population of one dataset: a superpopulation of 1000
//! Genomes, a genetic-ancestry group of gnomAD, or a dataset relayed by
//! dbSNP. Everything downstream (class, support, score) is computed from
//! the strata of a variant and nothing else.

use crate::fetch::{fetch_vcf, merge_regions, Cache, Region};
use crate::freq::{self, GnomadSite, Key, KgpRelease, PopCount, Site, GNOMAD_GROUPS, POPS};
use crate::hap::Panel;
use crate::locus::{Placed, Status, Window};
use crate::{io, note, tabix};
use anyhow::Result;
use std::collections::{HashMap, HashSet};

/// A variant as (contig, (position, REF, ALT)).
pub type Lookup = (String, Key);

/// One population of one dataset, by its position in that dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stratum {
    /// Superpopulation of 1000 Genomes (index into `POPS`).
    Kgp(usize),
    /// Genetic-ancestry group of gnomAD (index into `GNOMAD_GROUPS`).
    Gnomad(usize),
    /// Dataset relayed by dbSNP (index into the datasets in use).
    Dbsnp(usize),
}

/// Where each source is read from; `None` = not used in this run.
pub struct Sources {
    pub kgp: Option<String>,
    /// Which release `kgp` is, and what to call it in the output.
    pub kgp_release: KgpRelease,
    pub kgp_name: &'static str,
    pub gnomad: Option<String>,
    pub gnomad_name: &'static str,
    pub dbsnp: Option<String>,
    pub phased: Option<String>,
    pub pedigree: String,
    /// dbSNP datasets left out by name.
    pub dbsnp_exclude: Vec<String>,
    /// The frequency a `common` variant must reach (`--common`).
    pub common: f64,
}

/// True if a stratum of `an` chromosomes is large enough to tell a common
/// variant from a rare one: not seeing the allele at all would bound its
/// frequency below `common`.
pub fn informative(an: u32, common: f64) -> bool {
    unseen(an).upper95() < common
}

fn unseen(an: u32) -> PopCount {
    PopCount { ac: 0, an, het: 0, hom: 0 }
}

/// Largest upper bound of the frequency of an allele seen in none of the
/// strata, whose sizes are `full`.
pub fn absent_upper(full: &[u32]) -> f64 {
    full.iter().map(|an| unseen(*an).upper95()).fold(f64::NAN, f64::max)
}

/// One frequency dataset relayed by dbSNP, and whether it is used.
pub struct Dataset {
    pub name: String,
    pub records: u64,
    pub smallest_af: f64,
    /// Estimated size: one copy among all chromosomes gives the smallest
    /// non-zero frequency a dataset can report.
    pub chromosomes: u32,
    pub used: bool,
    pub reason: &'static str,
}

/// What dbSNP holds for the regions read.
#[derive(Default)]
pub struct Dbsnp {
    pub rsids: HashMap<Lookup, String>,
    /// Frequencies per variant allele, as (dataset, frequency).
    freqs: HashMap<Lookup, Vec<(String, f64)>>,
    /// Alleles carrying the COMMON flag: of a site with several alternate
    /// alleles, only the most frequent one.
    pub common: HashSet<Lookup>,
    /// Per dataset: number of records, smallest non-zero frequency.
    seen: HashMap<String, (u64, f64)>,
    /// Decided once all lines are read, largest dataset first.
    pub datasets: Vec<Dataset>,
    /// The datasets in use, by name: positions in `datasets`.
    used: Vec<usize>,
}

impl Dbsnp {
    fn add_line(&mut self, contig: &str, line: &str) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 || f[2] == "." || f[4] == "." {
            return;
        }
        let Ok(pos) = f[1].parse::<u64>() else { return };
        let reported = f.get(7).map(|info| freq::parse_dbsnp_freq(info)).unwrap_or_default();
        for (dataset, values) in &reported {
            let entry = self.seen.entry(dataset.clone()).or_insert((0, f64::INFINITY));
            entry.0 += 1;
            for v in values.iter().skip(1).flatten().filter(|v| **v > 0.0) {
                entry.1 = entry.1.min(*v);
            }
        }
        let alts: Vec<&str> = f[4].split(',').collect();
        let lookup = |alt: &str| (contig.to_string(), freq::trim(pos, f[3].as_bytes(), alt.as_bytes()));
        if f.get(7).is_some_and(|info| info.split(';').any(|x| x == "COMMON")) {
            // Frequency of each ALT: from 1000 Genomes if reported there
            // (the flag is defined on it), otherwise the highest reported.
            let of_alt = |i: usize| -> f64 {
                let pick = |name: &str| {
                    reported
                        .iter()
                        .find(|(d, _)| d == name)
                        .and_then(|(_, v)| v.get(i + 1).copied().flatten())
                };
                pick("1000Genomes").or(pick("1000Genomes_30X")).unwrap_or_else(|| {
                    reported
                        .iter()
                        .filter_map(|(_, v)| v.get(i + 1).copied().flatten())
                        .fold(0.0, f64::max)
                })
            };
            let top = (0..alts.len())
                .max_by(|a, b| of_alt(*a).total_cmp(&of_alt(*b)).then(b.cmp(a)))
                .unwrap_or(0);
            self.common.insert(lookup(alts[top]));
        }
        for (i, alt) in alts.iter().enumerate() {
            let of_alt: Vec<(String, f64)> = reported
                .iter()
                .filter_map(|(dataset, values)| values.get(i + 1).copied().flatten().map(|v| (dataset.clone(), v)))
                .collect();
            self.rsids.insert(lookup(alt), f[2].to_string());
            self.freqs.insert(lookup(alt), of_alt);
        }
    }

    /// Decide which datasets are used. `duplicated`: datasets that re-report
    /// a source read directly in this run.
    fn select(&mut self, excluded: &[String], duplicated: &[&str], common: f64) {
        self.datasets = self
            .seen
            .iter()
            .map(|(name, (records, smallest_af))| {
                let chromosomes = if smallest_af.is_finite() {
                    (1.0 / smallest_af).round() as u32
                } else {
                    0
                };
                let reason = if excluded.contains(name) {
                    "excluded by --dbsnp-exclude"
                } else if duplicated.contains(&name.as_str()) {
                    "re-report of a source read directly"
                } else if !informative(chromosomes, common) {
                    "too small for --common"
                } else {
                    "."
                };
                Dataset {
                    name: name.clone(),
                    records: *records,
                    smallest_af: *smallest_af,
                    chromosomes,
                    used: reason == ".",
                    reason,
                }
            })
            .collect();
        self.datasets
            .sort_by(|a, b| a.smallest_af.total_cmp(&b.smallest_af).then(a.name.cmp(&b.name)));
        self.used = (0..self.datasets.len()).filter(|i| self.datasets[*i].used).collect();
        self.used.sort_by(|a, b| self.datasets[*a].name.cmp(&self.datasets[*b].name));
    }

    /// The datasets in use, in the order `Stratum::Dbsnp` numbers them.
    pub fn used(&self) -> impl Iterator<Item = &Dataset> {
        self.used.iter().map(|i| &self.datasets[*i])
    }

    /// Non-zero frequencies of the used datasets for one allele, highest first.
    pub fn relayed(&self, lookup: &Lookup) -> Vec<(String, f64)> {
        let mut out: Vec<(String, f64)> = self
            .freqs
            .get(lookup)
            .map(|all| {
                all.iter()
                    .filter(|(d, af)| *af > 0.0 && self.used().any(|u| &u.name == d))
                    .cloned()
                    .collect()
            })
            .unwrap_or_default();
        out.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        out
    }

    /// One stratum per used dataset that reports this allele. dbSNP relays
    /// a frequency without counts; the count is frequency x estimated size.
    fn strata<'a>(&'a self, lookup: &Lookup) -> impl Iterator<Item = (Stratum, PopCount)> + 'a {
        let reported = self.freqs.get(lookup);
        self.used().enumerate().filter_map(move |(k, d)| {
            let (_, af) = reported?.iter().find(|(name, _)| name == &d.name)?;
            let ac = (af * d.chromosomes as f64).round() as u32;
            let count = PopCount {
                ac,
                an: d.chromosomes,
                het: 0,
                hom: 0,
            };
            Some((Stratum::Dbsnp(k), count))
        })
    }
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

/// One tabix-indexed file a run reads, and the regions it reads of it.
pub struct Download {
    /// Path or URL; `{contig}` stands for the contig name.
    pub spec: String,
    pub regions: Vec<Region>,
}

/// Everything `Evidence::load` will read region by region, so that a
/// program able to download concurrently can fetch it ahead.
pub fn downloads(sources: &Sources, windows: &[Window], placed: &[Placed]) -> Vec<Download> {
    let frequency = frequency_regions(windows, placed);
    let mut out: Vec<Download> = [&sources.kgp, &sources.gnomad, &sources.dbsnp]
        .into_iter()
        .flatten()
        .map(|spec| Download {
            spec: spec.clone(),
            regions: frequency.clone(),
        })
        .collect();
    if let Some(spec) = &sources.phased {
        out.push(Download {
            spec: spec.clone(),
            regions: haplotype_regions(windows, placed),
        });
    }
    out.retain(|d| !d.regions.is_empty());
    out
}

pub struct Evidence {
    pub sources: Sources,
    /// The regions whose population data was read.
    regions: Vec<Region>,
    pub kgp: HashMap<Lookup, Site>,
    /// Chromosomes genotyped at a fully called site, per superpopulation:
    /// the denominator for a variant without a record.
    pub kgp_full: [u32; 5],
    pub gnomad: HashMap<Lookup, GnomadSite>,
    pub gnomad_full: [u32; 9],
    pub gnomad_full_total: u32,
    pub dbsnp: Dbsnp,
    pub panel: Option<Panel>,
}

impl Evidence {
    /// Read every source for the loci that need it. Frequencies: loci with
    /// at least one variant. Phased panel: loci with any placed allele.
    pub fn load(sources: Sources, windows: &[Window], placed: &[Placed], cache: Cache) -> Result<Self> {
        let regions = frequency_regions(windows, placed);
        let mut ev = Evidence {
            sources,
            regions,
            kgp: HashMap::new(),
            kgp_full: [0; 5],
            gnomad: HashMap::new(),
            gnomad_full: [0; 9],
            gnomad_full_total: 0,
            dbsnp: Dbsnp::default(),
            panel: None,
        };
        if ev.regions.is_empty() {
            (ev.sources.kgp, ev.sources.gnomad, ev.sources.dbsnp) = (None, None, None);
        }

        if let Some(spec) = &ev.sources.kgp {
            note!("{}: {spec}", ev.sources.kgp_name);
            for (contig, line) in fetch_vcf(spec, &ev.regions, cache)? {
                for (key, site) in ev.sources.kgp_release.parse(&line) {
                    ev.kgp.insert((contig.clone(), key), site);
                }
            }
            for site in ev.kgp.values() {
                for (full, pop) in ev.kgp_full.iter_mut().zip(&site.pops) {
                    *full = (*full).max(pop.an);
                }
            }
        }
        if let Some(spec) = &ev.sources.gnomad {
            note!("{}: {spec}", ev.sources.gnomad_name);
            for (contig, line) in fetch_vcf(spec, &ev.regions, cache)? {
                for (key, site) in freq::parse_gnomad_line(&line) {
                    ev.gnomad.insert((contig.clone(), key), site);
                }
            }
            for site in ev.gnomad.values() {
                ev.gnomad_full_total = ev.gnomad_full_total.max(site.total.an);
                for (full, group) in ev.gnomad_full.iter_mut().zip(&site.groups) {
                    *full = (*full).max(group.an);
                }
            }
        }
        if let Some(spec) = &ev.sources.dbsnp {
            note!("dbSNP: {spec}");
            // If dbSNP is unavailable the run continues without it.
            let lines = fetch_vcf(spec, &ev.regions, cache).unwrap_or_else(|e| {
                note!("warning: dbSNP could not be read; no rsIDs, no relayed frequencies: {e:#}");
                Vec::new()
            });
            for (contig, line) in lines {
                ev.dbsnp.add_line(&contig, &line);
            }
            let mut duplicated: Vec<&str> = Vec::new();
            if ev.sources.kgp.is_some() {
                duplicated.extend(["1000Genomes", "1000Genomes_30X"]);
            }
            if ev.sources.gnomad.is_some() {
                duplicated.push("GnomAD_genomes");
            }
            ev.dbsnp.select(&ev.sources.dbsnp_exclude, &duplicated, ev.sources.common);
        }

        let hap_regions = haplotype_regions(windows, placed);
        if hap_regions.is_empty() {
            ev.sources.phased = None;
        }
        if let Some(spec) = &ev.sources.phased {
            note!("phased panel: {spec}");
            // Additional evidence: if the panel is unavailable the run
            // continues without the haplotype check.
            match load_panel(spec, &ev.sources.pedigree, &hap_regions, cache) {
                Ok(panel) => {
                    note!("  {} founders, {} haplotypes", panel.founders, 2 * panel.founders);
                    ev.panel = Some(panel);
                }
                Err(e) => note!("warning: haplotype check skipped: {e:#}"),
            }
        }
        Ok(ev)
    }

    pub fn kgp_full_total(&self) -> u32 {
        self.kgp_full.iter().sum()
    }

    /// True if the 1000 Genomes / gnomAD data cover the loci of the run, so
    /// that a variant without a record there counts as "looked for, absent".
    pub fn kgp_read(&self) -> bool {
        self.sources.kgp.is_some() && self.kgp_full_total() > 0
    }

    pub fn gnomad_read(&self) -> bool {
        self.sources.gnomad.is_some() && self.gnomad_full_total > 0
    }

    /// True if at least one source with allele counts was consulted.
    pub fn counted(&self) -> bool {
        self.kgp_read() || self.gnomad_read() || self.dbsnp.used().next().is_some()
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

    /// Every stratum that genotyped the variant's position, with its counts.
    pub fn strata(&self, lookup: &Lookup) -> Vec<(Stratum, PopCount)> {
        let mut out = Vec::new();
        if let Some(site) = self.kgp.get(lookup) {
            out.extend(site.pops.iter().enumerate().map(|(i, c)| (Stratum::Kgp(i), *c)));
        }
        if let Some(site) = self.gnomad.get(lookup) {
            out.extend(site.groups.iter().enumerate().map(|(i, c)| (Stratum::Gnomad(i), *c)));
        }
        out.extend(self.dbsnp.strata(lookup));
        out.retain(|(_, c)| c.an > 0);
        out
    }

    /// A stratum as `dataset:population`.
    pub fn stratum_name(&self, stratum: Stratum) -> String {
        match stratum {
            Stratum::Kgp(i) => format!("{}:{}", self.sources.kgp_name, POPS[i]),
            Stratum::Gnomad(i) => format!("{}:{}", self.sources.gnomad_name, GNOMAD_GROUPS[i]),
            Stratum::Dbsnp(k) => format!("dbSNP:{}", self.dbsnp.used().nth(k).map_or("?", |d| d.name.as_str())),
        }
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

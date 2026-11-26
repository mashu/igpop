//! One run: place the alleles, read the population data, judge, and write
//! the tables and the report.

use crate::assess::{self, AlleleHaplotype, AlleleResult, HaplotypeCheck, ScoreModel, Yardstick};
use crate::config::Config;
use crate::evidence::{self, Evidence, Sources};
use crate::locus::{self, Placed, Status};
use crate::plan::{BuildData, Plan, BY_SEQUENCE};
use crate::{align, fetch, io, net, note, par, report, tables};
use anyhow::{bail, Result};
use std::collections::{HashMap, HashSet};

/// The result of a run.
pub struct Output {
    /// (file name without the prefix, content), the report first.
    pub files: Vec<(&'static str, String)>,
    /// A few lines saying what happened.
    pub summary: Vec<String>,
    /// What limits how far the results can be relied on.
    pub warnings: Vec<String>,
}

/// The alleles judged on one genome build, with its loci and data.
struct Stage {
    build: String,
    windows: Vec<locus::Window>,
    placed: Vec<Placed>,
    sources: Sources,
}

/// A run in two halves, for programs that can download concurrently:
/// after `start` (loci and placement) they may fetch `downloads()` ahead,
/// before `finish` reads the population data and judges the alleles.
pub struct Session<'a> {
    config: &'a Config,
    coords_label: String,
    /// Genes that could not be placed because their coordinates are on a
    /// build that is not used in this run.
    elsewhere: HashMap<String, String>,
    /// Alleles placed by sequence because their gene is not in the
    /// coordinates table, with the gene their name gives.
    by_sequence: HashMap<String, String>,
    /// The build that was asked for first, then any other build that holds
    /// genes the first one lacks.
    stages: Vec<Stage>,
}

/// One whole run.
pub fn run(config: &Config) -> Result<Output> {
    Session::start(config)?.finish()
}

impl<'a> Session<'a> {
    /// Step 1: loci and placement.
    pub fn start(config: &'a Config) -> Result<Self> {
        config.check()?;
        let (plan, mut aliases) = Plan::new(config)?;
        let Plan {
            builds,
            all_coords,
            coords_label,
            elsewhere,
        } = plan;
        let options = locus::PlaceOptions {
            min_identity: config.min_identity,
            paralog_check: config.paralog_check,
            scoring: align::Scoring::default(),
        };
        let alleles = io::parse_fasta(&config.alleles.content);
        if alleles.is_empty() {
            bail!("{} holds no FASTA records", config.alleles.name);
        }
        let on_main: HashSet<&str> = all_coords
            .iter()
            .filter(|c| c.build == builds[0])
            .map(|c| c.gene.as_str())
            .collect();
        // The genes of each build: all of the main build, and of any other
        // build those the main build lacks.
        let coords_of = |n: usize| -> Vec<io::Coord> {
            all_coords
                .iter()
                .filter(|c| c.build == builds[n] && (n == 0 || !on_main.contains(c.gene.as_str())))
                .cloned()
                .collect()
        };
        // An allele belongs to the first build that has its gene; one that
        // no build has stays with the main build.
        let home = |name: &str, aliases: &HashMap<String, String>| {
            (0..builds.len()).find(|n| coords_of(*n).iter().any(|c| c.gene == locus::gene_for(name, aliases)))
        };
        let homes: Vec<Option<usize>> = alleles.iter().map(|(name, _)| home(name, &aliases)).collect();
        let mut by_sequence: HashMap<String, String> = HashMap::new();

        let mut stages = Vec::new();
        for (n, build) in builds.iter().enumerate() {
            let main = n == 0;
            let coords = coords_of(n);
            let members: Vec<(String, Vec<u8>)> = alleles
                .iter()
                .zip(&homes)
                .filter(|(_, home)| home.unwrap_or(0) == n)
                .map(|(allele, _)| allele.clone())
                .collect();
            if members.is_empty() && !main {
                continue;
            }
            let data = BuildData::new(config, build, main);
            let Some(reference) = data.reference else {
                bail!("no default reference for build '{build}'; give one");
            };
            note!("{build} reference: {reference}");
            let windows = locus::load_windows(coords, &mut io::Reference::open(&reference)?, config.flank)?;
            if main {
                // An allele whose gene is not in the table goes to the
                // locus its sequence fits best.
                let unknown: Vec<&(String, Vec<u8>)> = alleles
                    .iter()
                    .zip(&homes)
                    .filter(|((name, _), home)| home.is_none() && (config.place_unknown || locus::gene_for(name, &aliases) == BY_SEQUENCE))
                    .map(|(allele, _)| allele)
                    .collect();
                if !unknown.is_empty() {
                    note!("finding the gene of {} alleles by sequence", unknown.len());
                }
                let found = par::map(&unknown, None, |(_, seq)| {
                    locus::best_locus(seq, &windows, &options.scoring).map(|(gene, _)| gene.to_string())
                });
                for (name, _) in &unknown {
                    by_sequence.insert(name.clone(), locus::gene_of(name).to_string());
                }
                for ((name, _), gene) in unknown.iter().zip(found) {
                    match gene {
                        Some(gene) => aliases.insert(name.clone(), gene),
                        None => by_sequence.remove(name),
                    };
                }
            }
            note!("aligning {} alleles to {} loci on {build}", members.len(), windows.len());
            let mut placed = locus::place_all(&members, &windows, &aliases, &options)?;
            placed.sort_by(|a, b| a.allele.cmp(&b.allele));
            stages.push(Stage {
                build: build.clone(),
                windows,
                placed,
                sources: data.sources,
            });
        }
        Ok(Session {
            config,
            coords_label,
            elsewhere,
            by_sequence,
            stages,
        })
    }

    /// The tabix-indexed files `finish` will read, with their regions.
    pub fn downloads(&self) -> Vec<evidence::Download> {
        let main = self.stages.first().map(|stage| stage.build.as_str());
        self.stages
            .iter()
            .flat_map(|stage| {
                // Two builds can read a source of the same name.
                let other = (Some(stage.build.as_str()) != main).then_some(&stage.build);
                evidence::downloads(&stage.sources, &stage.windows, &stage.placed)
                    .into_iter()
                    .map(move |download| evidence::Download {
                        title: match other {
                            Some(build) => format!("{} ({build})", download.title),
                            None => download.title,
                        },
                        ..download
                    })
            })
            .collect()
    }

    /// Steps 2 to 4: population evidence, assessment, tables and report.
    pub fn finish(self) -> Result<Output> {
        let Session {
            config,
            coords_label,
            elsewhere,
            by_sequence,
            stages,
        } = self;
        let cache = fetch::Cache(config.cache_dir.as_deref());
        // Step 2: the population data of each build, and its haplotype check.
        struct Loaded {
            build: String,
            windows: Vec<locus::Window>,
            placed: Vec<Placed>,
            ev: Evidence,
        }
        let mut loaded = Vec::new();
        for stage in stages {
            let ev = Evidence::load(stage.sources, &stage.windows, &stage.placed, cache)?;
            loaded.push(Loaded {
                build: stage.build,
                windows: stage.windows,
                placed: stage.placed,
                ev,
            });
        }
        let checks: Vec<HaplotypeCheck> = loaded.iter().map(|l| HaplotypeCheck::new(&l.ev, &l.windows, &l.placed)).collect();
        let main = &loaded[0];
        let (build, ev, windows, placed, haplotypes) = (&main.build, &main.ev, &main.windows, &main.placed, &checks[0]);

        // Step 3: calibration and assessment.
        let trusted_names: Option<HashSet<&str>> = config.trusted.as_ref().map(|text| {
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .collect()
        });
        let is_trusted = |p: &Placed| match &trusted_names {
            Some(names) => names.contains(p.allele.as_str()),
            None => !p.allele.contains(&config.candidate_marker),
        };
        // The haplotype check comes first: it also tells which alleles are
        // real beyond the named trusted ones.
        let checked: Vec<Vec<AlleleHaplotype>> = loaded
            .iter()
            .zip(&checks)
            .map(|(l, check)| {
                l.placed
                    .iter()
                    .enumerate()
                    .map(|(i, p)| check.of(i, p.variants.len(), config.limits.coincidence))
                    .collect()
            })
            .collect();
        // The yardstick for the score: alleles of the main build, cleanly
        // placed where every source was read, that are named as trusted or
        // whose two or more differences are established together as a
        // haplotype. Other builds have other data and do not calibrate.
        let yardstick: Vec<Yardstick> = placed
            .iter()
            .zip(&checked[0])
            .map(|(p, haplotype)| match p.status {
                Status::Ok if !ev.reads(p.contig()) => Yardstick::No,
                Status::Ok if is_trusted(p) => Yardstick::Named,
                Status::Ok if haplotype.support.is_some() && p.variants.len() >= 2 => Yardstick::Haplotype,
                _ => Yardstick::No,
            })
            .collect();
        let calibrating: Vec<(&Placed, Yardstick)> = placed
            .iter()
            .zip(&yardstick)
            .filter(|(_, y)| **y != Yardstick::No)
            .map(|(p, y)| (p, *y))
            .collect();
        let calibration = config.calibration.as_ref().map(|t| (t.content.as_str(), t.name.as_str()));
        let model = ScoreModel::new(ev, windows, &calibrating, calibration)?;
        if model.calibration.is_some() {
            note!("score: calibrated from {}", model.description.replace('`', ""));
        }
        // Every allele is judged with the data of its build and the score
        // of the main build.
        let mut results: Vec<(usize, AlleleResult)> = Vec::new();
        for (n, (l, haplotypes)) in loaded.iter().zip(checked).enumerate() {
            let judge = assess::Judge {
                ev: &l.ev,
                model: &model,
                limits: &config.limits,
            };
            for (i, (p, haplotype)) in l.placed.iter().zip(haplotypes).enumerate() {
                let marked = if n == 0 { yardstick[i] } else { Yardstick::No };
                let other_build = elsewhere.get(&p.gene).map(String::as_str);
                let mut result = assess::assess_allele(&judge, p, &l.build, is_trusted(p), marked, haplotype, other_build);
                if n > 0 {
                    result.notes.insert(
                        0,
                        format!(
                            "gene absent from {build}: judged on {} with {}",
                            l.build,
                            l.ev.datasets
                                .iter()
                                .map(|counted| counted.dataset.name)
                                .collect::<Vec<_>>()
                                .join(" and ")
                        ),
                    );
                } else if p.placement.is_some() && !evidence::whole_chromosome(p.contig()) {
                    result.notes.insert(
                        0,
                        format!(
                            "gene on {}, which is not part of a chromosome: no population data are read for it",
                            p.contig()
                        ),
                    );
                }
                if let Some(named) = by_sequence.get(&p.allele) {
                    result.notes.insert(
                        0,
                        format!(
                            "gene {named} is not in the gene list: placed on {} because its sequence fits that gene best",
                            p.gene
                        ),
                    );
                }
                results.push((n, result));
            }
        }
        results.sort_by(|a, b| a.1.placed.allele.cmp(&b.1.placed.allele));
        let stage_of: Vec<usize> = results.iter().map(|(n, _)| *n).collect();
        let results: Vec<AlleleResult> = results.into_iter().map(|(_, r)| r).collect();

        // Step 4: tables and report.
        let variant_rows: Vec<tables::VariantRow> = results
            .iter()
            .zip(&stage_of)
            .flat_map(|(allele, n)| {
                let ev = &loaded[*n].ev;
                allele
                    .variants
                    .iter()
                    .zip(&allele.placed.variants)
                    .map(move |(v, called)| tables::VariantRow { ev, allele, called, v })
            })
            .collect();
        let placements: Vec<tables::PlacedRow> = results
            .iter()
            .map(|r| tables::PlacedRow {
                build: r.build,
                placed: r.placed,
            })
            .collect();
        let mut alignments = tables::alignments_legend();
        for l in &loaded {
            alignments.push_str(&tables::alignments_text(&l.build, &l.windows, &l.placed)?);
        }
        let (requests, bytes) = net::totals();
        let summary_input = report::Input {
            prefix: &config.prefix,
            build,
            other_builds: loaded[1..].iter().map(|l| (l.build.as_str(), &l.ev, l.placed.len())).collect(),
            alleles: config.alleles.name.clone(),
            coords: coords_label,
            limits: &config.limits,
            results: &results,
            ev,
            model: &model,
            haplotypes,
            requests,
            bytes,
        };
        // The population data behind the verdicts: each difference once,
        // however many alleles carry it.
        let mut listed: HashSet<(usize, &crate::variant::Lookup)> = HashSet::new();
        let sites: Vec<tables::Site> = results
            .iter()
            .zip(&stage_of)
            .flat_map(|(allele, n)| allele.variants.iter().map(move |v| (*n, allele.build, v)))
            .filter(|(n, _, v)| listed.insert((*n, &v.lookup)))
            .map(|(n, build, v)| tables::Site {
                build,
                ev: &loaded[n].ev,
                lookup: &v.lookup,
                rsid: &v.rsid,
            })
            .collect();
        let population_rows: Vec<tables::PopulationRow> = sites.iter().flat_map(tables::PopulationRow::of).collect();
        let quality_rows: Vec<tables::QualityRow> = sites.iter().flat_map(tables::QualityRow::of).collect();
        let dataset_rows: Vec<tables::DatasetRow> = loaded.iter().flat_map(|l| tables::DatasetRow::of(&l.build, &l.ev)).collect();
        let alleles_table = tables::alleles_table();
        let mut files = vec![
            ("report.md", report::render(&summary_input)),
            ("final.tsv", alleles_table.render(results.iter().filter(|r| r.verdict.passes()))),
            ("alleles.tsv", alleles_table.render(&results)),
            ("variants.tsv", tables::variants_table().render(&variant_rows)),
            ("populations.tsv", tables::populations_table().render(&population_rows)),
            ("site_quality.tsv", tables::quality_table().render(&quality_rows)),
            ("datasets.tsv", tables::datasets_table().render(&dataset_rows)),
            ("alignments.tsv", tables::alignments_table().render(&placements)),
            ("alignments.txt", alignments),
        ];
        if let Some(calibration) = &model.calibration {
            files.push(("score_levels.tsv", tables::levels_table(&model)?));
            files.push(("score_model.tsv", calibration.format()));
        }
        if haplotypes.coincidence.is_some() {
            files.push((
                "haplotype_calibration.tsv",
                tables::haplotype_table(haplotypes, config.limits.coincidence)?,
            ));
        }

        let mut summary = Vec::new();
        if requests > 0 {
            summary.push(format!("downloaded {:.2} MB in {requests} requests", bytes as f64 / 1e6));
        }
        let count = |status: Status| results.iter().filter(|r| r.placed.status == status).count();
        summary.push(format!(
            "{} alleles: {}",
            results.len(),
            Status::ALL
                .iter()
                .map(|s| format!("{} {}", count(*s), s.as_str()))
                .collect::<Vec<_>>()
                .join(", ")
        ));
        let tally: Vec<String> = report::VERDICTS
            .iter()
            .map(|v| (v, results.iter().filter(|r| r.verdict == *v).count()))
            .filter(|(_, k)| *k > 0)
            .map(|(v, k)| format!("{k} {}", v.as_str()))
            .collect();
        summary.push(format!("verdicts: {}", tally.join(", ")));
        Ok(Output {
            files,
            summary,
            warnings: report::warnings(&summary_input),
        })
    }
}

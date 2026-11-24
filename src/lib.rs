//! igpop: population evidence for germline alleles at given loci.
//!
//! One run ([`run`]) is four steps, each in its own module:
//!
//! 1. `locus`:    place every allele on the locus of its gene, call variants;
//! 2. `evidence`: read population data for those loci (only the bytes needed);
//! 3. `assess`:   class, score, verdict and filters per variant and allele;
//! 4. `tables`, `report`: the result tables and the report, as text.
//!
//! Supporting modules: `align` and `variant` (alignment, normalised
//! variants), `freq` (counts, exact bounds and tests), `score` (the
//! calibration table), `hap` (phased haplotypes), `io`, `tabix`, `fetch`
//! and `net` (remote files), `render` (alignments as text), `log`
//! (progress messages), `par` (parallel map).
//!
//! The library reads its inputs as text and returns its outputs as text,
//! so the same code serves the command-line program (`main.rs`) and the
//! web app (`web/`).

pub mod align;
pub mod assess;
pub mod evidence;
pub mod fetch;
pub mod freq;
pub mod hap;
pub mod io;
pub mod locus;
pub mod log;
pub mod net;
pub mod par;
pub mod render;
pub mod report;
pub mod score;
pub mod tabix;
pub mod tables;
pub mod variant;

use anyhow::{bail, Result};
use assess::{AlleleHaplotype, AlleleResult, HaplotypeCheck, Limits, ScoreModel, Yardstick};
use evidence::{Evidence, Sources};
use freq::KgpRelease;
use locus::{Placed, Status};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

const BUILTIN_COORDS: &str = include_str!("../data/coordinates.tsv");
const BUILTIN_ALIASES: &str = include_str!("../data/gene_aliases.tsv");
const DEFAULT_PEDIGREE: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt";

/// Default remote files of one genome build. `{contig}` is replaced by the
/// contig name.
struct Defaults {
    reference: &'static str,
    kgp: &'static str,
    kgp_release: KgpRelease,
    kgp_name: &'static str,
    gnomad: &'static str,
    gnomad_name: &'static str,
    dbsnp: &'static str,
    phased: &'static str,
}

fn defaults(build: &str) -> Option<Defaults> {
    match build {
        "GRCh38" => Some(Defaults {
            reference: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
            kgp: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{contig}.recalibrated_variants.annotated.vcf.gz",
            kgp_release: KgpRelease::HighCoverage,
            kgp_name: "1kGP_30x",
            gnomad: "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{contig}.vcf.bgz",
            gnomad_name: "gnomAD_v4.1",
            dbsnp: "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
            phased: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{contig}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        }),
        // Phase 3 of 1000 Genomes is one phased genotype file per
        // chromosome: it serves both the frequencies and the haplotypes.
        "GRCh37" => Some(Defaults {
            reference: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
            kgp: PHASE3,
            kgp_release: KgpRelease::Phase3,
            kgp_name: "1kGP_phase3",
            gnomad: "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{contig}.vcf.bgz",
            gnomad_name: "gnomAD_v2.1",
            dbsnp: "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz",
            phased: PHASE3,
        }),
        _ => None,
    }
}

const PHASE3: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz";

/// Where one data source is read from.
#[derive(Clone, Default)]
pub enum Choice {
    /// The default for the build, if there is one.
    #[default]
    Default,
    /// This path or URL.
    Custom(String),
    /// Not used.
    Off,
}

impl Choice {
    /// `main`: for the build the run was asked for. A file given by the
    /// user belongs to that build; any other build uses its default.
    fn resolve(&self, default: Option<&str>, main: bool) -> Option<String> {
        match self {
            Choice::Custom(spec) if main => Some(spec.clone()),
            Choice::Default | Choice::Custom(_) => default.map(str::to_string),
            Choice::Off => None,
        }
    }
}

/// A named text: the content of an input file and what to call it in
/// messages.
#[derive(Clone)]
pub struct Text {
    pub name: String,
    pub content: String,
}

/// Everything a run needs. `Config::new` gives the defaults.
pub struct Config {
    /// Allele FASTA.
    pub alleles: Text,
    /// Coordinates table; `None` = the shipped human IGH table.
    pub coords: Option<Text>,
    /// Genome build, in any spelling `io::normalise_build` accepts.
    pub build: String,
    /// Judge a gene that has coordinates only on another build with the
    /// data of that build.
    pub other_builds: bool,
    /// `allele<TAB>gene` table; `None` = the shipped one, if the shipped
    /// coordinates are used.
    pub aliases: Option<Text>,
    /// Trusted allele names, one per line; `None` = by `candidate_marker`.
    pub trusted: Option<String>,
    /// Alleles whose name contains this are candidates, the others trusted.
    pub candidate_marker: String,
    /// A `score_calibration.tsv` of an earlier run.
    pub calibration: Option<Text>,
    pub limits: Limits,
    /// Share of variants expected to be real, for `probability_real`.
    pub prior: Option<f64>,
    /// Reference FASTA; `None` = the default for the build.
    pub reference: Option<String>,
    pub kgp: Choice,
    pub gnomad: Choice,
    pub dbsnp: Choice,
    pub phased: Choice,
    pub pedigree: Option<String>,
    /// dbSNP datasets left out by name.
    pub dbsnp_exclude: Vec<String>,
    /// Where downloaded regions are kept; `None` = no cache.
    pub cache_dir: Option<PathBuf>,
    /// Reference bases added on each side of a locus.
    pub flank: u64,
    pub min_identity: f64,
    pub paralog_check: bool,
    /// Output prefix, used in the report to name the other files.
    pub prefix: String,
}

impl Config {
    pub fn new(alleles: Text) -> Self {
        Config {
            alleles,
            coords: None,
            build: "38".to_string(),
            other_builds: true,
            aliases: None,
            trusted: None,
            candidate_marker: "_S".to_string(),
            calibration: None,
            limits: Limits {
                common: 0.01,
                supported: 10.0,
                confirm: 0.1,
                coincidence: 0.05,
            },
            prior: None,
            reference: None,
            kgp: Choice::Default,
            gnomad: Choice::Default,
            dbsnp: Choice::Default,
            phased: Choice::Default,
            pedigree: None,
            dbsnp_exclude: vec!["SGDP_PRJ".to_string()],
            cache_dir: None,
            flank: 100,
            min_identity: 0.9,
            paralog_check: true,
            prefix: "igpop".to_string(),
        }
    }

    fn check(&self) -> Result<()> {
        if self.prior.is_some_and(|p| !(p > 0.0 && p < 1.0)) {
            bail!("the prior must lie between 0 and 1");
        }
        if self.limits.confirm >= self.limits.supported {
            bail!("the confirm score must be below the supported score");
        }
        Ok(())
    }
}

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
        let main_build = io::normalise_build(&config.build);

        let (all_coords, coords_label) = match &config.coords {
            Some(text) => (io::parse_coords(&text.content, &text.name)?, format!("`{}`", text.name)),
            None => (
                io::parse_coords(BUILTIN_COORDS, "built-in coordinates")?,
                "built-in human IGH table".to_string(),
            ),
        };
        let on_main: HashSet<&str> = all_coords
            .iter()
            .filter(|c| c.build == main_build)
            .map(|c| c.gene.as_str())
            .collect();
        if on_main.is_empty() {
            bail!("no coordinate rows for build '{main_build}'");
        }
        // A gene is judged on the build that was asked for if it is there,
        // and otherwise on a build that has it and has default data.
        let mut builds = vec![main_build.clone()];
        let mut elsewhere: HashMap<String, String> = HashMap::new();
        for c in all_coords.iter().filter(|c| !on_main.contains(c.gene.as_str())) {
            if config.other_builds && defaults(&c.build).is_some() {
                if !builds.contains(&c.build) {
                    builds.push(c.build.clone());
                }
            } else {
                elsewhere.insert(c.gene.clone(), c.build.clone());
            }
        }

        let aliases = match &config.aliases {
            Some(text) => locus::parse_aliases(&text.content, &text.name)?,
            None if config.coords.is_none() => locus::parse_aliases(BUILTIN_ALIASES, "built-in aliases")?,
            None => HashMap::new(),
        };
        let options = locus::PlaceOptions {
            min_identity: config.min_identity,
            paralog_check: config.paralog_check,
            scoring: align::Scoring::default(),
        };
        let alleles = io::parse_fasta(&config.alleles.content);
        if alleles.is_empty() {
            bail!("{} holds no FASTA records", config.alleles.name);
        }
        let gene_of = |name: &str| -> String { aliases.get(name).cloned().unwrap_or_else(|| locus::gene_of(name).to_string()) };

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
        // no build has stays with the main build, as not placed.
        let home = |name: &str| {
            (0..builds.len())
                .find(|n| coords_of(*n).iter().any(|c| c.gene == gene_of(name)))
                .unwrap_or(0)
        };
        let homes: Vec<usize> = alleles.iter().map(|(name, _)| home(name)).collect();

        let mut stages = Vec::new();
        for (n, build) in builds.iter().enumerate() {
            let main = n == 0;
            let coords = coords_of(n);
            let members: Vec<(String, Vec<u8>)> = alleles
                .iter()
                .zip(&homes)
                .filter(|(_, home)| **home == n)
                .map(|(allele, _)| allele.clone())
                .collect();
            if members.is_empty() && !main {
                continue;
            }
            let defaults = defaults(build);
            let reference = match (&config.reference, &defaults) {
                (Some(spec), _) if main => spec.clone(),
                (_, Some(d)) => d.reference.to_string(),
                _ => bail!("no default reference for build '{build}'; give one"),
            };
            note!("{build} reference: {reference}");
            let windows = locus::load_windows(coords, &mut io::Reference::open(&reference)?, config.flank)?;
            note!("aligning {} alleles to {} loci on {build}", members.len(), windows.len());
            let mut placed = locus::place_all(&members, &windows, &aliases, &options)?;
            placed.sort_by(|a, b| a.allele.cmp(&b.allele));
            let d = defaults.as_ref();
            let sources = Sources {
                kgp: config.kgp.resolve(d.map(|d| d.kgp), main),
                kgp_release: match (&config.kgp, d) {
                    (Choice::Custom(_), _) if main => KgpRelease::HighCoverage,
                    (_, Some(d)) => d.kgp_release,
                    _ => KgpRelease::HighCoverage,
                },
                kgp_name: d.map_or("1kGP", |d| d.kgp_name),
                gnomad: config.gnomad.resolve(d.map(|d| d.gnomad), main),
                gnomad_name: d.map_or("gnomAD", |d| d.gnomad_name),
                dbsnp: config.dbsnp.resolve(d.map(|d| d.dbsnp), main),
                phased: config.phased.resolve(d.map(|d| d.phased), main),
                pedigree: config.pedigree.clone().unwrap_or_else(|| DEFAULT_PEDIGREE.to_string()),
                dbsnp_exclude: config.dbsnp_exclude.clone(),
                common: config.limits.common,
            };
            stages.push(Stage {
                build: build.clone(),
                windows,
                placed,
                sources,
            });
        }
        Ok(Session {
            config,
            coords_label,
            elsewhere,
            stages,
        })
    }

    /// The tabix-indexed files `finish` will read, with their regions.
    pub fn downloads(&self) -> Vec<evidence::Download> {
        self.stages
            .iter()
            .flat_map(|stage| evidence::downloads(&stage.sources, &stage.windows, &stage.placed))
            .collect()
    }

    /// Steps 2 to 4: population evidence, assessment, tables and report.
    pub fn finish(self) -> Result<Output> {
        let Session {
            config,
            coords_label,
            elsewhere,
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
                            "gene absent from {build}: judged on {} with {} and {}",
                            l.build, l.ev.sources.kgp_name, l.ev.sources.gnomad_name
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
                    .map(move |(v, called)| tables::VariantRow {
                        ev,
                        allele,
                        called,
                        v,
                        prior: config.prior,
                    })
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
            prior: config.prior,
            limits: &config.limits,
            results: &results,
            ev,
            model: &model,
            haplotypes,
            requests,
            bytes,
        };
        let alleles_table = tables::alleles_table();
        let mut files = vec![
            ("report.md", report::render(&summary_input)),
            ("final.tsv", alleles_table.render(results.iter().filter(|r| r.verdict.passes()))),
            ("alleles.tsv", alleles_table.render(&results)),
            ("variants.tsv", tables::variants_table().render(&variant_rows)),
            ("alignments.tsv", tables::alignments_table().render(&placements)),
            ("alignments.txt", alignments),
        ];
        if ev.counted() {
            files.push(("score_calibration.tsv", tables::calibration_table(&model)?));
        }
        if !ev.dbsnp.datasets.is_empty() {
            files.push(("dbsnp_datasets.tsv", tables::datasets_table(ev)?));
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

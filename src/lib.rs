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
    kgp: Option<&'static str>,
    gnomad: Option<&'static str>,
    dbsnp: Option<&'static str>,
    phased: Option<&'static str>,
}

fn defaults(build: &str) -> Option<Defaults> {
    match build {
        "GRCh38" => Some(Defaults {
            reference: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
            kgp: Some("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{contig}.recalibrated_variants.annotated.vcf.gz"),
            gnomad: Some("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{contig}.vcf.bgz"),
            dbsnp: Some("https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"),
            phased: Some("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{contig}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        }),
        "GRCh37" => Some(Defaults {
            reference: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
            kgp: None,
            gnomad: None,
            dbsnp: Some("https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"),
            phased: None,
        }),
        _ => None,
    }
}

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
    fn resolve(&self, default: Option<&str>) -> Option<String> {
        match self {
            Choice::Default => default.map(str::to_string),
            Choice::Custom(spec) => Some(spec.clone()),
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

pub fn run(config: &Config) -> Result<Output> {
    config.check()?;
    let build = io::normalise_build(&config.build);
    let defaults = defaults(&build);

    // Step 1: loci and placement.
    let (all_coords, coords_label) = match &config.coords {
        Some(text) => (io::parse_coords(&text.content, &text.name)?, format!("`{}`", text.name)),
        None => (
            io::parse_coords(BUILTIN_COORDS, "built-in coordinates")?,
            "built-in human IGH table".to_string(),
        ),
    };
    // Genes whose coordinates are on another build, to say so in the report.
    let elsewhere: HashMap<String, String> = all_coords
        .iter()
        .filter(|c| c.build != build)
        .map(|c| (c.gene.clone(), c.build.clone()))
        .collect();
    let coords: Vec<io::Coord> = all_coords.into_iter().filter(|c| c.build == build).collect();
    if coords.is_empty() {
        bail!("no coordinate rows for build '{build}'");
    }
    let Some(reference_spec) = config.reference.clone().or(defaults.as_ref().map(|d| d.reference.to_string())) else {
        bail!("no default reference for build '{build}'; give one");
    };
    note!("reference: {reference_spec}");
    let windows = locus::load_windows(coords, &mut io::Reference::open(&reference_spec)?, config.flank)?;
    note!(
        "{} loci on {build} ({} genes on other builds skipped)",
        windows.len(),
        elsewhere.len()
    );

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
    note!("aligning {} alleles", alleles.len());
    let mut placed = locus::place_all(&alleles, &windows, &aliases, &options)?;
    placed.sort_by(|a, b| a.allele.cmp(&b.allele));

    // Step 2: population evidence.
    let default = |pick: fn(&Defaults) -> Option<&'static str>| defaults.as_ref().and_then(pick);
    let sources = Sources {
        kgp: config.kgp.resolve(default(|d| d.kgp)),
        gnomad: config.gnomad.resolve(default(|d| d.gnomad)),
        dbsnp: config.dbsnp.resolve(default(|d| d.dbsnp)),
        phased: config.phased.resolve(default(|d| d.phased)),
        pedigree: config.pedigree.clone().unwrap_or_else(|| DEFAULT_PEDIGREE.to_string()),
        dbsnp_exclude: config.dbsnp_exclude.clone(),
        common: config.limits.common,
    };
    let ev = Evidence::load(sources, &windows, &placed, fetch::Cache(config.cache_dir.as_deref()))?;

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
    let haplotypes = HaplotypeCheck::new(&ev, &windows, &placed);
    let checked: Vec<AlleleHaplotype> = placed
        .iter()
        .enumerate()
        .map(|(i, p)| haplotypes.of(i, p.variants.len(), config.limits.coincidence))
        .collect();
    // The yardstick for the score: cleanly placed alleles that are named
    // as trusted, or whose two or more differences are established
    // together as a haplotype.
    let yardstick: Vec<Yardstick> = placed
        .iter()
        .zip(&checked)
        .map(|(p, haplotype)| match p.status {
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
    let model = ScoreModel::new(&ev, &windows, &calibrating, calibration)?;
    if model.calibration.is_some() {
        note!("score: calibrated from {}", model.description.replace('`', ""));
    }
    let judge = assess::Judge {
        ev: &ev,
        model: &model,
        limits: &config.limits,
    };
    let results: Vec<AlleleResult> = placed
        .iter()
        .zip(checked)
        .enumerate()
        .map(|(i, (p, haplotype))| {
            let other_build = elsewhere.get(&p.gene).map(String::as_str);
            assess::assess_allele(&judge, p, is_trusted(p), yardstick[i], haplotype, other_build)
        })
        .collect();

    // Step 4: tables and report.
    let variant_rows: Vec<tables::VariantRow> = results
        .iter()
        .flat_map(|allele| {
            allele
                .variants
                .iter()
                .zip(&allele.placed.variants)
                .map(|(v, called)| tables::VariantRow {
                    ev: &ev,
                    allele,
                    called,
                    v,
                    prior: config.prior,
                })
        })
        .collect();
    let (requests, bytes) = net::totals();
    let summary_input = report::Input {
        prefix: &config.prefix,
        build: &build,
        alleles: config.alleles.name.clone(),
        coords: coords_label,
        prior: config.prior,
        limits: &config.limits,
        results: &results,
        ev: &ev,
        model: &model,
        haplotypes: &haplotypes,
        requests,
        bytes,
    };
    let alleles_table = tables::alleles_table();
    let mut files = vec![
        ("report.md", report::render(&summary_input)),
        ("final.tsv", alleles_table.render(results.iter().filter(|r| r.verdict.passes()))),
        ("alleles.tsv", alleles_table.render(&results)),
        ("variants.tsv", tables::variants_table().render(&variant_rows)),
        ("alignments.tsv", tables::alignments_table(&build).render(&placed)),
        ("alignments.txt", tables::alignments_text(&build, &windows, &placed)?),
    ];
    if ev.counted() {
        files.push(("score_calibration.tsv", tables::calibration_table(&model)?));
    }
    if !ev.dbsnp.datasets.is_empty() {
        files.push(("dbsnp_datasets.tsv", tables::datasets_table(&ev)?));
    }
    if haplotypes.coincidence.is_some() {
        files.push((
            "haplotype_calibration.tsv",
            tables::haplotype_table(&haplotypes, config.limits.coincidence)?,
        ));
    }

    let mut summary = Vec::new();
    if requests > 0 {
        summary.push(format!("downloaded {:.2} MB in {requests} requests", bytes as f64 / 1e6));
    }
    let count = |status: Status| placed.iter().filter(|p| p.status == status).count();
    summary.push(format!(
        "{} alleles: {}",
        placed.len(),
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

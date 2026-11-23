//! igpop: population evidence for germline alleles at given loci.
//!
//! One run is four steps, each in its own module:
//!
//! 1. `locus`:    place every allele on the locus of its gene, call variants;
//! 2. `evidence`: read population data for those loci (only the bytes needed);
//! 3. `assess`:   class, score, verdict and filters per variant and allele;
//! 4. `tables`, `report`: write the result tables and the report.
//!
//! Supporting modules: `align` and `variant` (alignment, normalised
//! variants), `freq` (counts, exact bounds and tests), `score` (the
//! calibration table), `hap` (phased haplotypes), `io`, `tabix` and `fetch`
//! (remote files), `render` (alignments as text).

mod align;
mod assess;
mod evidence;
mod fetch;
mod freq;
mod hap;
mod io;
mod locus;
mod render;
mod report;
mod score;
mod tabix;
mod tables;
mod variant;

use anyhow::{bail, Context, Result};
use assess::{AlleleResult, Limits, ScoreModel};
use clap::{Parser, Subcommand};
use evidence::{Evidence, Sources};
use locus::{Placed, Status};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

const BUILTIN_COORDS: &str = include_str!("../data/coordinates.tsv");
const BUILTIN_ALIASES: &str = include_str!("../data/gene_aliases.tsv");

/// Default remote sources per build: (1000 Genomes 30x, gnomAD genomes,
/// dbSNP, phased panel). `{contig}` is replaced by the contig name.
fn default_sources(build: &str) -> [Option<&'static str>; 4] {
    match build {
        "GRCh38" => [
            Some("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{contig}.recalibrated_variants.annotated.vcf.gz"),
            Some("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{contig}.vcf.bgz"),
            Some("https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"),
            Some("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{contig}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        ],
        "GRCh37" => [None, None, Some("https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"), None],
        _ => [None; 4],
    }
}
const DEFAULT_PEDIGREE: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt";

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Align alleles to their loci, call variants, judge them against population data.
    Variants(Args),
}

#[derive(Parser)]
struct Args {
    /// Allele FASTA; names like IGHV1-69*01 (gene = text before '*').
    #[arg(long)]
    alleles: PathBuf,
    /// Coordinates TSV: gene, build, contig, start, end, strand.
    /// Default: the human IGH table shipped with igpop.
    #[arg(long)]
    coords: Option<PathBuf>,
    /// Genome build to process: 38 / GRCh38 / hg38 or 37 / GRCh37 / hg19.
    /// Only coordinate rows on this build are used.
    #[arg(long, default_value = "38")]
    build: String,
    /// TSV of `allele<TAB>gene` for alleles whose locus is not the gene in
    /// their name (e.g. IGHV4-4*01 -> IGHV4-4B). `#` lines are skipped.
    /// Default with the shipped coordinates: the shipped alias table.
    #[arg(long)]
    gene_aliases: Option<PathBuf>,
    /// Output prefix: <prefix>.report.md, <prefix>.{final,alleles,variants,alignments}.tsv, ...
    #[arg(long, short, default_value = "igpop")]
    out: String,

    /// Text that marks an allele name as a candidate rather than trusted.
    #[arg(long, default_value = "_S")]
    candidate_marker: String,
    /// File listing trusted allele names, one per line; replaces the marker
    /// rule. Trusted alleles show what real variation looks like and
    /// calibrate the score.
    #[arg(long)]
    trusted: Option<PathBuf>,
    /// Use the scores of an earlier run's <prefix>.score_calibration.tsv.
    /// Default: calibrate on this run, or, with fewer than 30 trusted
    /// single-base variants, use the built-in human IGHV calibration.
    #[arg(long)]
    calibration: Option<PathBuf>,
    /// Allele frequency a variant must demonstrably reach, in at least one
    /// population, to be classed `common`.
    #[arg(long, default_value_t = 0.01)]
    common: f64,
    /// A variant with at least this score is `supported`.
    #[arg(long, default_value_t = 10.0)]
    supported_score: f64,
    /// A variant with at most this score is `needs_confirmation`.
    #[arg(long, default_value_t = 0.1)]
    confirm_score: f64,
    /// Share of variants like yours expected to be real before looking at
    /// population data (between 0 and 1); fills the probability_real column.
    #[arg(long)]
    prior: Option<f64>,

    /// Reference FASTA: a path or URL with .fai (and .gzi if bgzipped).
    /// Default: the 1000 Genomes reference for --build, read remotely.
    #[arg(long)]
    reference: Option<String>,
    /// 1000 Genomes 30x site-level VCF (bgzip + .tbi), path or URL;
    /// `{contig}` is replaced by the contig name. Default for GRCh38.
    #[arg(long)]
    vcf: Option<String>,
    /// gnomAD genomes sites VCF (bgzip + .tbi), path or URL with `{contig}`.
    /// Default for GRCh38: gnomAD v4.1 genomes.
    #[arg(long)]
    gnomad: Option<String>,
    /// Do not use gnomAD.
    #[arg(long)]
    no_gnomad: bool,
    /// dbSNP VCF (bgzip + .tbi): rsIDs, the COMMON flag and the frequencies
    /// relayed from other projects. Default: the latest release for --build.
    #[arg(long)]
    dbsnp: Option<String>,
    /// Do not use dbSNP.
    #[arg(long)]
    no_dbsnp: bool,
    /// dbSNP frequency datasets left out of the evidence, comma-separated.
    /// Default: SGDP_PRJ, whose per-site samples are often one or two
    /// genomes. Re-reports of a source read directly in the same run
    /// (1000 Genomes, gnomAD genomes) are left out automatically.
    #[arg(long, default_value = "SGDP_PRJ")]
    dbsnp_exclude: String,
    /// Phased 1000 Genomes panel (bgzip + .tbi) with genotypes, path or URL
    /// with `{contig}`, for the haplotype check. Default for GRCh38.
    #[arg(long)]
    phased_vcf: Option<String>,
    /// Sample table of the phased panel: family, sample, father, mother,
    /// sex, population, superpopulation. Default: the 1000 Genomes table.
    #[arg(long)]
    pedigree: Option<String>,
    /// Skip the haplotype check.
    #[arg(long)]
    no_haplotypes: bool,
    /// Read no population data at all: alignments and variants only.
    #[arg(long)]
    no_frequencies: bool,
    /// Directory where downloaded regions are kept for later runs.
    #[arg(long, default_value = ".igpop_cache")]
    cache_dir: PathBuf,
    /// Neither read nor write the cache.
    #[arg(long)]
    no_cache: bool,

    /// Reference bases added on each side of a locus.
    #[arg(long, default_value_t = 100)]
    flank: u64,
    /// Alignments below this identity are reported but yield no variants.
    #[arg(long, default_value_t = 0.9)]
    min_identity: f64,
    /// Skip scoring every allele against every other locus.
    #[arg(long)]
    no_paralog_check: bool,
}

impl Args {
    fn check(&self) -> Result<()> {
        if self.prior.is_some_and(|p| !(p > 0.0 && p < 1.0)) {
            bail!("--prior must lie between 0 and 1");
        }
        if self.confirm_score >= self.supported_score {
            bail!("--confirm-score must be below --supported-score");
        }
        Ok(())
    }

    /// Where each source is read from: the option if given, the default for
    /// the build otherwise, nothing if switched off.
    fn sources(&self, build: &str) -> Sources {
        let [kgp, gnomad, dbsnp, phased] = default_sources(build);
        let pick = |given: &Option<String>, default: Option<&str>, off: bool| {
            if off || self.no_frequencies {
                None
            } else {
                given.clone().or(default.map(str::to_string))
            }
        };
        Sources {
            kgp: pick(&self.vcf, kgp, false),
            gnomad: pick(&self.gnomad, gnomad, self.no_gnomad),
            dbsnp: pick(&self.dbsnp, dbsnp, self.no_dbsnp),
            phased: pick(&self.phased_vcf, phased, self.no_haplotypes),
            pedigree: self.pedigree.clone().unwrap_or_else(|| DEFAULT_PEDIGREE.to_string()),
            dbsnp_exclude: self
                .dbsnp_exclude
                .split(',')
                .map(|x| x.trim().to_string())
                .filter(|x| !x.is_empty())
                .collect(),
            common: self.common,
        }
    }
}

fn run(args: Args) -> Result<()> {
    args.check()?;
    let build = io::normalise_build(&args.build);

    // Step 1: loci and placement.
    let (all_coords, coords_label) = match &args.coords {
        Some(path) => (io::read_coords(path)?, format!("`{}`", path.display())),
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
    let reference_spec = match (&args.reference, build.as_str()) {
        (Some(r), _) => r.clone(),
        (None, "GRCh38") => io::DEFAULT_GRCH38.to_string(),
        (None, "GRCh37") => io::DEFAULT_GRCH37.to_string(),
        (None, other) => bail!("no default reference for build '{other}'; pass --reference"),
    };
    eprintln!("reference: {reference_spec}");
    let windows = locus::load_windows(coords, &mut io::Reference::open(&reference_spec)?, args.flank)?;
    eprintln!(
        "{} loci on {build} ({} rows on other builds skipped)",
        windows.len(),
        elsewhere.len()
    );

    let aliases = match &args.gene_aliases {
        Some(path) => {
            let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
            locus::parse_aliases(&text, &path.display().to_string())?
        }
        None if args.coords.is_none() => locus::parse_aliases(BUILTIN_ALIASES, "built-in aliases")?,
        None => HashMap::new(),
    };
    let options = locus::PlaceOptions {
        min_identity: args.min_identity,
        paralog_check: !args.no_paralog_check,
        scoring: align::Scoring::default(),
    };
    let mut placed: Vec<Placed> = io::read_fasta(&args.alleles)?
        .par_iter()
        .map(|(name, seq)| locus::place(name, seq, &windows, &aliases, &options))
        .collect::<Result<_>>()?;
    placed.sort_by(|a, b| a.allele.cmp(&b.allele));

    // Step 2: population evidence.
    let cache = fetch::Cache(if args.no_cache { None } else { Some(args.cache_dir.as_path()) });
    let ev = Evidence::load(args.sources(&build), &windows, &placed, cache)?;

    // Step 3: calibration and assessment.
    let trusted_names: Option<HashSet<String>> = match &args.trusted {
        Some(path) => Some(
            std::fs::read_to_string(path)
                .with_context(|| format!("reading {}", path.display()))?
                .lines()
                .map(|l| l.trim().to_string())
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .collect(),
        ),
        None => None,
    };
    let is_trusted = |p: &Placed| match &trusted_names {
        Some(names) => names.contains(&p.allele),
        None => !p.allele.contains(&args.candidate_marker),
    };
    // Only cleanly placed trusted alleles calibrate the score.
    let calibrating: Vec<&Placed> = placed.iter().filter(|p| p.status == Status::Ok && is_trusted(p)).collect();
    let model = ScoreModel::new(&ev, &windows, &calibrating, args.calibration.as_deref())?;
    if model.calibration.is_some() {
        eprintln!("score: calibrated from {}", model.description.replace('`', ""));
    }
    let limits = Limits {
        common: args.common,
        supported: args.supported_score,
        confirm: args.confirm_score,
    };
    let results: Vec<AlleleResult> = placed
        .iter()
        .map(|p| assess::assess_allele(&ev, p, is_trusted(p), elsewhere.get(&p.gene).map(String::as_str), &model, &limits))
        .collect();

    // Step 4: tables and report.
    let file = |name: &str| format!("{}.{name}", args.out);
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
                    prior: args.prior,
                })
        })
        .collect();
    tables::variants_table().write(&file("variants.tsv"), &variant_rows)?;
    let alleles = tables::alleles_table();
    alleles.write(&file("alleles.tsv"), &results)?;
    alleles.write(&file("final.tsv"), results.iter().filter(|r| r.verdict.passes()))?;
    tables::alignments_table(&build).write(&file("alignments.tsv"), &placed)?;
    tables::write_alignments(&file("alignments.txt"), &build, &windows, &placed)?;
    if ev.counted() {
        tables::write_calibration(&file("score_calibration.tsv"), &model)?;
    }
    if !ev.dbsnp.datasets.is_empty() {
        tables::write_datasets(&file("dbsnp_datasets.tsv"), &ev)?;
    }

    let (requests, bytes) = io::download_totals();
    if requests > 0 {
        eprintln!("downloaded {:.2} MB in {requests} requests", bytes as f64 / 1e6);
    }
    let summary = report::Input {
        prefix: &args.out,
        build: &build,
        alleles: args.alleles.display().to_string(),
        coords: coords_label,
        prior: args.prior,
        limits: &limits,
        results: &results,
        ev: &ev,
        model: &model,
        requests,
        bytes,
    };
    std::fs::write(file("report.md"), report::render(&summary))?;

    let count = |status: Status| placed.iter().filter(|p| p.status == status).count();
    eprintln!(
        "{} alleles: {}",
        placed.len(),
        Status::ALL
            .iter()
            .map(|s| format!("{} {}", count(*s), s.as_str()))
            .collect::<Vec<_>>()
            .join(", ")
    );
    let tally: Vec<String> = report::VERDICTS
        .iter()
        .map(|v| (v, results.iter().filter(|r| r.verdict == *v).count()))
        .filter(|(_, k)| *k > 0)
        .map(|(v, k)| format!("{k} {}", v.as_str()))
        .collect();
    eprintln!("verdicts: {}", tally.join(", "));
    for w in report::warnings(&summary) {
        eprintln!("warning: {}", w.replace('`', ""));
    }
    eprintln!("read {} first", file("report.md"));
    Ok(())
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Variants(args) => run(args),
    }
}

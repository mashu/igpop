//! The command-line program: options and files in, files out. Everything
//! else is in the library (`lib.rs`).

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use igpop::{Choice, Config, Text};
use std::path::{Path, PathBuf};

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
    /// The exact sequence of an allele, found on phased chromosomes, makes
    /// it `supported` if fewer than this share of arbitrary single-base
    /// changes of real chromosomes are found as often.
    #[arg(long, default_value_t = 0.05)]
    coincidence: f64,
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

fn read(path: &Path) -> Result<Text> {
    Ok(Text {
        name: path.display().to_string(),
        content: std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?,
    })
}

impl Args {
    fn config(&self) -> Result<Config> {
        let optional = |path: &Option<PathBuf>| path.as_deref().map(read).transpose();
        let choice = |given: &Option<String>, off: bool| match given {
            _ if off || self.no_frequencies => Choice::Off,
            Some(spec) => Choice::Custom(spec.clone()),
            None => Choice::Default,
        };
        let mut config = Config::new(read(&self.alleles)?);
        config.coords = optional(&self.coords)?;
        config.build = self.build.clone();
        config.aliases = optional(&self.gene_aliases)?;
        config.trusted = optional(&self.trusted)?.map(|t| t.content);
        config.candidate_marker = self.candidate_marker.clone();
        config.calibration = optional(&self.calibration)?;
        config.limits.common = self.common;
        config.limits.supported = self.supported_score;
        config.limits.confirm = self.confirm_score;
        config.limits.coincidence = self.coincidence;
        config.prior = self.prior;
        config.reference = self.reference.clone();
        config.kgp = choice(&self.vcf, false);
        config.gnomad = choice(&self.gnomad, self.no_gnomad);
        config.dbsnp = choice(&self.dbsnp, self.no_dbsnp);
        config.phased = choice(&self.phased_vcf, self.no_haplotypes);
        config.pedigree = self.pedigree.clone();
        config.dbsnp_exclude = self
            .dbsnp_exclude
            .split(',')
            .map(|x| x.trim().to_string())
            .filter(|x| !x.is_empty())
            .collect();
        config.cache_dir = (!self.no_cache).then(|| self.cache_dir.clone());
        config.flank = self.flank;
        config.min_identity = self.min_identity;
        config.paralog_check = !self.no_paralog_check;
        config.prefix = self.out.clone();
        Ok(config)
    }
}

fn main() -> Result<()> {
    let Command::Variants(args) = Cli::parse().command;
    let output = igpop::run(&args.config()?)?;
    for (name, content) in &output.files {
        let path = format!("{}.{name}", args.out);
        std::fs::write(&path, content).with_context(|| format!("writing {path}"))?;
    }
    for line in &output.summary {
        eprintln!("{line}");
    }
    for warning in &output.warnings {
        eprintln!("warning: {}", warning.replace('`', ""));
    }
    eprintln!("read {}.report.md first", args.out);
    Ok(())
}

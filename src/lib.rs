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

/// The shipped tables: human IGH loci, and the alleles whose gene is not
/// the part of their name before `*`.
pub const BUILTIN_COORDS: &str = include_str!("../data/coordinates.tsv");
pub const BUILTIN_ALIASES: &str = include_str!("../data/gene_aliases.tsv");
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
    /// Place an allele whose gene is not in the coordinates table on the
    /// locus its sequence fits best, instead of leaving it unplaced. The
    /// alias table asks the same for single alleles or genes with the
    /// gene `?` (`BY_SEQUENCE`).
    pub place_unknown: bool,
    /// Trusted allele names, one per line; `None` = by `candidate_marker`.
    pub trusted: Option<String>,
    /// Alleles whose name contains this are candidates, the others trusted.
    pub candidate_marker: String,
    /// A `score_model.tsv` of an earlier run.
    pub calibration: Option<Text>,
    pub limits: Limits,
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
            place_unknown: false,
            trusted: None,
            candidate_marker: "_S".to_string(),
            calibration: None,
            limits: Limits {
                common: 0.01,
                supported: 2.0,
                confirm: 0.1,
                coincidence: 0.05,
            },
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

/// Which genes a run knows, and on which build each is judged.
struct Plan {
    /// The build that was asked for, then any other build that holds
    /// genes the first one lacks.
    builds: Vec<String>,
    all_coords: Vec<io::Coord>,
    coords_label: String,
    /// Genes whose only coordinates are on a build this run does not use.
    elsewhere: HashMap<String, String>,
}

impl Plan {
    /// The plan, and the alias table of the run.
    fn new(config: &Config) -> Result<(Self, HashMap<String, String>)> {
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
        let plan = Plan {
            builds,
            all_coords,
            coords_label,
            elsewhere,
        };
        Ok((plan, aliases))
    }

    /// Where a gene is judged: the first build of the run that has it.
    fn place_of(&self, gene: &str) -> Option<&io::Coord> {
        self.builds
            .iter()
            .find_map(|build| self.all_coords.iter().find(|c| &c.build == build && c.gene == gene))
    }
}

/// The alias-table gene that means "find the gene by sequence".
pub const BY_SEQUENCE: &str = "?";

/// One allele of the file: its name, the gene its name gives (the part
/// before `*`), and the gene it is assigned to after the alias table.
pub struct AllelePreview {
    pub name: String,
    pub named: String,
    pub gene: String,
}

/// One gene the alleles of the file are assigned to.
pub struct GenePreview {
    pub gene: String,
    /// Where it is judged: build, contig, first and last position; `None`
    /// if the gene is not in the coordinates table of this run.
    pub place: Option<(String, String, u64, u64)>,
    /// Population data are read for it: it lies on a whole chromosome.
    pub evidence: bool,
}

/// A data source a run will read for the genes of one build.
pub struct PlannedSource {
    pub build: String,
    /// The name the output uses, e.g. `1kGP_30x`.
    pub name: &'static str,
    pub title: &'static str,
    /// `frequencies` or `haplotypes`.
    pub gives: &'static str,
    /// People in the dataset, where that is one number.
    pub people: Option<u32>,
    pub about: &'static str,
}

/// What a run will do, known before anything is downloaded.
pub struct Preview {
    pub alleles: Vec<AllelePreview>,
    pub genes: Vec<GenePreview>,
    /// The sources of every build that has a gene of the file.
    pub sources: Vec<PlannedSource>,
}

/// The sources a run reads for one build, given what the user switched off
/// or replaced.
fn planned(config: &Config, build: &str, main: bool) -> Vec<PlannedSource> {
    let Some(d) = defaults(build) else { return Vec::new() };
    let phase3 = d.kgp_release == KgpRelease::Phase3;
    let own = |choice: &Choice| main && matches!(choice, Choice::Custom(_));
    let mut out = Vec::new();
    let mut add = |choice: &Choice, default: &str, source: PlannedSource| {
        if choice.resolve(Some(default), main).is_some() {
            out.push(match own(choice) {
                true => PlannedSource {
                    title: "your own file",
                    people: None,
                    about: "given on the command line",
                    ..source
                },
                false => source,
            });
        }
    };
    let source = |name, title, gives, people, about| PlannedSource {
        build: build.to_string(),
        name,
        title,
        gives,
        people,
        about,
    };
    add(
        &config.kgp,
        d.kgp,
        match phase3 {
            true => source(
                d.kgp_name,
                "1000 Genomes, phase 3",
                "frequencies",
                Some(2504),
                "five continental groups; low sequencing depth, no site quality measures",
            ),
            false => source(
                d.kgp_name,
                "1000 Genomes, 30x",
                "frequencies",
                Some(2504),
                "unrelated people in five continental groups, with site quality measures",
            ),
        },
    );
    add(
        &config.gnomad,
        d.gnomad,
        match phase3 {
            true => source(
                d.gnomad_name,
                "gnomAD v2.1 genomes",
                "frequencies",
                Some(15708),
                "six genetic ancestry groups",
            ),
            false => source(
                d.gnomad_name,
                "gnomAD v4.1 genomes",
                "frequencies",
                Some(76215),
                "nine genetic ancestry groups",
            ),
        },
    );
    add(
        &config.dbsnp,
        d.dbsnp,
        source(
            "dbSNP",
            "Studies collected by dbSNP",
            "frequencies",
            None,
            "rsIDs, and frequencies from further studies of a few hundred to several hundred thousand people",
        ),
    );
    add(
        &config.phased,
        d.phased,
        match phase3 {
            true => source(
                "phased",
                "1000 Genomes phase 3, phased",
                "haplotypes",
                Some(2504),
                "whole sequences on single chromosomes, for the haplotype check",
            ),
            false => source(
                "phased",
                "1000 Genomes phased panel",
                "haplotypes",
                Some(2598),
                "whole sequences on single chromosomes, for the haplotype check",
            ),
        },
    );
    out
}

/// The genes of the allele file, whether the run knows them and which
/// data it will read for them: a check to make before a run, without any
/// download.
pub fn preview(config: &Config) -> Result<Preview> {
    let (plan, aliases) = Plan::new(config)?;
    let alleles: Vec<AllelePreview> = io::parse_fasta(&config.alleles.content)
        .into_iter()
        .map(|(name, _)| AllelePreview {
            named: locus::gene_of(&name).to_string(),
            gene: locus::gene_for(&name, &aliases).to_string(),
            name,
        })
        .collect();
    if alleles.is_empty() {
        bail!("{} holds no FASTA records", config.alleles.name);
    }
    let mut genes: Vec<GenePreview> = Vec::new();
    for allele in &alleles {
        if genes.iter().any(|g| g.gene == allele.gene) {
            continue;
        }
        let place = plan.place_of(&allele.gene);
        genes.push(GenePreview {
            gene: allele.gene.clone(),
            place: place.map(|c| (c.build.clone(), c.contig.clone(), c.start, c.end)),
            evidence: place.is_some_and(|c| evidence::whole_chromosome(&c.contig)),
        });
    }
    let sources = plan
        .builds
        .iter()
        .enumerate()
        .filter(|(_, build)| genes.iter().any(|g| g.evidence && g.place.as_ref().is_some_and(|p| &p.0 == *build)))
        .flat_map(|(n, build)| planned(config, build, n == 0))
        .collect();
    Ok(Preview { alleles, genes, sources })
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
            let defaults = defaults(build);
            let reference = match (&config.reference, &defaults) {
                (Some(spec), _) if main => spec.clone(),
                (_, Some(d)) => d.reference.to_string(),
                _ => bail!("no default reference for build '{build}'; give one"),
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
            by_sequence,
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
            files.push(("score_levels.tsv", tables::calibration_table(&model)?));
            if let Some(calibration) = &model.calibration {
                files.push(("score_model.tsv", calibration.format()));
            }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preview_tells_known_genes_from_unknown_ones() {
        let fasta = ">IGHV1-18*01\nACGT\n>IGHV4-4*02_S1\nACGT\n>IGHV3-9*01\nACGT\n>IGHVX-1*01\nACGT\n>IGHVX-1*02\nACGT\n";
        let config = Config::new(Text {
            name: "test.fa".to_string(),
            content: fasta.to_string(),
        });
        let Preview { alleles, genes, sources } = preview(&config).unwrap();
        // Both builds are needed, each with its four default sources.
        assert_eq!(sources.iter().filter(|s| s.build == "GRCh38").count(), 4);
        assert_eq!(sources.iter().filter(|s| s.build == "GRCh37").count(), 4);
        let of = |gene: &str| genes.iter().find(|g| g.gene == gene).unwrap();
        assert_eq!(of("IGHV1-18").place.as_ref().unwrap().0, "GRCh38");
        // Named IGHV4-4 after a database allele of the gene the table calls IGHV4-4B.
        let renamed = alleles.iter().find(|a| a.name == "IGHV4-4*02_S1").unwrap();
        assert_eq!((renamed.named.as_str(), renamed.gene.as_str()), ("IGHV4-4", "IGHV4-4B"));
        assert!(of("IGHV4-4B").place.is_some());
        // Absent from GRCh38: judged on the build that has it.
        assert_eq!(of("IGHV3-9").place.as_ref().unwrap().0, "GRCh37");
        assert!(of("IGHVX-1").place.is_none());
        assert_eq!(alleles.iter().filter(|a| a.gene == "IGHVX-1").count(), 2);
    }
}

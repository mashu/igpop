//! Before anything is downloaded: which genes a run knows, on which build
//! each is judged, and which data will be read for it.

use crate::catalogue::{self, Slot};
use crate::config::{Choice, Config};
use crate::evidence::{self, Sources};
use crate::population::Dataset;
use crate::{io, locus};
use anyhow::{bail, Result};
use std::collections::{HashMap, HashSet};

/// Which genes a run knows, and on which build each is judged.
pub(crate) struct Plan {
    /// The build that was asked for, then any other build that holds
    /// genes the first one lacks.
    pub builds: Vec<String>,
    pub all_coords: Vec<io::Coord>,
    pub coords_label: String,
    /// Genes whose only coordinates are on a build this run does not use.
    pub elsewhere: HashMap<String, String>,
}

impl Plan {
    /// The plan, and the alias table of the run.
    pub fn new(config: &Config) -> Result<(Self, HashMap<String, String>)> {
        let main_build = io::normalise_build(&config.build);
        let (all_coords, coords_label) = match &config.coords {
            Some(text) => (io::parse_coords(&text.content, &text.name)?, format!("`{}`", text.name)),
            None => (
                io::parse_coords(catalogue::COORDINATES, "built-in coordinates")?,
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
            if config.other_builds && catalogue::build(&c.build).is_some() {
                if !builds.contains(&c.build) {
                    builds.push(c.build.clone());
                }
            } else {
                elsewhere.insert(c.gene.clone(), c.build.clone());
            }
        }
        let aliases = match &config.aliases {
            Some(text) => locus::parse_aliases(&text.content, &text.name)?,
            None if config.coords.is_none() => locus::parse_aliases(catalogue::ALIASES, "built-in aliases")?,
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

/// What a run reads for one build, given what the user switched off or
/// replaced: the reference and the sources.
pub(crate) struct BuildData {
    pub reference: Option<String>,
    pub sources: Sources,
    /// The phased panel in words, if it is read.
    phased: Option<(&'static str, Option<u32>)>,
}

impl BuildData {
    /// `main`: the build the run was asked for. A file given by the user
    /// belongs to that build; any other build uses its defaults.
    pub fn new(config: &Config, build: &str, main: bool) -> Self {
        let defaults = catalogue::build(build);
        let d = defaults.as_ref();
        // A file of the user's takes the place of the default dataset and
        // is read the same way.
        let own = |choice: &Choice, dataset: Dataset| match choice {
            Choice::Custom(spec) if main => Some(Dataset {
                title: "your own file",
                about: "given on the command line",
                people: None,
                spec: spec.clone(),
                ..dataset
            }),
            Choice::Off => None,
            _ => Some(dataset),
        };
        let datasets = d
            .map(|d| d.datasets.clone())
            .unwrap_or_default()
            .into_iter()
            .filter_map(|(slot, dataset)| match slot {
                Slot::Kgp => own(&config.kgp, dataset),
                Slot::Gnomad => own(&config.gnomad, dataset),
                Slot::More => config.more_datasets.then_some(dataset),
            })
            .collect();
        let phased = config.phased.resolve(d.map(|d| d.phased.spec), main);
        BuildData {
            reference: match &config.reference {
                Some(spec) if main => Some(spec.clone()),
                _ => d.map(|d| d.reference.to_string()),
            },
            phased: phased.as_ref().map(|_| match (&config.phased, d) {
                (Choice::Custom(_), _) if main => ("your own file", None),
                (_, Some(d)) => (d.phased.title, Some(d.phased.people)),
                _ => ("phased panel", None),
            }),
            sources: Sources {
                datasets,
                dbsnp: config.dbsnp.resolve(d.map(|d| d.dbsnp), main),
                dbsnp_exclude: config.dbsnp_exclude.clone(),
                phased,
                pedigree: config.pedigree.clone().unwrap_or_else(|| catalogue::PEDIGREE.to_string()),
                common: config.limits.common,
            },
        }
    }

    /// The sources as a list for the preview.
    fn planned(&self, build: &str) -> Vec<PlannedSource> {
        let source = |name, title, gives, people, about| PlannedSource {
            build: build.to_string(),
            name,
            title,
            gives,
            people,
            about,
        };
        let mut out: Vec<PlannedSource> = self
            .sources
            .datasets
            .iter()
            .map(|d| source(d.name, d.title, "frequencies", d.people, d.about))
            .collect();
        if self.sources.dbsnp.is_some() {
            // dbSNP comes after the two datasets every build has.
            let at = out.len().min(2);
            out.insert(
                at,
                source("dbSNP", catalogue::DBSNP_TITLE, "frequencies", None, catalogue::DBSNP_ABOUT),
            );
        }
        if let Some((title, people)) = self.phased {
            out.push(source("phased", title, "haplotypes", people, catalogue::PHASED_ABOUT));
        }
        out
    }
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
        .flat_map(|(n, build)| BuildData::new(config, build, n == 0).planned(build))
        .collect();
    Ok(Preview { alleles, genes, sources })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Text;

    #[test]
    fn preview_tells_known_genes_from_unknown_ones() {
        let fasta = ">IGHV1-18*01\nACGT\n>IGHV4-4*02_S1\nACGT\n>IGHV3-9*01\nACGT\n>IGHVX-1*01\nACGT\n>IGHVX-1*02\nACGT\n";
        let config = Config::new(Text {
            name: "test.fa".to_string(),
            content: fasta.to_string(),
        });
        let Preview { alleles, genes, sources } = preview(&config).unwrap();
        // Both builds are needed, each with its four default sources.
        assert_eq!(sources.iter().filter(|s| s.build == "GRCh38").count(), 8);
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

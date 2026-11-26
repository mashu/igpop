//! What igpop reads by default for each genome build: the reference, the
//! population datasets, dbSNP and the phased panel; and the shipped gene
//! tables. Adding a dataset to igpop is adding an entry here.

use crate::population::{Dataset, Layout, ALFA_GROUPS, GENOME_ASIA_GROUPS, GNOMAD_GROUPS, KGP_GROUPS, TOMMO_GROUPS, WBBC_GROUPS};

/// The shipped tables: human IGH loci, and the alleles whose gene is not
/// the part of their name before `*`.
pub const COORDINATES: &str = include_str!("../data/coordinates.tsv");
pub const ALIASES: &str = include_str!("../data/gene_aliases.tsv");

/// Who is whose parent in 1000 Genomes: the haplotype check counts
/// founders only.
pub const PEDIGREE: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt";

/// Phase 3 of 1000 Genomes is one phased genotype file per chromosome: it
/// serves both the frequencies and the haplotypes of GRCh37.
const PHASE3: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{contig}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz";

/// Which option of a run switches a dataset off or replaces its file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Slot {
    /// The 1000 Genomes frequencies.
    Kgp,
    Gnomad,
    /// The further datasets, switched together.
    More,
}

/// The phased panel of a build, for the haplotype check.
pub struct Phased {
    pub spec: &'static str,
    pub title: &'static str,
    pub people: u32,
}

/// The default data of one genome build.
pub struct Build {
    pub reference: &'static str,
    pub datasets: Vec<(Slot, Dataset)>,
    pub dbsnp: &'static str,
    pub phased: Phased,
}

/// What dbSNP gives, for a list of sources.
pub const DBSNP_TITLE: &str = "Studies collected by dbSNP";
pub const DBSNP_ABOUT: &str = "rsIDs, and frequencies from further studies of a few hundred to several hundred thousand people";
pub const PHASED_ABOUT: &str = "whole sequences on single chromosomes, for the haplotype check";

fn kgp(name: &'static str, title: &'static str, about: &'static str, spec: &str, layout: Layout) -> (Slot, Dataset) {
    (
        Slot::Kgp,
        Dataset {
            name,
            title,
            about,
            people: Some(2504),
            spec: spec.to_string(),
            groups: &KGP_GROUPS,
            layout,
            complete: true,
            replaces: &["1000Genomes", "1000Genomes_30X"],
            required: true,
        },
    )
}

fn gnomad(name: &'static str, title: &'static str, about: &'static str, people: u32, spec: &str) -> (Slot, Dataset) {
    (
        Slot::Gnomad,
        Dataset {
            name,
            title,
            about,
            people: Some(people),
            spec: spec.to_string(),
            groups: &GNOMAD_GROUPS,
            layout: Layout::Gnomad,
            complete: true,
            replaces: &["GnomAD_genomes"],
            required: true,
        },
    )
}

/// The further datasets of GRCh38: populations the first two cover thinly
/// or as one block.
fn more() -> Vec<(Slot, Dataset)> {
    let ucsc = |file: &str| format!("https://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/{file}");
    [
        Dataset {
            name: "ALFA",
            title: "ALFA (NCBI), by population",
            about: "people of studies held at NCBI, in nine populations; how many were examined differs from position to position",
            people: None,
            spec: "https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz".to_string(),
            groups: &ALFA_GROUPS,
            // The last three sample columns are sums of the others.
            layout: Layout::Samples {
                columns: &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            },
            complete: false,
            replaces: &["dbGaP_PopFreq"],
            required: false,
        },
        Dataset {
            name: "GenomeAsia",
            title: "GenomeAsia pilot",
            about: "219 population groups of Asia and a worldwide comparison set, in seven groups including Southeast Asian",
            people: Some(1739),
            spec: ucsc("ga100k/ga100k.subst.vcf.gz"),
            groups: &GENOME_ASIA_GROUPS,
            layout: Layout::Counts {
                keys: &[
                    ("AC_NEA", "AN_NEA"),
                    ("AC_SEA", "AN_SEA"),
                    ("AC_SAS", "AN_SAS"),
                    ("AC_OCE", "AN_OCE"),
                    ("AC_AFR", "AN_AFR"),
                    ("AC_AMR", "AN_AMR"),
                    ("AC_WER", "AN_WER"),
                ],
            },
            complete: true,
            replaces: &[],
            required: false,
        },
        Dataset {
            name: "ToMMo_61k",
            title: "ToMMo 61KJPN",
            about: "Japanese people, sequenced as whole genomes",
            people: Some(61_000),
            spec: ucsc("tommo61kjpn/tommo-61kjpn-20250616-GRCh38-snvindel-af-autosome.vcf.gz"),
            groups: &TOMMO_GROUPS,
            layout: Layout::Counts { keys: &[("AC", "AN")] },
            complete: true,
            replaces: &["TOMMO"],
            required: false,
        },
        Dataset {
            name: "WBBC",
            title: "Westlake BioBank for Chinese",
            about: "Chinese people from four regions of China",
            people: Some(4480),
            spec: ucsc("wbbc/wbbc.vcf.gz"),
            groups: &WBBC_GROUPS,
            layout: Layout::Frequencies {
                keys: &[
                    ("North_AF", "North_AN"),
                    ("Central_AF", "Central_AN"),
                    ("South_AF", "South_AN"),
                    ("Lingnan_AF", "Lingnan_AN"),
                ],
            },
            complete: true,
            replaces: &[],
            required: false,
        },
    ]
    .into_iter()
    .map(|dataset| (Slot::More, dataset))
    .collect()
}

/// The default data of a build, if igpop has any for it.
pub fn build(name: &str) -> Option<Build> {
    match name {
        "GRCh38" => {
            let mut datasets = vec![
                kgp(
                    "1kGP_30x",
                    "1000 Genomes, 30x",
                    "unrelated people in five continental groups, with site quality measures",
                    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{contig}.recalibrated_variants.annotated.vcf.gz",
                    Layout::Kgp30x,
                ),
                gnomad(
                    "gnomAD_v4.1",
                    "gnomAD v4.1 genomes",
                    "nine genetic ancestry groups",
                    76_215,
                    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{contig}.vcf.bgz",
                ),
            ];
            datasets.extend(more());
            Some(Build {
                reference: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
                datasets,
                dbsnp: "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
                phased: Phased {
                    spec: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{contig}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
                    title: "1000 Genomes phased panel",
                    people: 2598,
                },
            })
        }
        "GRCh37" => Some(Build {
            reference: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
            datasets: vec![
                kgp(
                    "1kGP_phase3",
                    "1000 Genomes, phase 3",
                    "five continental groups; low sequencing depth, no site quality measures",
                    PHASE3,
                    Layout::KgpPhase3,
                ),
                gnomad(
                    "gnomAD_v2.1",
                    "gnomAD v2.1 genomes",
                    "six genetic ancestry groups",
                    15_708,
                    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{contig}.vcf.bgz",
                ),
            ],
            dbsnp: "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz",
            phased: Phased {
                spec: PHASE3,
                title: "1000 Genomes phase 3, phased",
                people: 2504,
            },
        }),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A layout names one key or column per group: the two lists of an
    /// entry must have the same length, and output names must be unique.
    #[test]
    fn every_entry_is_consistent() {
        for name in ["GRCh38", "GRCh37"] {
            let build = build(name).unwrap();
            let mut names: Vec<&str> = build.datasets.iter().map(|(_, d)| d.name).collect();
            names.sort_unstable();
            names.dedup();
            assert_eq!(names.len(), build.datasets.len());
            for (_, dataset) in &build.datasets {
                let fields = match dataset.layout {
                    Layout::Samples { columns } => columns.len(),
                    Layout::Counts { keys } | Layout::Frequencies { keys } => keys.len(),
                    Layout::Kgp30x | Layout::KgpPhase3 | Layout::Gnomad => dataset.groups.len(),
                };
                assert_eq!(fields, dataset.groups.len(), "{}", dataset.name);
            }
        }
    }
}

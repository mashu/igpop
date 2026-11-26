//! Population datasets: who is in them, how their files are read, and
//! what one record says about one variant.
//!
//! A dataset is one tabix-indexed VCF with allele counts per group of
//! people. Datasets differ only in where a record keeps those counts
//! (`Layout`) and in which site-quality measures it carries (`Quality`).
//! Everything downstream sees the same thing for all of them: a list of
//! counts, one per group. Adding a dataset means adding an entry to the
//! catalogue (`crate::catalogue`), and a `Layout` only if its file is laid
//! out in a way none of the existing ones covers.
//!
//! dbSNP is different in kind, since it relays frequencies of studies
//! without counts, and has its own module.

pub mod dbsnp;
mod layout;

pub use layout::Layout;

use crate::stats::{self, PopCount};

/// The continental group a population belongs to. It gives populations of
/// different datasets a common colour in pictures; no statistic uses it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ancestry {
    African,
    American,
    EastAsian,
    European,
    MiddleEastern,
    SouthAsian,
    /// Mixed, or none of the above.
    Other,
}

impl Ancestry {
    pub fn as_str(self) -> &'static str {
        match self {
            Ancestry::African => "african",
            Ancestry::American => "american",
            Ancestry::EastAsian => "east_asian",
            Ancestry::European => "european",
            Ancestry::MiddleEastern => "middle_eastern",
            Ancestry::SouthAsian => "south_asian",
            Ancestry::Other => "other",
        }
    }
}

/// One group of people of a dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Group {
    /// As written in the output, e.g. `AFR`.
    pub code: &'static str,
    /// In words, e.g. `African`.
    pub title: &'static str,
    pub ancestry: Ancestry,
}

/// A dataset with allele counts per group.
#[derive(Clone, Debug, PartialEq)]
pub struct Dataset {
    /// As written in the output, e.g. `1kGP_30x`.
    pub name: &'static str,
    pub title: &'static str,
    /// Who is in it, in a few words.
    pub about: &'static str,
    /// People in the dataset, where that is one number.
    pub people: Option<u32>,
    /// Path or URL of the VCF; `{contig}` stands for the contig name.
    pub spec: String,
    pub groups: &'static [Group],
    pub layout: Layout,
    /// Every person was examined at every position, so a variant without
    /// a record was looked for and not seen. False where the number of
    /// people differs from position to position without a record saying so.
    pub complete: bool,
    /// Studies relayed by dbSNP that report the same people.
    pub replaces: &'static [&'static str],
    /// A run cannot do without it: failing to read it stops the run. Any
    /// other dataset that cannot be read is reported and left out.
    pub required: bool,
}

/// What one record of a dataset says about one alternate allele.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    /// One count per group, in the order of `Dataset::groups`.
    pub counts: Vec<PopCount>,
    pub quality: Quality,
}

impl Record {
    /// All people of the dataset: the published total where there is one,
    /// otherwise the sum of the groups.
    pub fn all(&self) -> PopCount {
        self.quality.all.unwrap_or_else(|| stats::total(&self.counts))
    }
}

/// Site-quality measures of a record, as far as its dataset publishes
/// them. They say how far the counts can be trusted, never how frequent
/// the allele is.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Quality {
    /// The FILTER column as written.
    pub filter: String,
    /// Region flags set on the record: `lcr`, `segdup`.
    pub flags: Vec<&'static str>,
    /// Counts over all people, where the groups do not add up to them.
    pub all: Option<PopCount>,
    /// People with one copy and with two.
    pub heterozygous: Option<u32>,
    pub homozygous: Option<u32>,
    /// Hardy-Weinberg exact-test p-value per group.
    pub hwe: Vec<f64>,
    /// Excess-heterozygosity p-value over all people (small = excess).
    pub excess_het: Option<f64>,
    /// Mendelian error rate in the trios.
    pub mendel: Option<f64>,
    /// GATK's variant quality score.
    pub vqslod: Option<f64>,
}

impl Quality {
    /// The site criteria this record fails, as fixed keywords. They are
    /// the criteria the 1000 Genomes consortium applied before phasing
    /// (panel README, 2022), each used where the dataset publishes the
    /// measure: FILTER is PASS; fewer than 5% of genotypes missing (`full`
    /// = chromosomes of a fully called site, for a dataset that examines
    /// everyone everywhere); Hardy-Weinberg p > 1e-10 in at least one
    /// group; Mendelian error rate at most 5%.
    pub fn failures(&self, examined: u32, full: Option<u32>) -> Vec<String> {
        let mut failed = Vec::new();
        if !matches!(self.filter.as_str(), "PASS" | "." | "") {
            failed.push(format!("FILTER={}", self.filter));
        }
        if full.is_some_and(|full| (examined as f64) < 0.95 * full as f64) {
            failed.push("missingness>5%".to_string());
        }
        if self.hwe.iter().any(|p| !p.is_nan()) && self.hwe.iter().all(|p| p.is_nan() || *p <= 1e-10) {
            failed.push("HWE".to_string());
        }
        if self.mendel.is_some_and(|rate| rate > 0.05) {
            failed.push("mendelian_errors>5%".to_string());
        }
        failed
    }
}

const fn group(code: &'static str, title: &'static str, ancestry: Ancestry) -> Group {
    Group { code, title, ancestry }
}

/// A group whose code is its title.
const fn named(title: &'static str, ancestry: Ancestry) -> Group {
    group(title, title, ancestry)
}

/// The five superpopulations of 1000 Genomes.
pub const KGP_GROUPS: [Group; 5] = [
    group("AFR", "African", Ancestry::African),
    group("AMR", "American", Ancestry::American),
    group("EAS", "East Asian", Ancestry::EastAsian),
    group("EUR", "European", Ancestry::European),
    group("SAS", "South Asian", Ancestry::SouthAsian),
];

/// The genetic-ancestry groups of gnomAD genomes ("remaining" is left
/// out: it is not a population). gnomAD v2.1 has six of them; the others
/// have no chromosomes there and drop out by themselves.
pub const GNOMAD_GROUPS: [Group; 9] = [
    group("afr", "African / African American", Ancestry::African),
    group("ami", "Amish", Ancestry::European),
    group("amr", "Admixed American", Ancestry::American),
    group("asj", "Ashkenazi Jewish", Ancestry::European),
    group("eas", "East Asian", Ancestry::EastAsian),
    group("fin", "Finnish", Ancestry::European),
    group("mid", "Middle Eastern", Ancestry::MiddleEastern),
    group("nfe", "European (non-Finnish)", Ancestry::European),
    group("sas", "South Asian", Ancestry::SouthAsian),
];

/// The populations of ALFA, in the order of its sample columns. The three
/// columns after them are sums of these (African, Asian, total).
pub const ALFA_GROUPS: [Group; 9] = [
    named("European", Ancestry::European),
    named("African Others", Ancestry::African),
    named("East Asian", Ancestry::EastAsian),
    named("African American", Ancestry::African),
    named("Latin American 1", Ancestry::American),
    named("Latin American 2", Ancestry::American),
    named("Other Asian", Ancestry::Other),
    named("South Asian", Ancestry::SouthAsian),
    named("Other", Ancestry::Other),
];

pub const GENOME_ASIA_GROUPS: [Group; 7] = [
    named("Northeast Asian", Ancestry::EastAsian),
    named("Southeast Asian", Ancestry::EastAsian),
    named("South Asian", Ancestry::SouthAsian),
    named("Oceanian", Ancestry::Other),
    named("African", Ancestry::African),
    named("American", Ancestry::American),
    named("West Eurasian", Ancestry::European),
];

pub const TOMMO_GROUPS: [Group; 1] = [named("Japanese", Ancestry::EastAsian)];

pub const WBBC_GROUPS: [Group; 4] = [
    named("North China", Ancestry::EastAsian),
    named("Central China", Ancestry::EastAsian),
    named("South China", Ancestry::EastAsian),
    named("Lingnan", Ancestry::EastAsian),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failures_use_only_the_measures_a_dataset_publishes() {
        let clean = Quality {
            filter: "PASS".to_string(),
            ..Quality::default()
        };
        assert!(clean.failures(5008, Some(5008)).is_empty());
        assert_eq!(clean.failures(4000, Some(5008)), ["missingness>5%"]);
        // Without a full size there is no missingness to speak of.
        assert!(clean.failures(4000, None).is_empty());
        let poor = Quality {
            filter: "AS_VQSR".to_string(),
            hwe: vec![1e-12, f64::NAN, 1e-30],
            mendel: Some(0.2),
            ..Quality::default()
        };
        assert_eq!(poor.failures(10, None), ["FILTER=AS_VQSR", "HWE", "mendelian_errors>5%"]);
        // One group in equilibrium is enough.
        let mixed = Quality {
            hwe: vec![1e-12, 0.3],
            ..Quality::default()
        };
        assert!(mixed.failures(10, None).is_empty());
    }
}

//! Population evidence for one variant from the 1000 Genomes 30x site-level
//! callset (NYGC, 3,202 samples). Counts come from the `*_unrel` INFO
//! fields, i.e. the 2,504 unrelated samples only, per superpopulation.
//! Bounds are exact (Clopper-Pearson) one-sided 95% bounds; tests are exact
//! binomial tests.

use statrs::distribution::{Beta, ChiSquared, ContinuousCDF};
use std::collections::HashMap;

pub const POPS: [&str; 5] = ["AFR", "AMR", "EAS", "EUR", "SAS"];

/// (1-based position, REF, ALT) after trimming shared bases.
pub type Key = (u64, Vec<u8>, Vec<u8>);

#[derive(Clone, Copy, Debug, Default)]
pub struct PopCount {
    pub ac: u32,
    pub an: u32,
    pub het: u32,
    pub hom: u32,
}

impl PopCount {
    pub fn af(&self) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else {
            self.ac as f64 / self.an as f64
        }
    }
    /// Exact (Clopper-Pearson) one-sided 95% lower bound of the frequency.
    pub fn lower95(&self) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else if self.ac == 0 {
            0.0
        } else {
            beta_quantile(self.ac as f64, (self.an - self.ac) as f64 + 1.0, 0.05)
        }
    }
    /// Exact (Clopper-Pearson) one-sided 95% upper bound of the frequency.
    pub fn upper95(&self) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else if self.ac >= self.an {
            1.0
        } else {
            beta_quantile(self.ac as f64 + 1.0, (self.an - self.ac) as f64, 0.95)
        }
    }
    /// Exact binomial p-value for "the frequency is below `threshold`":
    /// the probability of seeing at least `ac` copies in `an` chromosomes
    /// if the frequency were exactly `threshold`.
    pub fn tail_p(&self, threshold: f64) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else if self.ac == 0 {
            1.0
        } else {
            Beta::new(self.ac as f64, (self.an - self.ac.min(self.an)) as f64 + 1.0)
                .map(|b| b.cdf(threshold))
                .unwrap_or(f64::NAN)
        }
    }
}

fn beta_quantile(a: f64, b: f64, q: f64) -> f64 {
    Beta::new(a, b).map(|d| d.inverse_cdf(q)).unwrap_or(f64::NAN)
}

/// Likelihood-ratio (G) test of equal frequency across populations; NaN
/// when the allele is absent or fixed.
pub fn heterogeneity_p(pops: &[PopCount]) -> f64 {
    let (ac, an): (u32, u32) = pops.iter().fold((0, 0), |t, c| (t.0 + c.ac, t.1 + c.an));
    if ac == 0 || ac == an {
        return f64::NAN;
    }
    let p = ac as f64 / an as f64;
    let term = |o: f64, e: f64| if o > 0.0 { o * (o / e).ln() } else { 0.0 };
    let used: Vec<&PopCount> = pops.iter().filter(|c| c.an > 0).collect();
    let g: f64 = used
        .iter()
        .map(|c| {
            let (ac, an) = (c.ac as f64, c.an as f64);
            2.0 * (term(ac, an * p) + term(an - ac, an * (1.0 - p)))
        })
        .sum();
    match ChiSquared::new((used.len().max(2) - 1) as f64) {
        Ok(d) => d.sf(g.max(0.0)),
        Err(_) => f64::NAN,
    }
}

#[derive(Clone, Debug)]
pub struct Site {
    pub filter: String,
    pub vcf_ref: String,
    pub vcf_alt: String,
    pub pops: [PopCount; 5],
    /// Hardy-Weinberg exact-test p-value per superpopulation.
    pub hwe: [f64; 5],
    /// Excess-heterozygosity p-value over all samples (small = excess).
    pub exc_het: f64,
    /// Mendelian error rate in the trios.
    pub mendel: f64,
    pub vqslod: f64,
}

impl Site {
    pub fn total(&self) -> PopCount {
        self.pops.iter().fold(PopCount::default(), |t, p| PopCount {
            ac: t.ac + p.ac,
            an: t.an + p.an,
            het: t.het + p.het,
            hom: t.hom + p.hom,
        })
    }

    /// See [`heterogeneity_p`].
    pub fn heterogeneity_p(&self) -> f64 {
        heterogeneity_p(&self.pops)
    }

    /// The site filters the 1000 Genomes consortium applied before phasing
    /// (panel README, 2022): PASS, genotype missingness < 5%, HWE p > 1e-10
    /// in at least one superpopulation, Mendelian error rate <= 5%.
    /// Returns the criteria this site fails.
    pub fn qc_failures(&self, full_an: u32) -> Vec<String> {
        let mut failed = Vec::new();
        if self.filter != "PASS" {
            failed.push(format!("FILTER={}", self.filter));
        }
        if (self.total().an as f64) < 0.95 * full_an as f64 {
            failed.push("missingness>5%".into());
        }
        if self.hwe.iter().all(|p| p.is_nan() || *p <= 1e-10) && self.hwe.iter().any(|p| !p.is_nan()) {
            failed.push("HWE".into());
        }
        if self.mendel > 0.05 {
            failed.push("mendelian_errors>5%".into());
        }
        failed
    }

    /// Largest upper bound over the superpopulations.
    pub fn max_upper95(&self) -> f64 {
        max_upper95(&self.pops)
    }

    /// (population, counts) with the highest lower bound.
    pub fn strongest_pop(&self) -> (&'static str, PopCount) {
        let i = strongest(&self.pops);
        (POPS[i], self.pops[i])
    }
}

/// Index of the stratum with the highest lower bound.
pub fn strongest(strata: &[PopCount]) -> usize {
    (0..strata.len())
        .max_by(|a, b| {
            let (x, y) = (strata[*a].lower95(), strata[*b].lower95());
            x.partial_cmp(&y).unwrap_or(if x.is_nan() {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            })
        })
        .unwrap_or(0)
}

/// Largest upper bound over the strata.
pub fn max_upper95(strata: &[PopCount]) -> f64 {
    strata.iter().map(|p| p.upper95()).fold(f64::NAN, f64::max)
}

/// Genetic-ancestry groups reported by gnomAD v4 genomes ("remaining" is
/// left out: it is not a population).
pub const GNOMAD_GROUPS: [&str; 9] = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"];

#[derive(Clone, Debug)]
pub struct GnomadSite {
    pub filter: String,
    /// Region flags set on the record: lcr, segdup.
    pub flags: String,
    /// All genomes; `hom` is the number of homozygous individuals.
    pub total: PopCount,
    pub groups: [PopCount; 9],
}

/// One gnomAD sites line (one ALT per line in v4; any further ALT is split).
pub fn parse_gnomad_line(line: &str) -> Vec<(Key, GnomadSite)> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return Vec::new();
    }
    let Ok(pos) = f[1].parse::<u64>() else {
        return Vec::new();
    };
    let mut info: HashMap<&str, &str> = HashMap::new();
    let mut flags = Vec::new();
    for item in f[7].split(';') {
        match item.split_once('=') {
            Some((k, v)) => {
                info.insert(k, v);
            }
            None if item == "lcr" || item == "segdup" => flags.push(item),
            None => {}
        }
    }
    let int = |key: String, i: usize| -> u32 {
        info.get(key.as_str())
            .and_then(|v| v.split(',').nth(i).or(v.split(',').next()))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    };
    let mut out = Vec::new();
    for (i, alt) in f[4].split(',').enumerate() {
        if alt == "*" || alt.starts_with('<') {
            continue;
        }
        let mut groups = [PopCount::default(); 9];
        for (k, g) in GNOMAD_GROUPS.iter().enumerate() {
            groups[k] = PopCount {
                ac: int(format!("AC_{g}"), i),
                an: int(format!("AN_{g}"), 0),
                het: 0,
                hom: int(format!("nhomalt_{g}"), i),
            };
        }
        let site = GnomadSite {
            filter: f[6].to_string(),
            flags: if flags.is_empty() {
                ".".into()
            } else {
                let mut f = flags.clone();
                f.sort_unstable();
                f.join(",")
            },
            total: PopCount {
                ac: int("AC".into(), i),
                an: int("AN".into(), 0),
                het: 0,
                hom: int("nhomalt".into(), i),
            },
            groups,
        };
        out.push((trim(pos, f[3].as_bytes(), alt.as_bytes()), site));
    }
    out
}

/// The FREQ field of a dbSNP record: `dataset:ref,alt1,alt2|dataset:...`.
/// One vector per dataset, indexed like [REF, ALT1, ALT2, ...]; `None` where
/// dbSNP prints `.` (the dataset did not report that allele).
pub fn parse_dbsnp_freq(info: &str) -> Vec<(String, Vec<Option<f64>>)> {
    let Some(field) = info.split(';').find_map(|kv| kv.strip_prefix("FREQ=")) else {
        return Vec::new();
    };
    field
        .split('|')
        .filter_map(|part| {
            let (dataset, values) = part.split_once(':')?;
            Some((
                dataset.to_string(),
                values.split(',').map(|v| v.trim().parse::<f64>().ok()).collect(),
            ))
        })
        .collect()
}

/// Trim bases shared by REF and ALT: suffix first, then prefix, keeping one.
pub fn trim(mut pos: u64, reference: &[u8], alternate: &[u8]) -> Key {
    let (mut r, mut a) = (reference.to_vec(), alternate.to_vec());
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    (pos, r, a)
}

/// One VCF line -> one entry per ALT allele (spanning-deletion `*` skipped).
pub fn parse_line(line: &str) -> Vec<(Key, Site)> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return Vec::new();
    }
    let Ok(pos) = f[1].parse::<u64>() else {
        return Vec::new();
    };
    let info: HashMap<&str, &str> = f[7].split(';').filter_map(|kv| kv.split_once('=')).collect();
    let alts: Vec<&str> = f[4].split(',').collect();
    let pick = |key: &str, i: usize| -> Option<&str> {
        let v = info.get(key)?;
        let parts: Vec<&str> = v.split(',').collect();
        parts.get(i).or(parts.first()).copied()
    };
    let int = |key: String, i: usize| pick(&key, i).and_then(|v| v.parse::<u32>().ok()).unwrap_or(0);
    let float = |key: String, i: usize| pick(&key, i).and_then(|v| v.parse::<f64>().ok()).unwrap_or(f64::NAN);
    let mut out = Vec::new();
    for (i, alt) in alts.iter().enumerate() {
        if *alt == "*" || alt.starts_with('<') {
            continue;
        }
        let mut pops = [PopCount::default(); 5];
        let mut hwe = [f64::NAN; 5];
        for (k, pop) in POPS.iter().enumerate() {
            pops[k] = PopCount {
                ac: int(format!("AC_{pop}_unrel"), i),
                an: int(format!("AN_{pop}_unrel"), 0),
                het: int(format!("AC_Het_{pop}_unrel"), i),
                hom: int(format!("AC_Hom_{pop}_unrel"), i) / 2,
            };
            hwe[k] = float(format!("HWE_{pop}"), i);
        }
        let site = Site {
            filter: f[6].to_string(),
            vcf_ref: f[3].to_string(),
            vcf_alt: f[4].to_string(),
            pops,
            hwe,
            exc_het: float("ExcHet".into(), i),
            mendel: float("ME".into(), 0),
            vqslod: float("VQSLOD".into(), 0),
        };
        out.push((trim(pos, f[3].as_bytes(), alt.as_bytes()), site));
    }
    out
}

/// The two releases of 1000 Genomes that igpop reads.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum KgpRelease {
    /// 30x call set on GRCh38 (NYGC): counts per superpopulation and site
    /// quality in INFO.
    HighCoverage,
    /// Phase 3 on GRCh37: 2,504 unrelated samples, every genotype called.
    Phase3,
}

impl KgpRelease {
    pub fn parse(self, line: &str) -> Vec<(Key, Site)> {
        match self {
            KgpRelease::HighCoverage => parse_line(line),
            KgpRelease::Phase3 => parse_phase3_line(line),
        }
    }
}

/// Chromosomes per superpopulation in phase 3, in the order of `POPS`.
/// Phase 3 has no missing genotypes, so these are the denominators.
pub const PHASE3_CHROMOSOMES: [u32; 5] = [1322, 694, 1008, 1006, 978];

/// One line of a phase 3 genotype file. INFO gives a frequency per
/// superpopulation to four decimals, which determines the count exactly
/// because no superpopulation has more than 1,322 chromosomes. Phase 3
/// publishes no site-quality measures.
pub fn parse_phase3_line(line: &str) -> Vec<(Key, Site)> {
    let f: Vec<&str> = line.splitn(9, '\t').collect();
    if f.len() < 8 {
        return Vec::new();
    }
    let Ok(pos) = f[1].parse::<u64>() else {
        return Vec::new();
    };
    let info: HashMap<&str, &str> = f[7].split(';').filter_map(|kv| kv.split_once('=')).collect();
    let mut out = Vec::new();
    for (i, alt) in f[4].split(',').enumerate() {
        if alt == "*" || alt.starts_with('<') {
            continue;
        }
        let mut pops = [PopCount::default(); 5];
        for (k, pop) in POPS.iter().enumerate() {
            let af = info
                .get(format!("{pop}_AF").as_str())
                .and_then(|v| v.split(',').nth(i))
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(0.0);
            pops[k] = PopCount {
                ac: (af * PHASE3_CHROMOSOMES[k] as f64).round() as u32,
                an: PHASE3_CHROMOSOMES[k],
                het: 0,
                hom: 0,
            };
        }
        let site = Site {
            filter: f[6].to_string(),
            vcf_ref: f[3].to_string(),
            vcf_alt: f[4].to_string(),
            pops,
            hwe: [f64::NAN; 5],
            exc_het: f64::NAN,
            mendel: f64::NAN,
            vqslod: f64::NAN,
        };
        out.push((trim(pos, f[3].as_bytes(), alt.as_bytes()), site));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase3_counts_follow_from_frequencies() {
        let line =
            "14\t106539088\t.\tG\tA\t100\tPASS\tAC=4;AF=0.000798722;AN=5008;EAS_AF=0;AMR_AF=0;AFR_AF=0.0008;EUR_AF=0.001;SAS_AF=0.002\tGT";
        let sites = parse_phase3_line(line);
        let site = &sites[0].1;
        assert_eq!(site.pops.map(|p| p.ac), [1, 0, 0, 1, 2]);
        assert_eq!(site.total().ac, 4);
        assert_eq!(site.total().an, 5008);
        assert!(site.qc_failures(5008).is_empty());
    }

    #[test]
    fn trimming_matches_normalised_form() {
        assert_eq!(trim(100, b"GAC", b"G"), (100, b"GAC".to_vec(), b"G".to_vec()));
        assert_eq!(trim(100, b"GAC", b"GAT"), (102, b"C".to_vec(), b"T".to_vec()));
        assert_eq!(trim(100, b"CAA", b"CA"), (100, b"CA".to_vec(), b"C".to_vec()));
    }

    #[test]
    fn dbsnp_freq_is_indexed_by_allele() {
        let f = parse_dbsnp_freq("RS=1;FREQ=TOPMED:0.93,.,0.066|KOREAN:0.99,0,0.0068;COMMON");
        assert_eq!(f.len(), 2);
        assert_eq!(f[0].0, "TOPMED");
        assert_eq!(f[0].1, vec![Some(0.93), None, Some(0.066)]);
        assert_eq!(f[1].1[1], Some(0.0));
        assert!(parse_dbsnp_freq("RS=1;VC=SNV").is_empty());
    }

    #[test]
    fn exact_test_adapts_to_sample_size() {
        let c = |ac, an| PopCount { ac, an, het: 0, hom: 0 };
        // 1.5% observed: three copies in 200 is not evidence for >= 1%, 300 in 20,000 is.
        assert!(c(3, 200).tail_p(0.01) > 0.05);
        assert!(c(300, 20_000).tail_p(0.01) < 1e-6);
        assert!(c(3, 200).lower95() < 0.01 && c(300, 20_000).lower95() > 0.01);
        assert_eq!(c(0, 500).lower95(), 0.0);
        assert!((c(0, 5008).upper95() - 0.000598).abs() < 1e-5);
    }

    #[test]
    fn parses_gnomad_groups_and_flags() {
        let line = "chr14\t100\t.\tG\tT\t.\tAS_VQSR\tAC=7;AN=150000;nhomalt=1;AC_afr=5;AN_afr=40000;AC_nfe=2;AN_nfe=60000;segdup;lcr";
        let sites = parse_gnomad_line(line);
        assert_eq!(sites.len(), 1);
        let site = &sites[0].1;
        assert_eq!((site.total.ac, site.total.an, site.total.hom), (7, 150000, 1));
        assert_eq!(site.groups[0].ac, 5);
        assert_eq!(site.groups[7].an, 60000);
        assert_eq!(site.flags, "lcr,segdup");
        assert_eq!(site.filter, "AS_VQSR");
    }

    #[test]
    fn parses_multiallelic_unrelated_counts() {
        let line =
            "chr14\t100\t.\tG\tC,T\t50\tPASS\tAN_AFR_unrel=1000;AC_AFR_unrel=5,7;AN_EUR_unrel=800;AC_EUR_unrel=0,1;HWE_AFR=1,0.5;ME=0.01";
        let sites = parse_line(line);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[1].0, (100, b"G".to_vec(), b"T".to_vec()));
        assert_eq!(sites[1].1.pops[0].ac, 7);
        assert_eq!(sites[1].1.pops[3].ac, 1);
        assert_eq!(sites[1].1.total().an, 1800);
    }
}

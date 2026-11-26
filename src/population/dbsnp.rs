//! dbSNP: the rsID of a variant, the COMMON flag, and the frequencies of
//! the studies dbSNP relays.
//!
//! dbSNP passes on a frequency per study without counts. To treat a study
//! like any other population, its size is estimated from the smallest
//! frequency it reports anywhere in the regions read (one copy among all
//! its chromosomes), and the count is frequency times size. This is an
//! estimate, and the output marks such counts as estimated.

use super::Ancestry;
use crate::stats::{informative, PopCount};
use crate::variant::{trim, Lookup};
use std::collections::{HashMap, HashSet};

/// One study relayed by dbSNP, and whether it is used.
pub struct Study {
    pub name: String,
    pub records: u64,
    pub smallest_af: f64,
    /// Estimated size: one copy among all chromosomes gives the smallest
    /// non-zero frequency a dataset can report.
    pub chromosomes: u32,
    pub used: bool,
    pub reason: &'static str,
}

/// One allele of a dbSNP record with the frequencies reported for it, as
/// (dataset, frequency).
type Reported = (Lookup, Vec<(String, f64)>);

/// What dbSNP holds for the regions read.
#[derive(Default)]
pub struct Dbsnp {
    pub rsids: HashMap<Lookup, String>,
    /// Frequencies per variant allele, as (dataset, frequency).
    freqs: HashMap<Lookup, Vec<(String, f64)>>,
    /// Alleles carrying the COMMON flag: of a site with several alternate
    /// alleles, only the most frequent one. Decided in `select`.
    pub common: HashSet<Lookup>,
    /// Records with the COMMON flag: their alleles, each with the
    /// frequencies reported for it.
    flagged: Vec<Vec<Reported>>,
    /// Per dataset: number of records, smallest non-zero frequency.
    seen: HashMap<String, (u64, f64)>,
    /// Decided once all lines are read, largest dataset first.
    pub studies: Vec<Study>,
    /// The datasets in use, by name: positions in `datasets`.
    used: Vec<usize>,
}

impl Dbsnp {
    pub fn add_line(&mut self, contig: &str, line: &str) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 || f[2] == "." || f[4] == "." {
            return;
        }
        let Ok(pos) = f[1].parse::<u64>() else { return };
        let reported = f.get(7).map(|info| parse_freq(info)).unwrap_or_default();
        for (dataset, values) in &reported {
            let entry = self.seen.entry(dataset.clone()).or_insert((0, f64::INFINITY));
            entry.0 += 1;
            for v in values.iter().skip(1).flatten().filter(|v| **v > 0.0) {
                entry.1 = entry.1.min(*v);
            }
        }
        let alts: Vec<&str> = f[4].split(',').collect();
        let lookup = |alt: &str| (contig.to_string(), trim(pos, f[3].as_bytes(), alt.as_bytes()));
        let flagged = f.get(7).is_some_and(|info| info.split(';').any(|x| x == "COMMON"));
        let mut alleles = Vec::new();
        for (i, alt) in alts.iter().enumerate() {
            let of_alt: Vec<(String, f64)> = reported
                .iter()
                .filter_map(|(dataset, values)| values.get(i + 1).copied().flatten().map(|v| (dataset.clone(), v)))
                .collect();
            self.rsids.insert(lookup(alt), f[2].to_string());
            if flagged {
                alleles.push((lookup(alt), of_alt.clone()));
            }
            self.freqs.insert(lookup(alt), of_alt);
        }
        if flagged {
            self.flagged.push(alleles);
        }
    }

    /// Decide which datasets are used. `duplicated`: datasets that re-report
    /// a source read directly in this run.
    pub fn select(&mut self, excluded: &[String], duplicated: &[&str], common: f64) {
        self.studies = self
            .seen
            .iter()
            .map(|(name, (records, smallest_af))| {
                let chromosomes = if smallest_af.is_finite() {
                    (1.0 / smallest_af).round() as u32
                } else {
                    0
                };
                let reason = if excluded.contains(name) {
                    "excluded by --dbsnp-exclude"
                } else if duplicated.contains(&name.as_str()) {
                    "re-report of a source read directly"
                } else if !informative(chromosomes, common) {
                    "too small for --common"
                } else {
                    "."
                };
                Study {
                    name: name.clone(),
                    records: *records,
                    smallest_af: *smallest_af,
                    chromosomes,
                    used: reason == ".",
                    reason,
                }
            })
            .collect();
        self.studies
            .sort_by(|a, b| a.smallest_af.total_cmp(&b.smallest_af).then(a.name.cmp(&b.name)));
        self.used = (0..self.studies.len()).filter(|i| self.studies[*i].used).collect();
        self.used.sort_by(|a, b| self.studies[*a].name.cmp(&self.studies[*b].name));
        // dbSNP flags a record, not an allele. The flag is given to the
        // allele of the record with the highest frequency in any dataset
        // large enough to show a frequency of `common`: a dataset of a few
        // people reports half of everything it sees as frequent.
        let large: HashSet<&str> = self
            .studies
            .iter()
            .filter(|d| informative(d.chromosomes, common) && !excluded.contains(&d.name))
            .map(|d| d.name.as_str())
            .collect();
        let highest = |freqs: &[(String, f64)]| {
            freqs
                .iter()
                .filter(|(dataset, _)| large.contains(dataset.as_str()))
                .map(|(_, af)| *af)
                .fold(0.0, f64::max)
        };
        self.common = self
            .flagged
            .iter()
            .filter_map(|alleles| {
                alleles
                    .iter()
                    .map(|(lookup, freqs)| (highest(freqs), lookup))
                    .filter(|(af, _)| *af > 0.0)
                    .max_by(|a, b| a.0.total_cmp(&b.0))
                    .map(|(_, lookup)| lookup.clone())
            })
            .collect();
    }

    /// The studies in use, in the order their strata are numbered.
    pub fn used(&self) -> impl Iterator<Item = &Study> {
        self.used.iter().map(|i| &self.studies[*i])
    }

    /// Non-zero frequencies of the used datasets for one allele, highest first.
    pub fn relayed(&self, lookup: &Lookup) -> Vec<(String, f64)> {
        let mut out: Vec<(String, f64)> = self
            .freqs
            .get(lookup)
            .map(|all| {
                all.iter()
                    .filter(|(d, af)| *af > 0.0 && self.used().any(|u| &u.name == d))
                    .cloned()
                    .collect()
            })
            .unwrap_or_default();
        out.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        out
    }

    /// The estimated counts of this allele in every used study that
    /// reports it, as (position among the used studies, counts).
    pub fn counts<'a>(&'a self, lookup: &Lookup) -> impl Iterator<Item = (usize, PopCount)> + 'a {
        let reported = self.freqs.get(lookup);
        self.used().enumerate().filter_map(move |(k, study)| {
            let (_, af) = reported?.iter().find(|(name, _)| name == &study.name)?;
            Some((k, PopCount::new((af * study.chromosomes as f64).round() as u32, study.chromosomes)))
        })
    }
}

impl Study {
    /// The continental group of the studies known to be of one.
    pub fn ancestry(&self) -> Ancestry {
        match self.name.as_str() {
            "KOREAN" | "Korea1K" | "Korea4K" | "TOMMO" | "Vietnamese" => Ancestry::EastAsian,
            "GoNL" | "TWINSUK" | "ALSPAC" | "Estonian" | "NorthernSweden" | "FINRISK" | "GENOME_DK" => Ancestry::European,
            "Qatari" => Ancestry::MiddleEastern,
            _ => Ancestry::Other,
        }
    }
}

/// The FREQ field of a dbSNP record: `dataset:ref,alt1,alt2|dataset:...`.
/// One vector per dataset, indexed like [REF, ALT1, ALT2, ...]; `None` where
/// dbSNP prints `.` (the dataset did not report that allele).
fn parse_freq(info: &str) -> Vec<(String, Vec<Option<f64>>)> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn freq_is_indexed_by_allele() {
        let f = parse_freq("RS=1;FREQ=TOPMED:0.93,.,0.066|KOREAN:0.99,0,0.0068;COMMON");
        assert_eq!(f.len(), 2);
        assert_eq!(f[0].0, "TOPMED");
        assert_eq!(f[0].1, vec![Some(0.93), None, Some(0.066)]);
        assert_eq!(f[1].1[1], Some(0.0));
        assert!(parse_freq("RS=1;VC=SNV").is_empty());
    }

    /// An rsID and its frequencies belong to one exact change. A base
    /// change recorded at the position of a deletion is another variant,
    /// and of a record with several alternate alleles only the one that
    /// is the same change counts.
    #[test]
    fn dbsnp_matches_the_exact_change_only() {
        let mut db = Dbsnp::default();
        db.add_line("chr14", "NC_000014.9\t106185188\trs1\tG\tA\t.\t.\tRS=1;VC=SNV;FREQ=TOPMED:0.5,0.5");
        db.add_line(
            "chr14",
            "NC_000014.9\t106185188\trs2\tGAAC\tG\t.\t.\tRS=2;VC=DEL;FREQ=TOPMED:0.99,0.01",
        );
        // The repeat written out in full, as dbSNP does: three alleles.
        db.add_line(
            "chr14",
            "NC_000014.9\t106116762\trs3\tAACTACT\tA,AACT,AACTACTACT\t.\t.\tRS=3;VC=INDEL;FREQ=TOPMED:0.9,0.01,0.02,0.07",
        );
        let at = |pos: u64, r: &str, a: &str| ("chr14".to_string(), (pos, r.as_bytes().to_vec(), a.as_bytes().to_vec()));
        let rs = |pos, r, a| db.rsids.get(&at(pos, r, a)).map(String::as_str);
        // The deletion of AAC after 106185188 is rs2, never the base change rs1.
        assert_eq!(rs(106185188, "GAAC", "G"), Some("rs2"));
        assert_eq!(rs(106185188, "G", "A"), Some("rs1"));
        assert_eq!(db.freqs[&at(106185188, "GAAC", "G")], vec![("TOPMED".to_string(), 0.01)]);
        // Another deletion at the same place has no record at all.
        assert_eq!(rs(106185188, "GA", "G"), None);
        // One repeat unit more, one less, two less: each its own allele.
        assert_eq!(db.freqs[&at(106116762, "A", "AACT")], vec![("TOPMED".to_string(), 0.07)]);
        assert_eq!(db.freqs[&at(106116762, "AACT", "A")], vec![("TOPMED".to_string(), 0.02)]);
        assert_eq!(db.freqs[&at(106116762, "AACTACT", "A")], vec![("TOPMED".to_string(), 0.01)]);
        assert_eq!(rs(106116762, "A", "C"), None);
    }

    /// dbSNP flags a record, not an allele. The flag goes to the allele
    /// that is most frequent in a dataset of adequate size, not to one
    /// that a study of a few people happens to report at one half.
    #[test]
    fn common_flag_goes_to_the_frequent_allele() {
        let mut db = Dbsnp::default();
        // Sizes of the studies, as other records reveal them.
        db.add_line("chr14", "NC_000014.9\t1\trs9\tA\tC\t.\t.\tRS=9;FREQ=1000Genomes_30X:0.9998,0.000156|TOMMO:0.9999,0.0000129|GoESP:0.9999,0.0000772|GnomAD_genomes:0.9999,0.0000067");
        db.add_line(
            "chr14",
            "NC_000014.9\t106714751\trs1\tG\tA,C\t.\t.\tRS=1;FREQ=1000Genomes_30X:0.9861,0.0139,.|SGDP_PRJ:0.5,.,0.5|TOMMO:0.9553,0.04469,.;COMMON",
        );
        // The frequent allele has no 1000 Genomes frequency at all.
        db.add_line(
            "chr14",
            "NC_000014.9\t106770643\trs2\tG\tA,C,T\t.\t.\tRS=2;FREQ=1000Genomes_30X:0.9948,0.005192,.,.|GnomAD_genomes:0.9489,.,0.05109,.|TOMMO:0.9364,.,0.06,.;COMMON",
        );
        db.add_line(
            "chr14",
            "NC_000014.9\t106770746\trs3\tC\tA,G,T\t.\t.\tRS=3;FREQ=1000Genomes_30X:0.9992,.,.,0.0007987|GoESP:0.9913,0.008701,.,.|SGDP_PRJ:0.5,.,.,0.5|TOMMO:0.8824,.,.,0.1176;COMMON",
        );
        db.select(&[], &[], 0.01);
        let common = |pos: u64, r: &str, a: &str| {
            db.common
                .contains(&("chr14".to_string(), (pos, r.as_bytes().to_vec(), a.as_bytes().to_vec())))
        };
        assert!(common(106714751, "G", "A") && !common(106714751, "G", "C"));
        assert!(common(106770643, "G", "C") && !common(106770643, "G", "A"));
        assert!(common(106770746, "C", "T") && !common(106770746, "C", "A"));
    }
}

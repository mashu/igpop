//! Where a dataset's VCF keeps the counts of each group. One line goes in,
//! one `Record` per alternate allele comes out, keyed by the trimmed
//! variant so that it compares directly with a called variant.

use super::{Quality, Record, GNOMAD_GROUPS, KGP_GROUPS};
use crate::stats::PopCount;
use crate::variant::{trim, Key};
use std::collections::HashMap;

/// Chromosomes per superpopulation in 1000 Genomes phase 3, in the order
/// of `KGP_GROUPS`. Phase 3 has no missing genotypes, so these are the
/// denominators everywhere.
pub const PHASE3_CHROMOSOMES: [u32; 5] = [1322, 694, 1008, 1006, 978];

/// The layout of a dataset's records.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// 1000 Genomes 30x (NYGC): INFO `AC_<POP>_unrel` / `AN_<POP>_unrel`,
    /// the unrelated samples only, with genotype counts and site quality.
    Kgp30x,
    /// 1000 Genomes phase 3: INFO `<POP>_AF` to four decimals, which gives
    /// the count exactly because no superpopulation has more than 1,322
    /// chromosomes. No site-quality measures.
    KgpPhase3,
    /// gnomAD genomes: INFO `AC_<group>` / `AN_<group>` / `nhomalt`, the
    /// totals over all genomes, and the region flags.
    Gnomad,
    /// One sample column per population, each `AN:AC[,AC...]:...` (ALFA).
    /// `columns`: which sample columns are the groups, in their order.
    Samples { columns: &'static [usize] },
    /// INFO fields: per group the key of the allele count and of the
    /// number of chromosomes.
    Counts { keys: &'static [(&'static str, &'static str)] },
    /// INFO fields: per group the key of the frequency and of the number
    /// of chromosomes; the count is worked out from the two.
    Frequencies { keys: &'static [(&'static str, &'static str)] },
}

/// The fixed columns of a VCF line and its INFO field by key.
struct Line<'a> {
    pos: u64,
    reference: &'a str,
    alternates: Vec<&'a str>,
    filter: &'a str,
    info: HashMap<&'a str, &'a str>,
    /// INFO entries without a value.
    flags: Vec<&'a str>,
    /// The sample columns, unsplit: a phased file has thousands.
    samples: &'a str,
}

impl<'a> Line<'a> {
    fn parse(line: &'a str) -> Option<Self> {
        let f: Vec<&str> = line.splitn(10, '\t').collect();
        if f.len() < 8 {
            return None;
        }
        let (mut info, mut flags) = (HashMap::new(), Vec::new());
        for item in f[7].split(';') {
            match item.split_once('=') {
                Some((key, value)) => {
                    info.insert(key, value);
                }
                None => flags.push(item),
            }
        }
        Some(Line {
            pos: f[1].parse().ok()?,
            reference: f[3],
            alternates: f[4].split(',').collect(),
            filter: f[6],
            info,
            flags,
            samples: f.get(9).copied().unwrap_or(""),
        })
    }

    /// The value of `key` for alternate allele `i`; a field with one value
    /// holds for every allele.
    fn of_allele(&self, key: &str, i: usize) -> Option<&'a str> {
        let value = self.info.get(key)?;
        value.split(',').nth(i).or(value.split(',').next())
    }

    fn int(&self, key: &str, i: usize) -> u32 {
        self.of_allele(key, i).and_then(|v| v.parse().ok()).unwrap_or(0)
    }

    fn float(&self, key: &str, i: usize) -> Option<f64> {
        self.of_allele(key, i).and_then(|v| v.parse().ok())
    }

    /// The value of `key` for alternate allele `i` of a field that must
    /// list every allele: nothing if it does not.
    fn listed<T: std::str::FromStr>(&self, key: &str, i: usize) -> Option<T> {
        self.info.get(key)?.split(',').nth(i)?.parse().ok()
    }

    /// One record per alternate allele that is a sequence (not `*`, `.`
    /// or a symbolic allele), built by `record` from the allele's index.
    fn records(&self, record: impl Fn(usize) -> Record) -> Vec<(Key, Record)> {
        self.alternates
            .iter()
            .enumerate()
            .filter(|(_, alt)| !matches!(**alt, "*" | ".") && !alt.starts_with('<'))
            .map(|(i, alt)| (trim(self.pos, self.reference.as_bytes(), alt.as_bytes()), record(i)))
            .collect()
    }
}

impl Layout {
    /// One VCF line -> per alternate allele its key and its record.
    pub fn parse(&self, line: &str) -> Vec<(Key, Record)> {
        let Some(line) = Line::parse(line) else { return Vec::new() };
        let filter = line.filter.to_string();
        match *self {
            Layout::Kgp30x => line.records(|i| {
                let of = |field: &str, pop: &str| line.int(&format!("{field}_{pop}_unrel"), i);
                let people = |field: &str| KGP_GROUPS.iter().map(|g| of(field, g.code)).sum::<u32>();
                Record {
                    counts: KGP_GROUPS
                        .iter()
                        .map(|g| PopCount::new(of("AC", g.code), line.int(&format!("AN_{}_unrel", g.code), 0)))
                        .collect(),
                    quality: Quality {
                        filter: filter.clone(),
                        heterozygous: Some(people("AC_Het")),
                        // AC_Hom counts chromosomes: two per person.
                        homozygous: Some(KGP_GROUPS.iter().map(|g| of("AC_Hom", g.code) / 2).sum()),
                        hwe: KGP_GROUPS
                            .iter()
                            .map(|g| line.float(&format!("HWE_{}", g.code), i).unwrap_or(f64::NAN))
                            .collect(),
                        excess_het: line.float("ExcHet", i),
                        mendel: line.float("ME", 0),
                        vqslod: line.float("VQSLOD", 0),
                        ..Quality::default()
                    },
                }
            }),
            Layout::KgpPhase3 => line.records(|i| Record {
                counts: KGP_GROUPS
                    .iter()
                    .zip(PHASE3_CHROMOSOMES)
                    .map(|(g, an)| {
                        let af: f64 = line.listed(&format!("{}_AF", g.code), i).unwrap_or(0.0);
                        PopCount::new((af * an as f64).round() as u32, an)
                    })
                    .collect(),
                quality: Quality {
                    filter: filter.clone(),
                    ..Quality::default()
                },
            }),
            Layout::Gnomad => {
                let mut flags: Vec<&'static str> = ["lcr", "segdup"].into_iter().filter(|flag| line.flags.contains(flag)).collect();
                flags.sort_unstable();
                line.records(|i| Record {
                    counts: GNOMAD_GROUPS
                        .iter()
                        .map(|g| PopCount::new(line.int(&format!("AC_{}", g.code), i), line.int(&format!("AN_{}", g.code), 0)))
                        .collect(),
                    quality: Quality {
                        filter: filter.clone(),
                        flags: flags.clone(),
                        all: Some(PopCount::new(line.int("AC", i), line.int("AN", 0))),
                        homozygous: Some(line.int("nhomalt", i)),
                        ..Quality::default()
                    },
                })
            }
            Layout::Samples { columns } => {
                let samples: Vec<&str> = line.samples.split('\t').collect();
                line.records(|i| Record {
                    counts: columns
                        .iter()
                        .map(|column| {
                            let mut parts = samples.get(*column).copied().unwrap_or("").split(':');
                            let an: u32 = parts.next().and_then(|x| x.parse().ok()).unwrap_or(0);
                            let ac: u32 = parts
                                .next()
                                .and_then(|x| x.split(',').nth(i))
                                .and_then(|x| x.parse().ok())
                                .unwrap_or(0);
                            PopCount::new(ac.min(an), an)
                        })
                        .collect(),
                    quality: Quality::default(),
                })
            }
            Layout::Counts { keys } => line.records(|i| Record {
                counts: keys
                    .iter()
                    .map(|(count, size)| {
                        let an = line.int(size, 0);
                        PopCount::new(line.listed::<u32>(count, i).unwrap_or(0).min(an), an)
                    })
                    .collect(),
                quality: Quality::default(),
            }),
            Layout::Frequencies { keys } => line.records(|i| Record {
                counts: keys
                    .iter()
                    .map(|(frequency, size)| {
                        let an = line.int(size, 0);
                        let af: f64 = line.listed(frequency, i).unwrap_or(0.0);
                        PopCount::new(((af * an as f64).round() as u32).min(an), an)
                    })
                    .collect(),
                quality: Quality::default(),
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::total;

    const ALFA: Layout = Layout::Samples {
        columns: &[0, 1, 2, 3, 4, 5, 6, 7, 8],
    };

    #[test]
    fn phase3_counts_follow_from_frequencies() {
        let line =
            "14\t106539088\t.\tG\tA\t100\tPASS\tAC=4;AF=0.000798722;AN=5008;EAS_AF=0;AMR_AF=0;AFR_AF=0.0008;EUR_AF=0.001;SAS_AF=0.002\tGT";
        let record = &Layout::KgpPhase3.parse(line)[0].1;
        assert_eq!(record.counts.iter().map(|c| c.ac).collect::<Vec<_>>(), [1, 0, 0, 1, 2]);
        assert_eq!(record.all(), PopCount::new(4, 5008));
        assert!(record.quality.failures(5008, Some(5008)).is_empty());
    }

    #[test]
    fn gnomad_groups_totals_and_flags() {
        let line = "chr14\t100\t.\tG\tT\t.\tAS_VQSR\tAC=7;AN=150000;nhomalt=1;AC_afr=5;AN_afr=40000;AC_nfe=2;AN_nfe=60000;segdup;lcr";
        let sites = Layout::Gnomad.parse(line);
        assert_eq!(sites.len(), 1);
        let record = &sites[0].1;
        // The total is the published one, not the sum of the groups.
        assert_eq!(record.all(), PopCount::new(7, 150000));
        assert_eq!(record.quality.homozygous, Some(1));
        assert_eq!(record.counts[0].ac, 5);
        assert_eq!(record.counts[7].an, 60000);
        assert_eq!(record.quality.flags, ["lcr", "segdup"]);
        assert_eq!(record.quality.filter, "AS_VQSR");
    }

    #[test]
    fn kgp_counts_are_those_of_the_unrelated_samples_per_allele() {
        let line =
            "chr14\t100\t.\tG\tC,T\t50\tPASS\tAN_AFR_unrel=1000;AC_AFR_unrel=5,7;AN_EUR_unrel=800;AC_EUR_unrel=0,1;HWE_AFR=1,0.5;ME=0.01";
        let sites = Layout::Kgp30x.parse(line);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[1].0, (100, b"G".to_vec(), b"T".to_vec()));
        assert_eq!(sites[1].1.counts[0].ac, 7);
        assert_eq!(sites[1].1.counts[3].ac, 1);
        assert_eq!(total(&sites[1].1.counts).an, 1800);
        assert_eq!(sites[1].1.quality.mendel, Some(0.01));
    }

    #[test]
    fn alfa_counts_come_from_the_population_columns() {
        let line = "NC_000014.9\t100\trs1\tG\tA,C\t.\t.\t.\tAN:AC:HWEP:GR:GV:GA\t9428:93,2:0:1:1:0\t86:0,0:-1:43:0:0\t86:1,0:0:0:0:0\t2352:2,0:0:0:0:0\t146:0,0:0:0:0:0\t610:1,0:0:0:0:0\t24:0,0:0:0:0:0\t98:0,5:0:0:0:0\t640:4,0:0:0:0:0\t2438:2,0:0:0:0:0\t110:1,0:0:0:0:0\t13470:101,7:0:0:0:0";
        let sites = ALFA.parse(line);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].0, (100, b"G".to_vec(), b"A".to_vec()));
        // Nine populations; the three sums at the end are left out.
        assert_eq!(sites[0].1.counts.len(), 9);
        assert_eq!(sites[0].1.counts[0], PopCount::new(93, 9428));
        assert_eq!(sites[1].1.counts[7], PopCount::new(5, 98));
        assert_eq!(sites[0].1.counts[8], PopCount::new(4, 640));
    }

    #[test]
    fn info_counts_and_frequencies() {
        let counts = Layout::Counts {
            keys: &[("AC_NEA", "AN_NEA"), ("AC_SEA", "AN_SEA"), ("AC_OCE", "AN_OCE")],
        };
        let asia = "chr14\t100\t.\tT\tA\t9\tPASS\tAC=473;AN=3230;AN_SEA=600;AC_SEA=12;AN_NEA=800;AC_NEA=0;AN_OCE=0;AC_OCE=0";
        let sites = counts.parse(asia);
        assert_eq!(sites[0].1.counts[1], PopCount::new(12, 600));
        assert_eq!(sites[0].1.counts[2].an, 0);
        let frequencies = Layout::Frequencies {
            keys: &[("North_AF", "North_AN"), ("South_AF", "South_AN")],
        };
        let china = "chr14\t100\t.\tA\tG\t9\tPASS\tAC=2;AN=8960;North_AF=0.0;North_AN=448;South_AF=0.000247831;South_AN=8070";
        let sites = frequencies.parse(china);
        assert_eq!(sites[0].1.counts[1], PopCount::new(2, 8070));
        assert_eq!(sites[0].1.counts[0], PopCount::new(0, 448));
        let one = Layout::Counts { keys: &[("AC", "AN")] };
        let japan = "chr14\t100\trs1\tT\tA,C\t180\tPASS\tAC=1,7;AF=0.000181,0.001;AN=5530";
        assert_eq!(one.parse(japan)[1].1.counts[0], PopCount::new(7, 5530));
        // An allele that is not a sequence has no record.
        assert!(one.parse("chr14\t100\t.\tT\t*\t.\tPASS\tAC=1;AN=10").is_empty());
    }
}

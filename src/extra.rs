//! Further frequency datasets with counts per group, read beside 1000
//! Genomes and gnomAD: NCBI's ALFA with its populations, the GenomeAsia
//! pilot, the Japanese ToMMo panel and the Chinese WBBC.
//!
//! Each is one tabix-indexed VCF. What differs is where the counts of a
//! group stand in a record, which `Layout` describes.

use crate::freq::{trim, Key, PopCount};
use std::collections::HashMap;

/// Where the counts of each group stand in a record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// One sample column per population, each `AN:AC[,AC...]:...`.
    /// `columns`: which sample columns are the groups, in their order.
    Samples { columns: &'static [usize] },
    /// INFO fields: per group the key of the allele count and of the
    /// number of chromosomes.
    Counts { keys: &'static [(&'static str, &'static str)] },
    /// INFO fields: per group the key of the frequency and of the number
    /// of chromosomes; the count is worked out from the two.
    Frequencies { keys: &'static [(&'static str, &'static str)] },
}

/// One further dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Extra {
    /// Short name, as in `support_from` (`ALFA:European`).
    pub name: &'static str,
    pub title: &'static str,
    pub about: &'static str,
    pub people: Option<u32>,
    pub spec: &'static str,
    pub groups: &'static [&'static str],
    pub layout: Layout,
    /// Datasets relayed by dbSNP that report the same people.
    pub replaces: &'static [&'static str],
}

/// The further datasets of GRCh38.
pub const GRCH38: [Extra; 4] = [
    Extra {
        name: "ALFA",
        title: "ALFA (NCBI), by population",
        about: "people of studies held at NCBI, in nine populations; how many were examined differs from position to position",
        people: None,
        spec: "https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz",
        groups: &[
            "European",
            "African Others",
            "East Asian",
            "African American",
            "Latin American 1",
            "Latin American 2",
            "Other Asian",
            "South Asian",
            "Other",
        ],
        // The last three columns are sums of the others: African, Asian, total.
        layout: Layout::Samples {
            columns: &[0, 1, 2, 3, 4, 5, 6, 7, 8],
        },
        replaces: &["dbGaP_PopFreq"],
    },
    Extra {
        name: "GenomeAsia",
        title: "GenomeAsia pilot",
        about: "219 population groups of Asia and a worldwide comparison set, in seven groups including Southeast Asian",
        people: Some(1_739),
        spec: "https://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/ga100k/ga100k.subst.vcf.gz",
        groups: &[
            "Northeast Asian",
            "Southeast Asian",
            "South Asian",
            "Oceanian",
            "African",
            "American",
            "West Eurasian",
        ],
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
        replaces: &[],
    },
    Extra {
        name: "ToMMo_61k",
        title: "ToMMo 61KJPN",
        about: "Japanese people, sequenced as whole genomes",
        people: Some(61_000),
        spec: "https://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/tommo61kjpn/tommo-61kjpn-20250616-GRCh38-snvindel-af-autosome.vcf.gz",
        groups: &["Japanese"],
        layout: Layout::Counts { keys: &[("AC", "AN")] },
        replaces: &["TOMMO"],
    },
    Extra {
        name: "WBBC",
        title: "Westlake BioBank for Chinese",
        about: "Chinese people from four regions of China",
        people: Some(4_480),
        spec: "https://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/wbbc/wbbc.vcf.gz",
        groups: &["North China", "Central China", "South China", "Lingnan"],
        layout: Layout::Frequencies {
            keys: &[
                ("North_AF", "North_AN"),
                ("Central_AF", "Central_AN"),
                ("South_AF", "South_AN"),
                ("Lingnan_AF", "Lingnan_AN"),
            ],
        },
        replaces: &[],
    },
];

impl Extra {
    /// One VCF line -> per alternate allele the counts of every group.
    pub fn parse(&self, line: &str) -> Vec<(Key, Vec<PopCount>)> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Vec::new();
        }
        let Ok(pos) = f[1].parse::<u64>() else {
            return Vec::new();
        };
        let alts: Vec<&str> = f[4].split(',').collect();
        let count = |ac: u32, an: u32| PopCount {
            ac: ac.min(an),
            an,
            het: 0,
            hom: 0,
        };
        // Per group: the chromosomes examined and one value per allele.
        let groups: Vec<(u32, Vec<f64>, bool)> = match self.layout {
            Layout::Samples { columns } => columns
                .iter()
                .map(|c| {
                    let mut parts = f.get(9 + c).copied().unwrap_or("").split(':');
                    let an = parts.next().and_then(|x| x.parse().ok()).unwrap_or(0);
                    let acs = parts.next().unwrap_or("").split(',').map(|x| x.parse().unwrap_or(0.0)).collect();
                    (an, acs, false)
                })
                .collect(),
            Layout::Counts { keys } | Layout::Frequencies { keys } => {
                let info: HashMap<&str, &str> = f[7].split(';').filter_map(|kv| kv.split_once('=')).collect();
                let share = matches!(self.layout, Layout::Frequencies { .. });
                keys.iter()
                    .map(|(value, size)| {
                        let an = info.get(size).and_then(|x| x.parse().ok()).unwrap_or(0);
                        let values = info
                            .get(value)
                            .map(|x| x.split(',').map(|v| v.parse().unwrap_or(0.0)).collect())
                            .unwrap_or_default();
                        (an, values, share)
                    })
                    .collect()
            }
        };
        alts.iter()
            .enumerate()
            .filter(|(_, alt)| !matches!(**alt, "*" | ".") && !alt.starts_with('<'))
            .map(|(i, alt)| {
                let counts = groups
                    .iter()
                    .map(|(an, values, share)| {
                        let value = values.get(i).copied().unwrap_or(0.0);
                        count(
                            if *share {
                                (value * *an as f64).round() as u32
                            } else {
                                value as u32
                            },
                            *an,
                        )
                    })
                    .collect();
                (trim(pos, f[3].as_bytes(), alt.as_bytes()), counts)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alfa_counts_come_from_the_population_columns() {
        let line = "NC_000014.9\t100\trs1\tG\tA,C\t.\t.\t.\tAN:AC:HWEP:GR:GV:GA\t9428:93,2:0:1:1:0\t86:0,0:-1:43:0:0\t86:1,0:0:0:0:0\t2352:2,0:0:0:0:0\t146:0,0:0:0:0:0\t610:1,0:0:0:0:0\t24:0,0:0:0:0:0\t98:0,5:0:0:0:0\t640:4,0:0:0:0:0\t2438:2,0:0:0:0:0\t110:1,0:0:0:0:0\t13470:101,7:0:0:0:0";
        let sites = GRCH38[0].parse(line);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].0, (100, b"G".to_vec(), b"A".to_vec()));
        // Nine populations; the three sums at the end are left out.
        assert_eq!(sites[0].1.len(), 9);
        assert_eq!((sites[0].1[0].ac, sites[0].1[0].an), (93, 9428));
        assert_eq!((sites[1].1[7].ac, sites[1].1[7].an), (5, 98));
        assert_eq!((sites[0].1[8].ac, sites[0].1[8].an), (4, 640));
    }

    #[test]
    fn info_counts_and_frequencies() {
        let asia = "chr14\t100\t.\tT\tA\t9\tPASS\tAC=473;AN=3230;AN_SEA=600;AC_SEA=12;AN_NEA=800;AC_NEA=0;AN_OCE=0;AC_OCE=0";
        let sites = GRCH38[1].parse(asia);
        // Southeast Asian is the second group.
        assert_eq!((sites[0].1[1].ac, sites[0].1[1].an), (12, 600));
        assert_eq!(sites[0].1[3].an, 0);
        let china = "chr14\t100\t.\tA\tG\t9\tPASS\tAC=2;AN=8960;North_AF=0.0;North_AN=448;South_AF=0.000247831;South_AN=8070";
        let sites = GRCH38[3].parse(china);
        assert_eq!((sites[0].1[2].ac, sites[0].1[2].an), (2, 8070));
        assert_eq!((sites[0].1[0].ac, sites[0].1[0].an), (0, 448));
        let japan = "chr14\t100\trs1\tT\tA,C\t180\tPASS\tAC=1,7;AF=0.000181,0.001;AN=5530";
        let sites = GRCH38[2].parse(japan);
        assert_eq!((sites[1].1[0].ac, sites[1].1[0].an), (7, 5530));
    }
}

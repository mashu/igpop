//! Haplotype check against a phased panel: is the allele's combination of
//! variants carried together on one chromosome, and in which populations?
//!
//! Only founders (samples without a parent in the panel) are counted, so
//! that a chromosome transmitted in a trio is not counted twice.

use crate::freq::{self, Key, PopCount, POPS};
use std::collections::{HashMap, HashSet};

/// One alternate allele of the panel with the haplotypes that carry it.
struct Carried {
    /// First and last reference position affected.
    from: u64,
    to: u64,
    /// `None` for symbolic alleles (structural variants).
    key: Option<Key>,
    bits: Vec<u64>,
}

pub struct Panel {
    words: usize,
    /// Founder haplotypes per superpopulation.
    masks: [Vec<u64>; 5],
    pub founders: usize,
    sites: HashMap<String, Vec<Carried>>,
}

pub struct Assessment {
    /// Indices of the allele's variants that are not sites of the panel.
    pub untested: Vec<usize>,
    /// Haplotypes carrying every tested variant, per superpopulation.
    pub carrying: [PopCount; 5],
    /// Haplotypes carrying every variant and no other panel variant within
    /// the allele; `None` when some variant could not be tested.
    pub exact: Option<[PopCount; 5]>,
}

impl Assessment {
    pub fn status(&self, n_variants: usize) -> &'static str {
        let any = |c: &[PopCount; 5]| c.iter().any(|p| p.ac > 0);
        match &self.exact {
            None if self.untested.len() == n_variants => "untestable",
            None => "incomplete",
            Some(exact) if any(exact) => "found_exact",
            Some(_) if any(&self.carrying) => "found_with_extra",
            Some(_) => "not_found",
        }
    }
}

pub fn total(counts: &[PopCount; 5]) -> PopCount {
    counts.iter().fold(PopCount::default(), |t, p| PopCount {
        ac: t.ac + p.ac,
        an: t.an + p.an,
        het: 0,
        hom: 0,
    })
}

impl Panel {
    /// `header`: the `#CHROM` line of the phased VCF. `pedigree`: lines of
    /// `family sample father mother sex population superpopulation`.
    pub fn new(header: &str, pedigree: &str) -> anyhow::Result<Self> {
        let samples: Vec<&str> = header.trim_end().split('\t').skip(9).collect();
        if samples.is_empty() {
            anyhow::bail!("the phased VCF header lists no samples");
        }
        let present: HashSet<&str> = samples.iter().copied().collect();
        let mut group: HashMap<&str, usize> = HashMap::new();
        for line in pedigree.lines().skip(1) {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 7 {
                continue;
            }
            let founder = !present.contains(f[2]) && !present.contains(f[3]);
            if let (true, Some(pop)) = (founder, POPS.iter().position(|p| *p == f[6])) {
                group.insert(f[1], pop);
            }
        }
        let words = (2 * samples.len()).div_ceil(64);
        let mut masks: [Vec<u64>; 5] = std::array::from_fn(|_| vec![0u64; words]);
        let mut founders = 0;
        for (i, sample) in samples.iter().enumerate() {
            if let Some(pop) = group.get(sample) {
                founders += 1;
                for h in [2 * i, 2 * i + 1] {
                    masks[*pop][h / 64] |= 1 << (h % 64);
                }
            }
        }
        if founders == 0 {
            anyhow::bail!("no sample of the phased VCF is a founder with a superpopulation in the pedigree file");
        }
        Ok(Panel {
            words,
            masks,
            founders,
            sites: HashMap::new(),
        })
    }

    /// Add one data line of the phased VCF.
    pub fn add_line(&mut self, contig: &str, line: &str) {
        // The nine fixed columns, then one genotype per sample. A line
        // holds thousands of genotypes, nearly all "0|0": those are
        // recognised without parsing.
        let mut columns = line.splitn(10, '\t');
        let f: Vec<&str> = columns.by_ref().take(9).collect();
        let Some(genotypes) = columns.next() else { return };
        let Ok(pos) = f[1].parse::<u64>() else { return };
        let alts: Vec<&str> = f[4].split(',').collect();
        let mut bits: Vec<Vec<u64>> = vec![vec![0u64; self.words]; alts.len()];
        for (i, field) in genotypes.split('\t').enumerate() {
            if field == "0|0" || 2 * i + 1 >= 64 * self.words {
                continue;
            }
            let gt = field.split(':').next().unwrap_or("");
            for (side, code) in gt.split(['|', '/']).take(2).enumerate() {
                if let Ok(a) = code.parse::<usize>() {
                    if a >= 1 && a <= alts.len() {
                        let h = 2 * i + side;
                        bits[a - 1][h / 64] |= 1 << (h % 64);
                    }
                }
            }
        }
        let out = self.sites.entry(contig.to_string()).or_default();
        for (alt, bits) in alts.iter().zip(bits) {
            if bits.iter().all(|w| *w == 0) {
                continue;
            }
            if alt.starts_with('<') || alt.contains('[') || alt.contains(']') || *alt == "*" {
                let end = f[7]
                    .split(';')
                    .find_map(|x| x.strip_prefix("END="))
                    .and_then(|x| x.parse::<u64>().ok())
                    .unwrap_or(pos);
                out.push(Carried {
                    from: pos,
                    to: end.max(pos),
                    key: None,
                    bits,
                });
            } else {
                let key = freq::trim(pos, f[3].as_bytes(), alt.as_bytes());
                let to = key.0 + key.1.len().max(1) as u64 - 1;
                out.push(Carried {
                    from: key.0,
                    to,
                    key: Some(key),
                    bits,
                });
            }
        }
    }

    fn count(&self, bits: &[u64]) -> [PopCount; 5] {
        std::array::from_fn(|p| PopCount {
            ac: bits.iter().zip(&self.masks[p]).map(|(b, m)| (b & m).count_ones()).sum(),
            an: self.masks[p].iter().map(|m| m.count_ones()).sum(),
            het: 0,
            hom: 0,
        })
    }

    /// `span`: first and last reference position covered by the allele.
    pub fn assess(&self, contig: &str, span: (u64, u64), variants: &[&Key]) -> Assessment {
        let empty = Vec::new();
        let sites = self.sites.get(contig).unwrap_or(&empty);
        let mut all = vec![u64::MAX; self.words];
        let mut untested = Vec::new();
        for (i, key) in variants.iter().enumerate() {
            match sites.iter().find(|s| s.key.as_ref() == Some(*key)) {
                Some(site) => all.iter_mut().zip(&site.bits).for_each(|(a, b)| *a &= b),
                None => untested.push(i),
            }
        }
        let carrying = self.count(&all);
        let exact = untested.is_empty().then(|| {
            let mut clean = all.clone();
            for site in sites.iter().filter(|s| s.from <= span.1 && s.to >= span.0) {
                if site.key.as_ref().is_some_and(|k| variants.contains(&k)) {
                    continue;
                }
                clean.iter_mut().zip(&site.bits).for_each(|(c, b)| *c &= !b);
            }
            self.count(&clean)
        });
        Assessment { untested, carrying, exact }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn panel() -> Panel {
        // s1 and s2 are founders, s3 is their child.
        let header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3";
        let ped = "FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n\
                   f s1 0 0 1 YRI AFR\nf s2 0 0 2 GBR EUR\nf s3 s1 s2 1 YRI AFR\n";
        let mut p = Panel::new(header, ped).unwrap();
        p.add_line("c", "c\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1|0\t0|0\t1|0");
        p.add_line("c", "c\t110\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|1\t1|1");
        p.add_line("c", "c\t120\t.\tG\tA\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0");
        p
    }

    #[test]
    fn combination_is_counted_on_founder_haplotypes() {
        let p = panel();
        assert_eq!(p.founders, 2);
        let (a, b): (Key, Key) = ((100, b"A".to_vec(), b"G".to_vec()), (110, b"C".to_vec(), b"T".to_vec()));
        let r = p.assess("c", (90, 130), &[&a, &b]);
        assert_eq!(r.status(2), "found_exact");
        assert_eq!((total(&r.carrying).ac, total(&r.carrying).an), (1, 4));
        assert_eq!(r.exact.unwrap()[0].ac, 1);
        // 110 alone: three founder haplotypes carry it, one of them with 100, one with 120.
        let r = p.assess("c", (90, 130), &[&b]);
        assert_eq!((total(&r.carrying).ac, total(&r.exact.unwrap()).ac), (3, 1));
        // A variant that is not a site of the panel.
        let missing: Key = (105, b"T".to_vec(), b"C".to_vec());
        let r = p.assess("c", (90, 130), &[&a, &missing]);
        assert_eq!((r.status(2), r.untested.clone()), ("incomplete", vec![1]));
        // No variants: haplotypes identical to the reference over the span.
        let r = p.assess("c", (90, 130), &[]);
        assert_eq!(total(&r.exact.unwrap()).ac, 1);
    }
}

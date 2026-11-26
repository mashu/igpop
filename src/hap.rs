//! Haplotype check against a phased panel: is the allele's combination of
//! variants carried together on one chromosome, and in which populations?
//!
//! Only founders (samples without a parent in the panel) are counted, so
//! that a chromosome transmitted in a trio is not counted twice.

use crate::population::KGP_GROUPS;
use crate::stats::PopCount;
use crate::variant::{trim, Key};
use std::collections::{BTreeMap, HashMap, HashSet};

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
    /// Founder haplotypes per superpopulation, and all of them.
    masks: [Vec<u64>; 5],
    founder_mask: Vec<u64>,
    pub founders: usize,
    sites: HashMap<String, Vec<Carried>>,
}

#[derive(Clone)]
pub struct Assessment {
    /// Indices of the allele's variants that are not sites of the panel.
    pub untested: Vec<usize>,
    /// Haplotypes carrying every tested variant, per superpopulation.
    pub carrying: [PopCount; 5],
    /// Haplotypes carrying every variant and no other panel variant within
    /// the allele; `None` when some variant could not be tested.
    pub exact: Option<[PopCount; 5]>,
    /// Founders (people) carrying such a chromosome.
    pub exact_donors: Option<u32>,
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
            if let (true, Some(pop)) = (founder, KGP_GROUPS.iter().position(|g| g.code == f[6])) {
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
        let founder_mask = (0..words).map(|w| masks.iter().fold(0, |all, m| all | m[w])).collect();
        Ok(Panel {
            words,
            masks,
            founder_mask,
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
                let key = trim(pos, f[3].as_bytes(), alt.as_bytes());
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
        std::array::from_fn(|p| {
            PopCount::new(
                bits.iter().zip(&self.masks[p]).map(|(b, m)| (b & m).count_ones()).sum(),
                self.masks[p].iter().map(|m| m.count_ones()).sum(),
            )
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
        let clean = untested.is_empty().then(|| {
            let mut clean = all.clone();
            for site in sites.iter().filter(|s| s.from <= span.1 && s.to >= span.0) {
                if site.key.as_ref().is_some_and(|k| variants.contains(&k)) {
                    continue;
                }
                clean.iter_mut().zip(&site.bits).for_each(|(c, b)| *c &= !b);
            }
            clean
        });
        // A person's two chromosomes are neighbouring bits.
        let donors = |bits: &[u64]| -> u32 {
            bits.iter()
                .zip(&self.founder_mask)
                .map(|(b, m)| (((b & m) | ((b & m) >> 1)) & 0x5555_5555_5555_5555).count_ones())
                .sum()
        };
        Assessment {
            untested,
            carrying,
            exact: clean.as_deref().map(|bits| self.count(bits)),
            exact_donors: clean.as_deref().map(donors),
        }
    }
}

/// How often an arbitrary single-base change of a real chromosome gives a
/// sequence that itself exists in the panel, and on how many chromosomes.
///
/// For every founder chromosome and every possible single-base change
/// within a locus, the changed sequence is looked up among the founder
/// chromosomes. This is what a sequencing or PCR error on a real allele
/// would look like in the haplotype check, and so the yardstick for how
/// much a "found" is worth.
#[derive(Default)]
pub struct Coincidence {
    /// Changed sequences found on exactly `k` chromosomes: `k -> number`.
    found: BTreeMap<u32, u64>,
    /// All (chromosome, change) pairs tried.
    tried: u64,
}

impl Coincidence {
    /// Add one locus; `span` = first and last reference position.
    pub fn add_locus(&mut self, panel: &Panel, contig: &str, span: (u64, u64)) {
        let empty = Vec::new();
        let sites: Vec<&Carried> = panel
            .sites
            .get(contig)
            .unwrap_or(&empty)
            .iter()
            .filter(|s| s.from <= span.1 && s.to >= span.0)
            .collect();
        // The sequence of each founder chromosome, as the sites it carries.
        let mut carried: HashMap<usize, Vec<u16>> = HashMap::new();
        for (local, site) in sites.iter().enumerate() {
            for (w, word) in site.bits.iter().zip(&panel.founder_mask).map(|(b, m)| b & m).enumerate() {
                let mut rest = word;
                while rest != 0 {
                    carried
                        .entry(64 * w + rest.trailing_zeros() as usize)
                        .or_default()
                        .push(local as u16);
                    rest &= rest - 1;
                }
            }
        }
        let chromosomes = 2 * panel.founders as u32;
        let mut sequences: HashMap<Vec<u16>, u32> = HashMap::new();
        sequences.insert(Vec::new(), chromosomes - carried.len() as u32);
        for signature in carried.into_values() {
            *sequences.entry(signature).or_default() += 1;
        }
        sequences.retain(|_, count| *count > 0);
        // Single-base sites inside the locus, with their position.
        let single: Vec<(u16, u64)> = sites
            .iter()
            .enumerate()
            .filter_map(|(local, site)| match &site.key {
                Some((pos, r, a)) if r.len() == 1 && a.len() == 1 && *pos >= span.0 && *pos <= span.1 => Some((local as u16, *pos)),
                _ => None,
            })
            .collect();
        for (sequence, weight) in &sequences {
            for (site, pos) in &single {
                // The base at this position changes: to the reference if
                // the chromosome carries this allele, otherwise to this
                // allele (replacing another one it may carry there).
                let mut changed: Vec<u16> = sequence
                    .iter()
                    .copied()
                    .filter(|s| !single.iter().any(|(other, at)| other == s && at == pos))
                    .collect();
                if !sequence.contains(site) {
                    changed.push(*site);
                    changed.sort_unstable();
                }
                if let Some(count) = sequences.get(&changed) {
                    *self.found.entry(*count).or_default() += *weight as u64;
                }
            }
        }
        self.tried += chromosomes as u64 * 3 * (span.1 - span.0 + 1);
    }

    /// Share of arbitrary changes found on at least `count` chromosomes.
    pub fn rate(&self, count: u32) -> f64 {
        if self.tried == 0 {
            return f64::NAN;
        }
        self.found.range(count..).map(|(_, n)| *n).sum::<u64>() as f64 / self.tried as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::total;

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
        assert_eq!((total(&r.exact.unwrap()).ac, r.exact_donors), (1, Some(1)));
    }

    #[test]
    fn coincidence_counts_changed_chromosomes_that_exist() {
        // Founder chromosomes: {100,110}, {110,120}, {}, {110}.
        let p = panel();
        let mut c = Coincidence::default();
        c.add_locus(&p, "c", (91, 130));
        // 4 chromosomes x 3 changes x 40 positions were tried.
        assert_eq!(c.tried, 480);
        // One base apart: {100,110}-{110}, {110,120}-{110}, {}-{110}; each
        // pair counts once from either side, every sequence exists once.
        assert_eq!(c.found.get(&1), Some(&6));
        assert!((c.rate(1) - 6.0 / 480.0).abs() < 1e-12);
        assert_eq!(c.rate(2), 0.0);
    }
}

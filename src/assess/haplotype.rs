//! The haplotype check of every allele, and how often an arbitrary change
//! would pass it.

use crate::evidence::Evidence;
use crate::hap;
use crate::locus::{Placed, Status, Window};
use crate::stats;

/// The haplotype check supports an allele.
#[derive(Clone, Copy)]
pub struct HaplotypeSupport {
    /// Founder chromosomes with exactly the allele's sequence.
    pub chromosomes: u32,
    /// The people carrying them.
    pub donors: u32,
    /// Share of arbitrary single-base changes of real chromosomes found on
    /// at least that many chromosomes.
    pub coincidence: f64,
}

/// The haplotype check of one allele.
#[derive(Clone)]
pub struct AlleleHaplotype {
    pub assessment: Option<hap::Assessment>,
    pub coincidence: Option<f64>,
    pub support: Option<HaplotypeSupport>,
}

/// The haplotype check of every allele, and its calibration.
pub struct HaplotypeCheck {
    /// Per allele, in the order of `placed`.
    pub alleles: Vec<Option<hap::Assessment>>,
    /// `None` without a phased panel.
    pub coincidence: Option<hap::Coincidence>,
}

impl HaplotypeCheck {
    pub fn new(ev: &Evidence, windows: &[Window], placed: &[Placed]) -> Self {
        let Some(panel) = &ev.panel else {
            return HaplotypeCheck {
                alleles: placed.iter().map(|_| None).collect(),
                coincidence: None,
            };
        };
        let alleles = placed
            .iter()
            .map(|p| match &p.placement {
                Some(at) if p.status != Status::LowIdentity => {
                    let keys: Vec<_> = p.variants.iter().map(|c| &c.key).collect();
                    Some(panel.assess(&at.contig, at.span, &keys))
                }
                _ => None,
            })
            .collect();
        // Calibrated on the loci of the run, each counted once.
        let mut loci: Vec<(&str, u64, u64)> = placed
            .iter()
            .filter_map(|p| p.placement.as_ref())
            .map(|at| &windows[at.window].coord)
            .map(|c| (c.contig.as_str(), c.start, c.end))
            .collect();
        loci.sort();
        loci.dedup();
        let mut coincidence = hap::Coincidence::default();
        for (contig, start, end) in loci {
            coincidence.add_locus(panel, contig, (start, end));
        }
        HaplotypeCheck {
            alleles,
            coincidence: Some(coincidence),
        }
    }

    /// The check of the allele at `index`. Its haplotype supports it if
    /// the exact sequence is on at least two founder chromosomes (one can
    /// be a genotyping error) and fewer than `limit` of arbitrary changes
    /// are found that often.
    pub fn of(&self, index: usize, n_variants: usize, limit: f64) -> AlleleHaplotype {
        let assessment = self.alleles[index].clone();
        let exact = assessment.as_ref().and_then(|a| a.exact.as_ref()).map(|e| stats::total(e).ac);
        let coincidence = match (exact, &self.coincidence) {
            (Some(chromosomes), Some(c)) if n_variants > 0 => Some(c.rate(chromosomes)),
            _ => None,
        };
        let support = match (exact, coincidence, &assessment) {
            (Some(chromosomes), Some(rate), Some(a)) if chromosomes >= 2 && rate <= limit => Some(HaplotypeSupport {
                chromosomes,
                donors: a.exact_donors.unwrap_or(0),
                coincidence: rate,
            }),
            _ => None,
        };
        AlleleHaplotype {
            assessment,
            coincidence,
            support,
        }
    }

    /// The smallest number of chromosomes that counts as support.
    pub fn required(&self, limit: f64) -> Option<u32> {
        let coincidence = self.coincidence.as_ref()?;
        (2..=100_000).find(|c| coincidence.rate(*c) <= limit)
    }
}

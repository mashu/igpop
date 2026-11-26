//! Allele counts and what they say about a frequency.
//!
//! Every statement igpop makes about a frequency starts from two numbers:
//! `ac` copies of an allele among `an` chromosomes examined. Bounds are
//! exact (Clopper-Pearson) one-sided 95% bounds, the test is the exact
//! binomial test; neither assumes a large sample.

use statrs::distribution::{Beta, ChiSquared, ContinuousCDF};

/// Copies of an allele (`ac`) among the chromosomes examined (`an`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PopCount {
    pub ac: u32,
    pub an: u32,
}

impl PopCount {
    pub fn new(ac: u32, an: u32) -> Self {
        PopCount { ac, an }
    }

    pub fn af(&self) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else {
            self.ac as f64 / self.an as f64
        }
    }

    /// Exact one-sided 95% lower bound of the frequency.
    pub fn lower95(&self) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else if self.ac == 0 {
            0.0
        } else {
            beta_quantile(self.ac as f64, (self.an - self.ac.min(self.an)) as f64 + 1.0, 0.05)
        }
    }

    /// Exact one-sided 95% upper bound of the frequency.
    pub fn upper95(&self) -> f64 {
        if self.an == 0 {
            f64::NAN
        } else if self.ac >= self.an {
            1.0
        } else {
            beta_quantile(self.ac as f64 + 1.0, (self.an - self.ac) as f64, 0.95)
        }
    }

    /// Exact binomial p-value against "the frequency is below
    /// `threshold`": the probability of at least `ac` copies in `an`
    /// chromosomes if the frequency were exactly `threshold`.
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

/// The sum of several counts.
pub fn total(counts: &[PopCount]) -> PopCount {
    counts
        .iter()
        .fold(PopCount::default(), |t, c| PopCount::new(t.ac + c.ac, t.an + c.an))
}

/// Index of the count with the highest lower bound.
pub fn strongest(counts: &[PopCount]) -> usize {
    (0..counts.len())
        .max_by(|a, b| {
            let (x, y) = (counts[*a].lower95(), counts[*b].lower95());
            x.partial_cmp(&y).unwrap_or(if x.is_nan() {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            })
        })
        .unwrap_or(0)
}

/// True if `an` chromosomes are enough to tell a frequency of `limit`
/// from a rare variant: not seeing the allele at all would bound its
/// frequency below `limit`.
pub fn informative(an: u32, limit: f64) -> bool {
    PopCount::new(0, an).upper95() < limit
}

/// Likelihood-ratio (G) test of equal frequency in all groups; NaN when
/// the allele is absent or fixed. Descriptive only: no verdict uses it.
pub fn heterogeneity_p(groups: &[PopCount]) -> f64 {
    let all = total(groups);
    if all.ac == 0 || all.ac == all.an {
        return f64::NAN;
    }
    let p = all.af();
    let term = |o: f64, e: f64| if o > 0.0 { o * (o / e).ln() } else { 0.0 };
    let used: Vec<&PopCount> = groups.iter().filter(|c| c.an > 0).collect();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_test_adapts_to_sample_size() {
        let c = PopCount::new;
        // 1.5% observed: three copies in 200 is not evidence for >= 1%, 300 in 20,000 is.
        assert!(c(3, 200).tail_p(0.01) > 0.05);
        assert!(c(300, 20_000).tail_p(0.01) < 1e-6);
        assert!(c(3, 200).lower95() < 0.01 && c(300, 20_000).lower95() > 0.01);
        assert_eq!(c(0, 500).lower95(), 0.0);
        assert!((c(0, 5008).upper95() - 0.000598).abs() < 1e-5);
    }

    /// The test and the bound are the same statement: the p-value against
    /// "below t" is at most 5% exactly when the 95% lower bound reaches t.
    #[test]
    fn the_test_and_the_bound_agree() {
        for (ac, an) in [(2, 126), (5, 400), (30, 2000), (300, 20_000)] {
            let c = PopCount::new(ac, an);
            let bound = c.lower95();
            assert!(c.tail_p(bound * 0.999) <= 0.05);
            assert!(c.tail_p(bound * 1.001) > 0.05);
        }
    }

    #[test]
    fn a_stratum_must_be_large_enough_to_tell() {
        assert!(!informative(200, 0.01));
        assert!(informative(400, 0.01));
    }
}

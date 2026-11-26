//! The main score: how much more typical a variant's population support
//! is of the variants of real alleles than of arbitrary base changes.
//!
//! Support is one number per variant: the largest frequency its counts
//! guarantee in some population (see `evidence::support`). Two groups of
//! single-base changes are compared in the same data: every change that
//! is possible in the loci of the run ("arbitrary"), and the changes found
//! in alleles accepted as real ("real").
//!
//! A change that no dataset has seen is a category of its own: its score
//! is the share of real variants never seen over the share of arbitrary
//! changes never seen.
//!
//! For a change that was seen, the score is a smooth curve of support. A
//! logistic regression on x = log10(support) estimates, at every x, the
//! odds that a seen change is one of the real variants rather than one of
//! the arbitrary changes; multiplied by the ratio of the group sizes this
//! is the ratio of the two densities at x, the likelihood ratio:
//!
//! ```text
//! log score(x) = b0 + b1 x + b2 x^2
//! ```
//!
//! More support must never mean a lower score. The curved term is used
//! only if the fit is better with it (Akaike's criterion) and the curve
//! still rises over the whole range; otherwise the line is used. A seen
//! change never scores below one that was never seen.

use anyhow::{bail, Result};

/// Upper edges of the support levels; level 0 is "no support at all". The
/// levels do not enter the score: they are a ruler for tables.
const EDGES: [f64; 5] = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
pub const BINS: usize = EDGES.len() + 2;
pub const LABELS: [&str; BINS] = ["none", "<0.001%", "0.001-0.01%", "0.01-0.1%", "0.1-1%", "1-10%", ">=10%"];

pub fn bin(support: f64) -> usize {
    if support.is_nan() || support <= 0.0 {
        0
    } else {
        1 + EDGES.iter().take_while(|edge| support >= **edge).count()
    }
}

/// The range of support of level `i` (1 or more), for display.
pub fn level_range(i: usize) -> (f64, f64) {
    (
        if i == 1 { FLOOR } else { EDGES[i - 2] },
        if i == BINS - 1 { 1.0 } else { EDGES[i - 1] },
    )
}

/// Fewer real variants than this cannot place a curve.
pub const MIN_TRUSTED: u64 = 30;

/// The smallest support the curve is evaluated at: one in a million.
const FLOOR: f64 = 1e-6;
/// Supports at which the share of arbitrary changes at least that well
/// supported is kept: eight per step of ten, from FLOOR to 1.
const GRID: usize = 49;
fn grid(k: usize) -> f64 {
    10f64.powf(-6.0 + k as f64 / 8.0)
}

/// The sizes of the two groups a calibration was made from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Counts {
    pub arbitrary: u64,
    pub arbitrary_seen: u64,
    pub trusted: u64,
    pub trusted_seen: u64,
}

#[derive(Clone, Debug)]
pub struct Calibration {
    /// Score of a change no dataset has seen.
    pub never: f64,
    /// log score = b0 + b1 x + b2 x^2 for a change seen somewhere.
    pub beta: [f64; 3],
    /// Covariance of `beta`, for the interval of a score.
    pub cov: [[f64; 3]; 3],
    /// Share of all arbitrary changes whose support is at least `grid(k)`.
    pub tail: [f64; GRID],
    pub counts: Counts,
    /// The support of each real variant that was seen, sorted: the points
    /// the curve was drawn through, kept so that they can be shown.
    pub real: Vec<f64>,
}

impl Calibration {
    fn curve(&self, support: f64) -> (f64, f64) {
        let x = support.clamp(FLOOR, 1.0).log10();
        let v = [1.0, x, x * x];
        let log = v.iter().zip(&self.beta).map(|(a, b)| a * b).sum::<f64>();
        let variance: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| (i, j)))
            .map(|(i, j)| v[i] * self.cov[i][j] * v[j])
            .sum();
        (log, variance.max(0.0).sqrt())
    }

    /// The score of a variant with this support.
    pub fn score(&self, support: f64) -> f64 {
        if support.is_nan() || support <= 0.0 {
            self.never
        } else {
            self.curve(support).0.exp().max(self.never)
        }
    }

    /// A 95% interval for the score: how well the calibration knows it.
    pub fn interval(&self, support: f64) -> (f64, f64) {
        if support.is_nan() || support <= 0.0 {
            let c = self.counts;
            let (t0, a0) = (
                (c.trusted - c.trusted_seen) as f64 + 0.5,
                (c.arbitrary - c.arbitrary_seen) as f64 + 0.5,
            );
            let se = (1.0 / t0 - 1.0 / (c.trusted as f64 + 1.0) + 1.0 / a0 - 1.0 / (c.arbitrary as f64 + 1.0))
                .max(0.0)
                .sqrt();
            return (self.never * (-1.96 * se).exp(), self.never * (1.96 * se).exp());
        }
        let (log, se) = self.curve(support);
        (
            (log - 1.96 * se).exp().max(self.never.min(self.score(support))),
            (log + 1.96 * se).exp().max(self.never),
        )
    }

    /// The share of all arbitrary changes that are supported at least as
    /// well as a variant with this support: how often chance does as well.
    pub fn chance(&self, support: f64) -> f64 {
        if support.is_nan() || support <= 0.0 {
            return 1.0;
        }
        // The nearest kept support at or below, so that the share is not
        // understated.
        let k = (((support.clamp(FLOOR, 1.0).log10() + 6.0) * 8.0) + 1e-9).floor() as usize;
        self.tail[k.min(GRID - 1)]
    }

    /// The smallest support whose score reaches `score`, if a seen change
    /// can score that much: where a verdict limit lies on the support axis.
    pub fn support_at(&self, score: f64) -> Option<f64> {
        if self.score(1.0) < score {
            return None;
        }
        let (mut low, mut high) = (FLOOR.log10(), 0.0);
        if self.score(FLOOR) >= score {
            return Some(FLOOR);
        }
        for _ in 0..60 {
            let middle = (low + high) / 2.0;
            if self.score(10f64.powf(middle)) >= score {
                high = middle;
            } else {
                low = middle;
            }
        }
        Some(10f64.powf(high))
    }

    /// The two groups counted per level of support, as (possible changes,
    /// real variants): the curve's data on a coarse ruler, for tables.
    pub fn levels(&self) -> [(u64, u64); BINS] {
        let c = self.counts;
        // Possible changes with a support of at least each edge. The edges
        // are points of the kept grid, eight steps apart.
        let at_least = |edge: usize| (self.tail[8 * (edge + 1)] * c.arbitrary as f64).round() as u64;
        let mut out = [(0, 0); BINS];
        out[0] = (c.arbitrary - c.arbitrary_seen, c.trusted - c.trusted_seen);
        out[1].0 = c.arbitrary_seen - at_least(0);
        for (edge, level) in out.iter_mut().skip(2).enumerate() {
            level.0 = at_least(edge) - if edge + 1 < EDGES.len() { at_least(edge + 1) } else { 0 };
        }
        for support in &self.real {
            out[bin(*support)].1 += 1;
        }
        out
    }

    /// The text of `score_model.tsv`, which `parse` reads back.
    pub fn format(&self) -> String {
        let c = &self.cov;
        let mut out = format!(
            "key\tvalue\tshare\nmodel\tlog(score) = b0 + b1*x + b2*x^2 with x = log10(support); never = score of a change never seen\n\
             never\t{:e}\nb0\t{:e}\nb1\t{:e}\nb2\t{:e}\ncov\t{:e},{:e},{:e},{:e},{:e},{:e}\n\
             arbitrary\t{}\narbitrary_seen\t{}\ntrusted\t{}\ntrusted_seen\t{}\n",
            self.never,
            self.beta[0],
            self.beta[1],
            self.beta[2],
            c[0][0],
            c[0][1],
            c[0][2],
            c[1][1],
            c[1][2],
            c[2][2],
            self.counts.arbitrary,
            self.counts.arbitrary_seen,
            self.counts.trusted,
            self.counts.trusted_seen
        );
        for (k, share) in self.tail.iter().enumerate() {
            out.push_str(&format!("tail\t{:e}\t{share:e}\n", grid(k)));
        }
        for support in &self.real {
            out.push_str(&format!("real\t{support:e}\n"));
        }
        out
    }
}

/// Solve the 3x3 system `a x = b` (Gaussian elimination with pivoting).
fn solve(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for col in 0..3 {
        let pivot = (col..3).max_by(|i, j| a[*i][col].abs().total_cmp(&a[*j][col].abs()))?;
        if a[pivot][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..3 {
            let f = a[row][col] / a[col][col];
            let pivot_row = a[col];
            for (cell, above) in a[row].iter_mut().zip(pivot_row).skip(col) {
                *cell -= f * above;
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = [0.0; 3];
    for row in (0..3).rev() {
        let rest: f64 = (row + 1..3).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - rest) / a[row][row];
    }
    Some(x)
}

/// A logistic regression of `y` (1 = real) on 1, x and, with `curved`,
/// x^2, by Newton's method: the coefficients, their covariance and the
/// negative log-likelihood.
fn logistic(xs: &[f64], ys: &[bool], curved: bool) -> Option<([f64; 3], [[f64; 3]; 3], f64)> {
    let terms = if curved { 3 } else { 2 };
    let mut beta = [0.0; 3];
    let likelihood = |beta: &[f64; 3]| -> f64 {
        xs.iter()
            .zip(ys)
            .map(|(x, y)| {
                let z = beta[0] + beta[1] * x + beta[2] * x * x;
                // log(1 + e^z) - y z, written to stay finite for large |z|.
                z.max(0.0) + (-z.abs()).exp().ln_1p() - if *y { z } else { 0.0 }
            })
            .sum()
    };
    let mut nll = likelihood(&beta);
    let mut hessian = [[0.0; 3]; 3];
    for _ in 0..100 {
        let mut gradient = [0.0; 3];
        hessian = [[0.0; 3]; 3];
        for (x, y) in xs.iter().zip(ys) {
            let v = [1.0, *x, x * x];
            let p = 1.0 / (1.0 + (-(beta[0] + beta[1] * x + beta[2] * x * x)).exp());
            let (residual, weight) = (if *y { 1.0 } else { 0.0 } - p, p * (1.0 - p));
            for i in 0..terms {
                gradient[i] += residual * v[i];
                for j in 0..terms {
                    hessian[i][j] += weight * v[i] * v[j];
                }
            }
        }
        if !curved {
            hessian[2][2] = 1.0;
        }
        let step = solve(hessian, gradient)?;
        // Halve the step until the fit does not get worse.
        let mut scale = 1.0;
        let (next, next_nll) = loop {
            let mut next = beta;
            for i in 0..terms {
                next[i] += scale * step[i];
            }
            let value = likelihood(&next);
            if value <= nll + 1e-12 || scale < 1e-6 {
                break (next, value);
            }
            scale /= 2.0;
        };
        let moved = (nll - next_nll).abs();
        (beta, nll) = (next, next_nll);
        if moved < 1e-10 {
            break;
        }
    }
    // The covariance is the inverse of the Hessian at the solution.
    let mut cov = [[0.0; 3]; 3];
    for col in 0..terms {
        let mut unit = [0.0; 3];
        unit[col] = 1.0;
        let column = solve(hessian, unit)?;
        for row in 0..terms {
            cov[row][col] = column[row];
        }
    }
    Some((beta, cov, nll))
}

/// Calibrate the score from the support of every arbitrary change and of
/// every real variant. `None` if the real variants are too few or do not
/// show more support than arbitrary changes.
pub fn calibrate(arbitrary: &[f64], trusted: &[f64]) -> Option<Calibration> {
    let seen = |values: &[f64]| -> Vec<f64> { values.iter().filter(|s| **s > 0.0).map(|s| s.clamp(FLOOR, 1.0).log10()).collect() };
    let (a, t) = (seen(arbitrary), seen(trusted));
    let counts = Counts {
        arbitrary: arbitrary.len() as u64,
        arbitrary_seen: a.len() as u64,
        trusted: trusted.len() as u64,
        trusted_seen: t.len() as u64,
    };
    if counts.trusted < MIN_TRUSTED || a.is_empty() || t.is_empty() {
        return None;
    }
    // Never seen: the ratio of the two shares, with half a count added so
    // that an empty group gives neither zero nor infinity.
    let never = ((counts.trusted - counts.trusted_seen) as f64 + 0.5)
        / (counts.trusted as f64 + 1.0)
        / (((counts.arbitrary - counts.arbitrary_seen) as f64 + 0.5) / (counts.arbitrary as f64 + 1.0));

    let xs: Vec<f64> = a.iter().chain(&t).copied().collect();
    let ys: Vec<bool> = a.iter().map(|_| false).chain(t.iter().map(|_| true)).collect();
    let lowest = xs.iter().copied().fold(0.0, f64::min);
    let rises = |beta: &[f64; 3]| beta[1] + 2.0 * beta[2] * lowest >= 0.0 && beta[1] >= 0.0;
    let line = logistic(&xs, &ys, false)?;
    let curve = logistic(&xs, &ys, true).filter(|(beta, _, nll)| rises(beta) && 2.0 * (line.2 - nll) > 2.0);
    let (mut beta, cov, _) = curve.unwrap_or(line);
    if !rises(&beta) || beta[1] <= 0.0 {
        return None;
    }
    // From the odds among seen changes to the ratio of the two densities:
    // times the ratio of the group sizes.
    beta[0] += (counts.arbitrary as f64 / counts.trusted as f64).ln();

    let mut tail = [0.0; GRID];
    for (k, share) in tail.iter_mut().enumerate() {
        let at = grid(k).log10() - 1e-9;
        *share = a.iter().filter(|x| **x >= at).count() as f64 / counts.arbitrary as f64;
    }
    let mut real: Vec<f64> = trusted.iter().copied().filter(|s| *s > 0.0).collect();
    real.sort_by(f64::total_cmp);
    Some(Calibration {
        never,
        beta,
        cov,
        tail,
        counts,
        real,
    })
}

/// The calibration of a run of 566 human IGHV alleles (KIARVA) on GRCh38
/// against every default dataset of that build: used when a run has too
/// few real variants of its own.
const BUILTIN: &str = include_str!("../data/score_model.tsv");

pub fn builtin() -> Calibration {
    parse(BUILTIN, "built-in score model").expect("the built-in score model is complete")
}

/// A calibration from the text of a `<prefix>.score_model.tsv` written by
/// an earlier run; `label` names it in error messages.
pub fn parse(text: &str, label: &str) -> Result<Calibration> {
    let value = |key: &str| -> Result<&str> {
        text.lines()
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .find(|f| f.first() == Some(&key) && f.len() >= 2)
            .map(|f| f[1])
            .ok_or_else(|| anyhow::anyhow!("{label} is not a score model of igpop: no line '{key}' (expected a score_model.tsv)"))
    };
    let number = |key: &str| -> Result<f64> { Ok(value(key)?.trim().parse::<f64>()?) };
    let never = number("never")?;
    let beta = [number("b0")?, number("b1")?, number("b2")?];
    let counts = Counts {
        arbitrary: number("arbitrary")? as u64,
        arbitrary_seen: number("arbitrary_seen")? as u64,
        trusted: number("trusted")? as u64,
        trusted_seen: number("trusted_seen")? as u64,
    };
    let c: Vec<f64> = value("cov")?.split(',').filter_map(|x| x.trim().parse().ok()).collect();
    if c.len() != 6 {
        bail!("{label}: the line 'cov' needs six numbers");
    }
    let cov = [[c[0], c[1], c[2]], [c[1], c[3], c[4]], [c[2], c[4], c[5]]];
    let shares: Vec<f64> = text
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|f| f.first() == Some(&"tail") && f.len() >= 3)
        .filter_map(|f| f[2].trim().parse().ok())
        .collect();
    if shares.len() != GRID {
        bail!("{label}: expected {GRID} lines 'tail', found {}", shares.len());
    }
    let mut tail = [0.0; GRID];
    tail.copy_from_slice(&shares);
    let real = text
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|f| f.first() == Some(&"real") && f.len() >= 2)
        .filter_map(|f| f[1].trim().parse().ok())
        .collect();
    Ok(Calibration {
        never,
        beta,
        cov,
        tail,
        counts,
        real,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_follow_decades() {
        assert_eq!(bin(0.0), 0);
        assert_eq!(bin(5e-6), 1);
        assert_eq!(bin(1e-5), 2);
        assert_eq!(bin(0.005), 4);
        assert_eq!(bin(0.02), 5);
        assert_eq!(bin(0.5), 6);
        assert_eq!(level_range(4), (1e-3, 1e-2));
    }

    /// Arbitrary changes mostly unseen or rare, real variants mostly frequent.
    fn example() -> (Vec<f64>, Vec<f64>) {
        let mut arbitrary = vec![0.0; 5000];
        let mut trusted = Vec::new();
        for i in 0..5000 {
            // Log-uniform between one in a million and one in ten.
            arbitrary.push(10f64.powf(-6.0 + 5.0 * (i as f64 + 0.5) / 5000.0));
        }
        for i in 0..200 {
            // Log-uniform between one in a thousand and one.
            trusted.push(10f64.powf(-3.0 + 3.0 * (i as f64 + 0.5) / 200.0));
        }
        (arbitrary, trusted)
    }

    /// The levels are recovered from what a calibration keeps: they
    /// equal a direct count of the two groups.
    #[test]
    fn levels_equal_a_direct_count() {
        let (arbitrary, trusted) = example();
        let levels = calibrate(&arbitrary, &trusted).unwrap().levels();
        for level in 0..BINS {
            let count = |values: &[f64]| values.iter().filter(|s| bin(**s) == level).count() as u64;
            assert_eq!(levels[level], (count(&arbitrary), count(&trusted)), "level {level}");
        }
    }

    #[test]
    fn score_rises_smoothly_with_support() {
        let (arbitrary, trusted) = example();
        let c = calibrate(&arbitrary, &trusted).unwrap();
        let supports = [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.5];
        let scores: Vec<f64> = supports.iter().map(|s| c.score(*s)).collect();
        assert!(scores.windows(2).all(|w| w[1] >= w[0]), "{scores:?}");
        // No jump at the edge of a level.
        assert!(c.score(0.0101) / c.score(0.0099) < 1.1);
        assert!(c.never < 0.1 && c.score(0.5) > 10.0);
        // The interval holds the score, and chance falls with support.
        let (low, high) = c.interval(0.01);
        assert!(low <= c.score(0.01) && c.score(0.01) <= high);
        assert_eq!(c.chance(0.0), 1.0);
        assert!(c.chance(0.05) < c.chance(0.001) && c.chance(0.001) < 0.5);
    }

    #[test]
    fn too_few_or_uninformative_variants_give_no_calibration() {
        let (arbitrary, trusted) = example();
        assert!(calibrate(&arbitrary, &trusted[..10]).is_none());
        // "Real" variants that are rarer than arbitrary changes: no rising curve.
        let rare: Vec<f64> = (0..100).map(|i| 10f64.powf(-6.0 + 2.0 * (i as f64 + 0.5) / 100.0)).collect();
        assert!(calibrate(&arbitrary, &rare).is_none());
    }

    #[test]
    fn a_model_survives_writing_and_reading() {
        let (arbitrary, trusted) = example();
        let c = calibrate(&arbitrary, &trusted).unwrap();
        let again = parse(&c.format(), "test").unwrap();
        for s in [0.0, 2e-5, 0.004, 0.3] {
            assert!((c.score(s) - again.score(s)).abs() < 1e-9 * c.score(s));
            assert!((c.chance(s) - again.chance(s)).abs() < 1e-12);
        }
        assert_eq!(c.counts, again.counts);
    }
}

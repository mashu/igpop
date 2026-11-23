//! The main score: how much more often a variant's level of population
//! support occurs among trusted variants than among arbitrary base changes.
//!
//! Support is one number per variant (the largest lower confidence bound of
//! its frequency over all populations of all count-based datasets). Its
//! distribution is tabulated twice, over every possible single-base change
//! in the loci of the run and over the changes found in trusted alleles;
//! the score is the ratio of the two proportions in the variant's bin.

/// Upper edges of the support bins; bin 0 is "no support at all".
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

pub struct Calibration {
    /// Score of each bin after pooling.
    pub score: [f64; BINS],
    /// Bins pooled together share a group number.
    pub group: [usize; BINS],
}

/// Fewer trusted variants than this cannot fill the bins.
pub const MIN_TRUSTED: u64 = 30;

/// Ratio of proportions per bin, with half a count added to every bin so
/// that an empty bin does not give zero or infinity. More support can never
/// mean a lower score, so neighbouring bins that break that order are pooled
/// (isotonic regression, pool-adjacent-violators).
pub fn calibrate(arbitrary: [u64; BINS], trusted: [u64; BINS]) -> Option<Calibration> {
    let (n_arbitrary, n_trusted): (u64, u64) = (arbitrary.iter().sum(), trusted.iter().sum());
    if n_trusted < MIN_TRUSTED || n_arbitrary == 0 {
        return None;
    }
    let total_a = n_arbitrary as f64 + 0.5 * BINS as f64;
    let total_t = n_trusted as f64 + 0.5 * BINS as f64;
    // Blocks of pooled bins: (first bin, last bin, arbitrary count, trusted count).
    let ratio = |a: f64, t: f64| (t / total_t) / (a / total_a);
    let mut blocks: Vec<(usize, usize, f64, f64)> = Vec::new();
    for i in 0..BINS {
        blocks.push((i, i, arbitrary[i] as f64 + 0.5, trusted[i] as f64 + 0.5));
        while blocks.len() > 1 {
            let (last, before) = (blocks[blocks.len() - 1], blocks[blocks.len() - 2]);
            if ratio(last.2, last.3) >= ratio(before.2, before.3) {
                break;
            }
            blocks.truncate(blocks.len() - 2);
            blocks.push((before.0, last.1, before.2 + last.2, before.3 + last.3));
        }
    }
    let (mut score, mut group) = ([0.0; BINS], [0; BINS]);
    for (g, (first, last, a, t)) in blocks.iter().enumerate() {
        for i in *first..=*last {
            score[i] = ratio(*a, *t);
            group[i] = g + 1;
        }
    }
    Some(Calibration { score, group })
}

/// Counts behind the built-in calibration: 559 human IGHV alleles (171
/// single-base variants of the trusted ones, 38,316 possible changes) on
/// GRCh38 against 1000 Genomes 30x, gnomAD v4.1 genomes and the datasets
/// relayed by dbSNP, --common 0.01.
pub const BUILTIN_ARBITRARY: [u64; BINS] = [18812, 5488, 7809, 3899, 1614, 332, 362];
pub const BUILTIN_TRUSTED: [u64; BINS] = [0, 0, 0, 1, 8, 26, 136];

/// The calibration used when a run has too few trusted variants of its own.
pub fn builtin() -> Calibration {
    calibrate(BUILTIN_ARBITRARY, BUILTIN_TRUSTED).expect("built-in counts suffice")
}

/// Probability that a variant is real, given the share of real variants
/// expected before looking at the data (`prior`) and its score.
pub fn posterior(prior: f64, score: f64) -> f64 {
    prior * score / (prior * score + 1.0 - prior)
}

/// Scores from the text of a `<prefix>.score_calibration.tsv` written by an
/// earlier run; `label` names it in error messages.
pub fn parse(text: &str, label: &str) -> anyhow::Result<Calibration> {
    let (mut score, mut group) = ([f64::NAN; BINS], [0; BINS]);
    for line in text.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if let (Some(i), Some(value)) = (
            LABELS.iter().position(|l| Some(l) == f.first()),
            f.get(5).and_then(|v| v.parse::<f64>().ok()),
        ) {
            score[i] = value;
            group[i] = f.get(6).and_then(|g| g.parse().ok()).unwrap_or(0);
        }
    }
    if score.iter().any(|x| x.is_nan()) {
        anyhow::bail!("{label} does not hold a score for every support bin");
    }
    Ok(Calibration { score, group })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_follow_decades() {
        assert_eq!(bin(0.0), 0);
        assert_eq!(bin(5e-6), 1);
        assert_eq!(bin(1e-5), 2);
        assert_eq!(bin(0.005), 4);
        assert_eq!(bin(0.02), 5);
        assert_eq!(bin(0.5), 6);
    }

    #[test]
    fn builtin_and_posterior() {
        let c = builtin();
        assert!((c.score[0] - 0.006).abs() < 1e-3 && (c.score[6] - 82.689).abs() < 1e-2);
        assert!((posterior(0.5, 1.0) - 0.5).abs() < 1e-12);
        assert!((posterior(0.5, 9.0) - 0.9).abs() < 1e-12);
    }

    #[test]
    fn score_rises_with_support_and_pools_reversals() {
        // Trusted variants concentrate at high support; one reversal between bins 2 and 3.
        let c = calibrate([7000, 1500, 800, 400, 200, 60, 40], [2, 1, 6, 2, 20, 60, 90]).unwrap();
        assert!(c.score.windows(2).all(|w| w[1] >= w[0]));
        assert_eq!(c.group[2], c.group[3]);
        assert!(c.score[0] < 0.1 && c.score[6] > 50.0);
        assert!(calibrate([100; BINS], [1; BINS]).is_none());
    }
}

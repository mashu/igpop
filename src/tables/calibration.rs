//! The two calibrations of a run as tables.

use super::num;
use crate::assess::{HaplotypeCheck, ScoreModel};
use crate::score;
use anyhow::Result;
use std::fmt::Write;

/// `score_levels.tsv`: the two groups counted per level of support, and
/// the score at the edges of each level. The score itself is the curve of
/// `score_model.tsv`; the levels are a ruler for reading it.
pub fn levels_table(model: &ScoreModel) -> Result<String> {
    let mut out = String::new();
    writeln!(
        out,
        "support\tarbitrary_changes\tarbitrary_fraction\ttrusted_variants\ttrusted_fraction\tratio_of_fractions\tscore_from\tscore_to\tscore_source"
    )?;
    let Some(calibration) = &model.calibration else { return Ok(out) };
    let levels = calibration.levels();
    let (na, nt) = (calibration.counts.arbitrary, calibration.counts.trusted);
    let fraction = |x: u64, n: u64| {
        if n == 0 {
            ".".to_string()
        } else {
            format!("{:.5}", x as f64 / n as f64)
        }
    };
    for (i, (arbitrary, trusted)) in levels.iter().enumerate() {
        let ratio = if na == 0 || nt == 0 || *arbitrary == 0 {
            ".".to_string()
        } else {
            format!("{:.3}", (*trusted as f64 / nt as f64) / (*arbitrary as f64 / na as f64))
        };
        let (from, to) = match i {
            0 => (calibration.never, calibration.never),
            _ => {
                let (low, high) = score::level_range(i);
                (calibration.score(low), calibration.score(high))
            }
        };
        writeln!(
            out,
            "{}\t{arbitrary}\t{}\t{trusted}\t{}\t{ratio}\t{from:.3}\t{to:.3}\t{}",
            score::LABELS[i],
            fraction(*arbitrary, na),
            fraction(*trusted, nt),
            model.source.as_str()
        )?;
    }
    Ok(out)
}

/// `haplotype_calibration.tsv`: how often an arbitrary single-base change
/// of a real chromosome is itself found on at least so many chromosomes.
pub fn haplotype_table(check: &HaplotypeCheck, limit: f64) -> Result<String> {
    let mut out = String::new();
    writeln!(out, "chromosomes\tcoincidence\tsupports")?;
    if let Some(coincidence) = &check.coincidence {
        for count in [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200, 500, 1000] {
            let rate = coincidence.rate(count);
            let supports = if count >= 2 && rate <= limit { "yes" } else { "no" };
            writeln!(out, "{count}\t{}\t{supports}", num(rate))?;
        }
    }
    Ok(out)
}

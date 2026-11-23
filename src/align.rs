//! Semi-global affine-gap alignment (Gotoh): the query (allele) is aligned
//! end to end, the reference window has free leading and trailing overhangs.
//! No seeding or mapping is involved: the locus is given by the caller.

#[derive(Clone, Copy, Debug)]
pub struct Scoring {
    pub matched: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Scoring {
            matched: 2,
            mismatch: -3,
            gap_open: -5,
            gap_extend: -2,
        }
    }
}

impl Scoring {
    #[inline]
    fn sub(&self, q: u8, r: u8) -> i32 {
        if q == b'N' || r == b'N' {
            0
        } else if q == r {
            self.matched
        } else {
            self.mismatch
        }
    }
}

/// One alignment column.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// Query base `q` against window base `r` (match or mismatch).
    Pair { q: usize, r: usize },
    /// Query base `q` with no reference base (insertion in the allele).
    Ins { q: usize },
    /// Window base `r` with no query base (deletion in the allele).
    Del { r: usize },
}

#[derive(Clone, Debug)]
pub struct Alignment {
    pub score: i32,
    pub ops: Vec<Op>,
    /// Half-open window interval covered by the query.
    pub ref_start: usize,
    pub ref_end: usize,
    pub matches: usize,
    pub columns: usize,
}

impl Alignment {
    pub fn identity(&self) -> f64 {
        if self.columns == 0 {
            0.0
        } else {
            self.matches as f64 / self.columns as f64
        }
    }
}

const NEG: i32 = i32::MIN / 4;

#[inline]
fn max3(a: i32, b: i32, c: i32) -> i32 {
    a.max(b).max(c)
}

/// Score only, O(n) memory. Used to compare an allele against other loci,
/// which is most of the alignment work of a run, so the rows are laid out
/// for the compiler to vectorise: everything that depends only on the
/// previous row is computed in one pass over plain slices, and only the
/// gap running along the row is carried from cell to cell.
pub fn score_only(query: &[u8], window: &[u8], s: &Scoring) -> i32 {
    let n = window.len();
    let (open, extend) = (s.gap_open + s.gap_extend, s.gap_extend);
    // Previous row: `best` = best of the three states, `gap_free` = best of
    // match and gap-in-window (the states a gap in the query can open from),
    // `up` = gap in the query.
    let mut best = vec![0i32; n + 1];
    let mut gap_free = vec![0i32; n + 1];
    let mut up = vec![NEG; n + 1];
    let (mut next_best, mut next_gap_free, mut next_up) = (vec![NEG; n + 1], vec![NEG; n + 1], vec![NEG; n + 1]);
    let mut diagonal = vec![NEG; n + 1];
    let mut substitution = vec![0i32; n];
    for (i, &q) in query.iter().enumerate() {
        for (cell, &r) in substitution.iter_mut().zip(window) {
            *cell = s.sub(q, r);
        }
        // From the previous row: the diagonal move and the gap in the query.
        for (j, (cell, sub)) in diagonal[1..].iter_mut().zip(&substitution).enumerate() {
            *cell = best[j] + sub;
        }
        for (cell, (free, gap)) in next_up[1..].iter_mut().zip(gap_free[1..].iter().zip(&up[1..])) {
            *cell = (free + open).max(gap + extend);
        }
        // Column 0: the whole query so far is one gap.
        next_up[0] = s.gap_open + s.gap_extend * (i as i32 + 1);
        next_best[0] = next_up[0];
        next_gap_free[0] = NEG;
        // Along the row: the gap in the window.
        let mut left = NEG;
        let mut opening = next_up[0];
        for j in 1..=n {
            left = (opening + open).max(left + extend);
            opening = diagonal[j].max(next_up[j]);
            next_gap_free[j] = diagonal[j].max(left);
            next_best[j] = opening.max(left);
        }
        std::mem::swap(&mut best, &mut next_best);
        std::mem::swap(&mut gap_free, &mut next_gap_free);
        std::mem::swap(&mut up, &mut next_up);
    }
    best.into_iter().max().unwrap_or(NEG)
}

/// Full alignment with traceback.
pub fn align(query: &[u8], window: &[u8], s: &Scoring) -> Alignment {
    let (m, n) = (query.len(), window.len());
    let w = n + 1;
    let oe = s.gap_open + s.gap_extend;
    let mut mm = vec![NEG; (m + 1) * w];
    let mut xx = vec![NEG; (m + 1) * w];
    let mut yy = vec![NEG; (m + 1) * w];
    mm[..=n].fill(0);
    for i in 1..=m {
        xx[i * w] = s.gap_open + s.gap_extend * i as i32;
        for j in 1..=n {
            let (c, u, l, d) = (i * w + j, (i - 1) * w + j, i * w + j - 1, (i - 1) * w + j - 1);
            mm[c] = max3(mm[d], xx[d], yy[d]) + s.sub(query[i - 1], window[j - 1]);
            xx[c] = max3(mm[u] + oe, xx[u] + s.gap_extend, yy[u] + oe);
            yy[c] = max3(mm[l] + oe, yy[l] + s.gap_extend, xx[l] + oe);
        }
    }

    #[derive(Clone, Copy, PartialEq)]
    enum St {
        M,
        X,
        Y,
    }
    let best_state = |i: usize, j: usize| -> (i32, St) {
        let c = i * w + j;
        let mut b = (mm[c], St::M);
        if xx[c] > b.0 {
            b = (xx[c], St::X);
        }
        if yy[c] > b.0 {
            b = (yy[c], St::Y);
        }
        b
    };

    // Free trailing overhang: best cell in the last row (leftmost on ties).
    let (mut j, mut score, mut st) = (0usize, NEG, St::M);
    for jj in 0..=n {
        let (sc, s_) = best_state(m, jj);
        if sc > score {
            score = sc;
            j = jj;
            st = s_;
        }
    }
    let ref_end = j;
    let mut i = m;
    let mut ops = Vec::with_capacity(m + 8);
    while i > 0 {
        let c = i * w + j;
        match st {
            St::M => {
                ops.push(Op::Pair { q: i - 1, r: j - 1 });
                let prev = mm[c] - s.sub(query[i - 1], window[j - 1]);
                let d = (i - 1) * w + j - 1;
                st = if prev == mm[d] {
                    St::M
                } else if prev == xx[d] {
                    St::X
                } else {
                    St::Y
                };
                i -= 1;
                j -= 1;
            }
            St::X => {
                ops.push(Op::Ins { q: i - 1 });
                if j == 0 {
                    i -= 1;
                    continue;
                }
                let u = (i - 1) * w + j;
                st = if xx[c] == mm[u] + oe {
                    St::M
                } else if xx[c] == xx[u] + s.gap_extend {
                    St::X
                } else {
                    St::Y
                };
                i -= 1;
            }
            St::Y => {
                ops.push(Op::Del { r: j - 1 });
                let l = i * w + j - 1;
                st = if yy[c] == mm[l] + oe {
                    St::M
                } else if yy[c] == yy[l] + s.gap_extend {
                    St::Y
                } else {
                    St::X
                };
                j -= 1;
            }
        }
    }
    ops.reverse();
    let matches = ops
        .iter()
        .filter(|o| matches!(o, Op::Pair { q, r } if query[*q] == window[*r]))
        .count();
    Alignment {
        score,
        columns: ops.len(),
        ops,
        ref_start: j,
        ref_end,
        matches,
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The textbook recurrence, cell by cell: what `score_only` must equal.
    fn score_reference(query: &[u8], window: &[u8], s: &Scoring) -> i32 {
        let (m, n) = (query.len(), window.len());
        let oe = s.gap_open + s.gap_extend;
        let mut pm = vec![0i32; n + 1];
        let mut px = vec![NEG; n + 1];
        let mut py = vec![NEG; n + 1];
        let (mut cm, mut cx, mut cy) = (vec![NEG; n + 1], vec![NEG; n + 1], vec![NEG; n + 1]);
        for i in 1..=m {
            cm[0] = NEG;
            cy[0] = NEG;
            cx[0] = s.gap_open + s.gap_extend * i as i32;
            for j in 1..=n {
                cm[j] = max3(pm[j - 1], px[j - 1], py[j - 1]) + s.sub(query[i - 1], window[j - 1]);
                cx[j] = max3(pm[j] + oe, px[j] + s.gap_extend, py[j] + oe);
                cy[j] = max3(cm[j - 1] + oe, cy[j - 1] + s.gap_extend, cx[j - 1] + oe);
            }
            std::mem::swap(&mut pm, &mut cm);
            std::mem::swap(&mut px, &mut cx);
            std::mem::swap(&mut py, &mut cy);
        }
        (0..=n).map(|j| max3(pm[j], px[j], py[j])).max().unwrap_or(NEG)
    }

    /// Deterministic pseudo-random sequence over ACGT with a few Ns.
    fn sequence(seed: &mut u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGTACGTACGTACGTN"[((*seed >> 33) % 17) as usize]
            })
            .collect()
    }

    #[test]
    fn fast_score_equals_the_reference_recurrence() {
        let s = Scoring::default();
        let mut seed = 7;
        for case in 0..300 {
            let window = sequence(&mut seed, 20 + case % 90);
            // Unrelated queries, and queries cut from the window with edits.
            let mut query = if case % 2 == 0 {
                sequence(&mut seed, 1 + case % 60)
            } else {
                window[case % 10..window.len() - case % 7].to_vec()
            };
            if case % 4 == 1 && query.len() > 12 {
                query.remove(5);
                query.insert(9, b'T');
                query[2] = b'G';
                query.drain(7..7 + case % 4);
            }
            assert_eq!(score_only(&query, &window, &s), score_reference(&query, &window, &s), "case {case}");
            assert_eq!(score_only(&query, &window, &s), align(&query, &window, &s).score, "case {case}");
        }
    }

    #[test]
    fn exact_match_inside_window() {
        let w = b"TTTTACGTACGGATCCATGCAAAA";
        let q = b"ACGTACGGATCCATGC";
        let a = align(q, w, &Scoring::default());
        assert_eq!((a.ref_start, a.ref_end), (4, 20));
        assert_eq!(a.matches, q.len());
        assert_eq!(a.score, score_only(q, w, &Scoring::default()));
    }

    #[test]
    fn deletion_and_score_agree() {
        let w = b"GGGGACGTACGGATCCATGCATTAGCCGGGG";
        let q = b"ACGTACGGATATGCATTAGCC"; // "CC" deleted
        let a = align(q, w, &Scoring::default());
        assert_eq!(a.ops.iter().filter(|o| matches!(o, Op::Del { .. })).count(), 2);
        assert_eq!(a.score, score_only(q, w, &Scoring::default()));
    }

    #[test]
    fn revcomp_roundtrip() {
        assert_eq!(revcomp(b"AACGT"), b"ACGTT");
    }
}

//! Conversion of an alignment into normalised, VCF-style variants:
//! atomised SNVs, and indels that are left-aligned and parsimonious
//! (Tan, Abecasis & Kang 2015), so they compare directly with split,
//! normalised VCF records.

use crate::align::{Alignment, Op};
use anyhow::{bail, Result};

/// A variant as population files name it: 1-based position, REF and ALT,
/// trimmed to the shortest form (see `trim`).
pub type Key = (u64, Vec<u8>, Vec<u8>);

/// A variant on a contig.
pub type Lookup = (String, Key);

/// The key of one ALT allele of a VCF record: bases shared by REF and ALT
/// are trimmed, the end first, then the start, keeping one. Records with
/// several ALT alleles are split this way, so that the same change has the
/// same key whatever else the record holds.
pub fn trim(mut pos: u64, reference: &[u8], alternate: &[u8]) -> Key {
    let (mut r, mut a) = (reference.to_vec(), alternate.to_vec());
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    (pos, r, a)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Variant {
    /// 0-based offset of the first REF base inside the window.
    pub offset: usize,
    pub reference: Vec<u8>,
    pub alternate: Vec<u8>,
    /// 0-based index, in the aligned (window-strand) query, of the first
    /// query base involved (for a deletion: the base following it).
    pub query_index: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Snv,
    Ins,
    Del,
}

impl Kind {
    pub fn as_str(self) -> &'static str {
        match self {
            Kind::Snv => "SNV",
            Kind::Ins => "INS",
            Kind::Del => "DEL",
        }
    }
}

impl Variant {
    pub fn kind(&self) -> Kind {
        match (self.reference.len(), self.alternate.len()) {
            (1, 1) => Kind::Snv,
            (r, a) if r > a => Kind::Del,
            _ => Kind::Ins,
        }
    }
}

/// Left-align and trim one raw variant against the window sequence.
pub fn normalise(mut v: Variant, window: &[u8]) -> Result<Variant> {
    loop {
        let mut changed = false;
        if let (Some(a), Some(b)) = (v.reference.last(), v.alternate.last()) {
            if a == b {
                v.reference.pop();
                v.alternate.pop();
                changed = true;
            }
        }
        if v.reference.is_empty() || v.alternate.is_empty() {
            if v.offset == 0 {
                bail!("indel reaches the window edge; increase --flank");
            }
            v.offset -= 1;
            v.reference.insert(0, window[v.offset]);
            v.alternate.insert(0, window[v.offset]);
            changed = true;
        }
        if !changed {
            break;
        }
    }
    while v.reference.len() > 1 && v.alternate.len() > 1 && v.reference[0] == v.alternate[0] {
        v.reference.remove(0);
        v.alternate.remove(0);
        v.offset += 1;
    }
    Ok(v)
}

/// All differences between the query and the window, normalised and sorted.
/// N in the query is "unknown", never a variant.
pub fn call(aln: &Alignment, query: &[u8], window: &[u8]) -> Result<Vec<Variant>> {
    let mut out = Vec::new();
    let mut next_ref = aln.ref_start; // window offset of the next unconsumed base
    let mut k = 0;
    while k < aln.ops.len() {
        match aln.ops[k] {
            Op::Pair { q, r } => {
                if query[q] != window[r] && query[q] != b'N' && window[r] != b'N' {
                    out.push(Variant {
                        offset: r,
                        reference: vec![window[r]],
                        alternate: vec![query[q]],
                        query_index: q,
                    });
                }
                next_ref = r + 1;
                k += 1;
            }
            Op::Ins { q } => {
                let mut alt = Vec::new();
                while let Some(Op::Ins { q: qq }) = aln.ops.get(k) {
                    alt.push(query[*qq]);
                    k += 1;
                }
                let raw = Variant {
                    offset: next_ref,
                    reference: vec![],
                    alternate: alt,
                    query_index: q,
                };
                out.push(normalise(raw, window)?);
            }
            Op::Del { r } => {
                let mut del = Vec::new();
                while let Some(Op::Del { r: rr }) = aln.ops.get(k) {
                    del.push(window[*rr]);
                    next_ref = *rr + 1;
                    k += 1;
                }
                let q = match aln.ops.get(k) {
                    Some(Op::Pair { q, .. }) | Some(Op::Ins { q }) => *q,
                    _ => query.len(),
                };
                let raw = Variant {
                    offset: r,
                    reference: del,
                    alternate: vec![],
                    query_index: q,
                };
                out.push(normalise(raw, window)?);
            }
        }
    }
    out.sort_by(|a, b| a.offset.cmp(&b.offset).then(a.reference.len().cmp(&b.reference.len())));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trimming_matches_normalised_form() {
        assert_eq!(trim(100, b"GAC", b"G"), (100, b"GAC".to_vec(), b"G".to_vec()));
        assert_eq!(trim(100, b"GAC", b"GAT"), (102, b"C".to_vec(), b"T".to_vec()));
        assert_eq!(trim(100, b"CAA", b"CA"), (100, b"CA".to_vec(), b"C".to_vec()));
    }
    use crate::align::{align, Scoring};

    fn variants(q: &[u8], w: &[u8]) -> Vec<(usize, String, String)> {
        let a = align(q, w, &Scoring::default());
        call(&a, q, w)
            .unwrap()
            .into_iter()
            .map(|v| {
                (
                    v.offset,
                    String::from_utf8(v.reference).unwrap(),
                    String::from_utf8(v.alternate).unwrap(),
                )
            })
            .collect()
    }

    #[test]
    fn snv() {
        let w = b"TTTTACGTACGGATCCATGCATTAGCAAAA";
        let q = b"ACGTACGGATGCATGCATTAGC";
        assert_eq!(variants(q, w), vec![(14, "C".into(), "G".into())]);
    }

    #[test]
    fn deletion_in_repeat_is_left_aligned() {
        //            0123456789
        let w = b"GGCTGACACACACTTGACCGTAGG";
        let q = b"CTGACACACTTGACCGTA"; // one AC unit removed
                                       // leftmost representation: anchor G at offset 4, delete "AC"
        assert_eq!(variants(q, w), vec![(4, "GAC".into(), "G".into())]);
    }

    #[test]
    fn insertion_in_homopolymer_is_left_aligned() {
        let w = b"GGCTGCAAAATCGGACCGTAGG";
        let q = b"CTGCAAAAATCGGACCGTA"; // one extra A
        assert_eq!(variants(q, w), vec![(5, "C".into(), "CA".into())]);
    }

    #[test]
    fn n_is_not_a_variant() {
        let w = b"TTTTACGTACGGATCCATGCATTAGCAAAA";
        let q = b"ACGTACGGATNCATGCATTAGC";
        assert!(variants(q, w).is_empty());
    }
}

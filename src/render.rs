//! Human-readable alignments, one block per gene: the reference on top and
//! every allele of that gene stacked under it, in gene (5'->3') orientation.

use crate::align::{Alignment, Op};

/// One allele laid over its locus window, in reference-strand coordinates:
/// `.` match, a base for a mismatch, `-` for a deleted reference base, and a
/// space where the allele does not reach.
pub fn row(aln: &Alignment, query: &[u8], window: &[u8]) -> Vec<u8> {
    let mut row = vec![b' '; window.len()];
    for op in &aln.ops {
        match *op {
            Op::Pair { q, r } => row[r] = if query[q] == window[r] { b'.' } else { query[q] },
            Op::Del { r } => row[r] = b'-',
            Op::Ins { .. } => {}
        }
    }
    row
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

pub struct Entry<'a> {
    pub name: &'a str,
    pub row: &'a [u8],
    /// Free text after the alignment: identity, status, insertions.
    pub note: String,
}

/// Text block for one gene. `window_start` is the 1-based reference
/// position of the first window base; `reverse` shows the minus strand.
pub fn gene_block(title: &str, window: &[u8], window_start: u64, reverse: bool, entries: &[Entry], width: usize) -> String {
    let covered = |i: usize| entries.iter().any(|e| e.row[i] != b' ');
    let Some(first) = (0..window.len()).find(|i| covered(*i)) else {
        return format!(">{title}\n  (no allele aligned)\n\n");
    };
    let last = (0..window.len()).rev().find(|i| covered(*i)).unwrap_or(first);
    // Display columns as window offsets, in display order.
    let columns: Vec<usize> = if reverse {
        (first..=last).rev().collect()
    } else {
        (first..=last).collect()
    };
    let show = |b: u8| if reverse { complement(b) } else { b } as char;
    let label = entries.iter().map(|e| e.name.len()).max().unwrap_or(0).max("reference".len());

    let mut out = format!(">{title}\n");
    for chunk in columns.chunks(width) {
        let (a, b) = (window_start + chunk[0] as u64, window_start + *chunk.last().unwrap() as u64);
        out.push_str(&format!("{:label$}  {a} .. {b}\n", ""));
        out.push_str(&format!(
            "{:label$}  {}\n",
            "reference",
            chunk.iter().map(|i| show(window[*i])).collect::<String>()
        ));
        for e in entries {
            let line: String = chunk.iter().map(|i| show(e.row[*i])).collect();
            if line.trim().is_empty() {
                continue;
            }
            out.push_str(&format!("{:label$}  {}\n", e.name, line.trim_end()));
        }
        out.push('\n');
    }
    for e in entries {
        out.push_str(&format!("  {:label$}  {}\n", e.name, e.note));
    }
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{align, Scoring};

    #[test]
    fn stacked_rows_show_only_differences() {
        let window = b"TTTTACGTACGGATCCATGCATTAGCAAAA";
        let query = b"ACGTACGGATGCATGCATTAGC";
        let aln = align(query, window, &Scoring::default());
        let r = row(&aln, query, window);
        assert_eq!(std::str::from_utf8(&r).unwrap(), "    ..........G...........    ");
        let block = gene_block(
            "G",
            window,
            101,
            false,
            &[Entry {
                name: "G*01",
                row: &r,
                note: "ok".into(),
            }],
            100,
        );
        assert!(block.contains("reference  ACGTACGGATCCATGCATTAGC"));
        assert!(block.contains("G*01       ..........G"));
        assert!(block.contains("105 .. 126"));
    }

    #[test]
    fn minus_strand_is_shown_in_gene_orientation() {
        let window = b"AACCGGTT";
        let r = b"  ...A  ".to_vec();
        let block = gene_block(
            "G",
            window,
            1,
            true,
            &[Entry {
                name: "a",
                row: &r,
                note: String::new(),
            }],
            100,
        );
        // columns 5..2 reversed and complemented: window CCGG -> CCGG, row "...A" -> "T..."
        assert!(block.contains("reference  CCGG"));
        assert!(block.contains("a          T..."));
    }
}

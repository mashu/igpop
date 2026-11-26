//! Where each allele was placed: `alignments.tsv`, and `alignments.txt`
//! to read by eye.

use super::Table;
use crate::locus::{Placed, Placement, Window};
use crate::render;
use anyhow::Result;
use std::fmt::Write;

/// One row of `alignments.tsv`.
pub struct PlacedRow<'a> {
    pub build: &'a str,
    pub placed: &'a Placed,
}

/// `alignments.tsv`: where each allele was placed.
pub fn alignments_table<'a>() -> Table<PlacedRow<'a>> {
    type R<'a> = PlacedRow<'a>;
    let placement = |cell: fn(&Placement) -> String| move |r: &R| r.placed.placement.as_ref().map_or(".".to_string(), cell);
    Table::new()
        .col("allele", |r: &R| r.placed.allele.clone())
        .col("gene", |r: &R| r.placed.gene.clone())
        .col("build", |r: &R| r.build.to_string())
        .col("contig", |r: &R| r.placed.contig().to_string())
        .col("ref_start", placement(|a| a.span.0.to_string()))
        .col("ref_end", placement(|a| a.span.1.to_string()))
        .col("strand", placement(|a| if a.reverse { "-" } else { "+" }.to_string()))
        .col("identity", placement(|a| format!("{:.4}", a.identity)))
        .col("score", placement(|a| a.score.to_string()))
        .col(
            "best_other_gene",
            placement(|a| a.rival.as_ref().map_or(".".to_string(), |(gene, _)| gene.clone())),
        )
        .col(
            "best_other_score",
            placement(|a| a.rival.as_ref().map_or(".".to_string(), |(_, s)| s.to_string())),
        )
        .col("shares_locus_with", |r: &R| r.placed.shared())
        .col("status", |r: &R| r.placed.status.as_str().to_string())
        .col("missing_5prime", placement(|a| a.missing.0.to_string()))
        .col("missing_3prime", placement(|a| a.missing.1.to_string()))
}

/// The head of `alignments.txt`: how to read the blocks.
pub fn alignments_legend() -> String {
    "# igpop alignments, one block per gene, shown 5'->3' in gene orientation.\n\
     # '.' = same as reference, letter = mismatch, '-' = reference base missing in the allele,\n\
     # blank = not covered by the allele. Inserted bases are listed under each block as\n\
     # <reference position before the insertion>+<bases on the reference strand>.\n\
     # Coordinates above each chunk are the reference positions of its first and last column.\n\n"
        .to_string()
}

/// The blocks of `alignments.txt` for the loci of one build: stacked
/// alignments, one block per gene, to read by eye.
pub fn alignments_text(build: &str, windows: &[Window], placed: &[Placed]) -> Result<String> {
    let mut out = String::new();
    let mut order: Vec<usize> = (0..windows.len()).collect();
    order.sort_by(|x, y| windows[*x].coord.gene.cmp(&windows[*y].coord.gene));
    for index in order {
        let members: Vec<(&Placed, &Placement)> = placed
            .iter()
            .filter_map(|p| p.placement.as_ref().filter(|a| a.window == index).map(|a| (p, a)))
            .collect();
        let Some((_, first)) = members.first() else { continue };
        let w = &windows[index];
        let reverse = first.reverse;
        let title = format!(
            "{}  {build} {}:{}-{}  strand {}  {} allele(s)",
            w.coord.gene,
            w.coord.contig,
            w.coord.start,
            w.coord.end,
            if reverse { '-' } else { '+' },
            members.len(),
        );
        let entries: Vec<render::Entry> = members
            .iter()
            .map(|(p, a)| render::Entry {
                name: &p.allele,
                row: &a.row,
                note: format!(
                    "identity {:.4}  status {}  variants {}  insertions {}{}",
                    a.identity,
                    p.status.as_str(),
                    p.variants.len(),
                    if a.insertions.is_empty() {
                        "none".to_string()
                    } else {
                        a.insertions.join(",")
                    },
                    if a.reverse == reverse {
                        ""
                    } else {
                        "  (aligned on the other strand)"
                    },
                ),
            })
            .collect();
        write!(out, "{}", render::gene_block(&title, &w.seq, w.start, reverse, &entries, 100))?;
    }
    Ok(out)
}

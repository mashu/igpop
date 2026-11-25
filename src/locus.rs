//! Step 1: place every allele on the locus of its gene and call its
//! differences from the reference.

use crate::align::{align, revcomp, score_only, Alignment, Op, Scoring};
use crate::fetch::merge_regions;
use crate::freq::Key;
use crate::io::{Coord, Reference, Strand};
use crate::variant::Kind;
use crate::{par, render, variant};
use anyhow::{bail, Context, Result};
use std::collections::HashMap;

/// A gene interval plus flanks, with its reference sequence.
pub struct Window {
    pub coord: Coord,
    /// 1-based reference position of the first window base.
    pub start: u64,
    pub seq: Vec<u8>,
}

impl Window {
    pub fn end(&self) -> u64 {
        self.start + self.seq.len() as u64 - 1
    }

    fn same_locus(&self, other: &Window) -> bool {
        self.coord.contig == other.coord.contig && self.coord.start == other.coord.start && self.coord.end == other.coord.end
    }

    /// The query in every orientation this locus allows, as (is_reverse, sequence).
    fn orientations(&self, seq: &[u8]) -> Vec<(bool, Vec<u8>)> {
        match self.coord.strand {
            Strand::Forward => vec![(false, seq.to_vec())],
            Strand::Reverse => vec![(true, revcomp(seq))],
            Strand::Auto => vec![(false, seq.to_vec()), (true, revcomp(seq))],
        }
    }
}

/// Reference windows for all loci. Neighbouring loci share one read, so a
/// gene-dense locus costs a single range request.
pub fn load_windows(coords: Vec<Coord>, reference: &mut Reference, flank: u64) -> Result<Vec<Window>> {
    let mut spans: Vec<(Coord, u64, u64)> = Vec::with_capacity(coords.len());
    for mut coord in coords {
        if spans.iter().any(|s| s.0.gene == coord.gene) {
            bail!("gene {} has more than one row for build {}", coord.gene, coord.build);
        }
        coord.contig = reference.resolve_contig(&coord.contig)?;
        let len = reference.contig_length(&coord.contig)?;
        let start = coord.start.saturating_sub(flank).max(1);
        let end = (coord.end + flank).min(len);
        spans.push((coord, start, end));
    }
    spans.sort_by(|a, b| (&a.0.contig, a.1).cmp(&(&b.0.contig, b.1)));
    let blocks = merge_regions(spans.iter().map(|s| (s.0.contig.clone(), s.1, s.2)).collect(), 250_000);
    let mut windows: Vec<Window> = Vec::with_capacity(spans.len());
    for (contig, block_start, block_end) in &blocks {
        let block = reference
            .fetch(contig, *block_start, *block_end)
            .with_context(|| format!("fetching {contig}:{block_start}-{block_end}"))?;
        for (coord, start, end) in spans
            .iter()
            .filter(|s| &s.0.contig == contig && s.1 >= *block_start && s.2 <= *block_end)
        {
            let seq = block[(start - block_start) as usize..=(end - block_start) as usize].to_vec();
            windows.push(Window {
                coord: coord.clone(),
                start: *start,
                seq,
            });
        }
    }
    Ok(windows)
}

/// `allele<TAB>gene` lines; `#` lines are skipped.
pub fn parse_aliases(text: &str, label: &str) -> Result<HashMap<String, String>> {
    let mut aliases = HashMap::new();
    for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        match line.split_once('\t') {
            Some((allele, gene)) => aliases.insert(allele.trim().to_string(), gene.trim().to_string()),
            None => bail!("{label}: expected allele<TAB>gene, found '{line}'"),
        };
    }
    Ok(aliases)
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// Aligned to the locus of its gene.
    Ok,
    /// Another locus fits at least as well.
    AmbiguousLocus,
    /// Too different from its locus; no variants are called.
    LowIdentity,
    /// Its gene has no coordinates on this build.
    NoCoordinates,
}

impl Status {
    pub const ALL: [Status; 4] = [Status::Ok, Status::AmbiguousLocus, Status::LowIdentity, Status::NoCoordinates];

    pub fn as_str(self) -> &'static str {
        match self {
            Status::Ok => "ok",
            Status::AmbiguousLocus => "ambiguous_locus",
            Status::LowIdentity => "low_identity",
            Status::NoCoordinates => "no_coordinates",
        }
    }
}

/// One difference between an allele and the reference.
pub struct Called {
    pub key: Key,
    pub kind: Kind,
    /// 1-based position in the allele as supplied.
    pub allele_pos: usize,
    /// The same base counted along the forward strand of the genome: for
    /// a gene on the reverse strand, from the other end of the allele.
    pub strand_pos: usize,
    /// The gene lies on the reverse strand.
    pub reverse: bool,
}

/// Where and how well an allele aligned.
pub struct Placement {
    /// Index of its window.
    pub window: usize,
    pub contig: String,
    /// First and last reference position covered by the allele.
    pub span: (u64, u64),
    pub reverse: bool,
    pub identity: f64,
    pub score: i32,
    /// Best competing locus and the allele's score there.
    pub rival: Option<(String, i32)>,
    /// Other genes annotated at exactly the same interval.
    pub shared: Vec<String>,
    /// Gene bases the allele does not reach at its 5' and 3' ends.
    pub missing: (u64, u64),
    /// Display row over the window, and inserted bases as `pos+SEQ`.
    pub row: Vec<u8>,
    pub insertions: Vec<String>,
}

pub struct Placed {
    pub allele: String,
    pub gene: String,
    pub status: Status,
    /// `None` when the gene has no coordinates.
    pub placement: Option<Placement>,
    pub variants: Vec<Called>,
}

impl Placed {
    pub fn contig(&self) -> &str {
        self.placement.as_ref().map_or(".", |p| p.contig.as_str())
    }

    pub fn shared(&self) -> String {
        match &self.placement {
            Some(p) if !p.shared.is_empty() => p.shared.join(","),
            _ => ".".to_string(),
        }
    }

    /// The variant as (contig, key), the form evidence is looked up by.
    pub fn lookup(&self, called: &Called) -> (String, Key) {
        (self.contig().to_string(), called.key.clone())
    }
}

pub fn gene_of(allele: &str) -> &str {
    allele.split('*').next().unwrap_or(allele)
}

/// The gene an allele belongs to. The alias table is asked first: for the
/// allele's own name, and otherwise for the longest name in the table that
/// the allele's name extends (`IGHV4-4*02_S1234` extends `IGHV4-4*02`,
/// and `IGHV4-4*05` extends `IGHV4-4`, so a table can map a whole gene).
/// Without a match the gene is the part of the name before `*`.
pub fn gene_for<'a>(allele: &'a str, aliases: &'a HashMap<String, String>) -> &'a str {
    if let Some(gene) = aliases.get(allele) {
        return gene;
    }
    aliases
        .iter()
        .filter(|(name, _)| {
            allele
                .strip_prefix(name.as_str())
                .is_some_and(|rest| rest.starts_with(['*', '_', '.']))
        })
        .max_by_key(|(name, _)| name.len())
        .map_or_else(|| gene_of(allele), |(_, gene)| gene.as_str())
}

/// The locus a sequence fits best, with its score there; of equal scores
/// the first locus wins.
pub fn best_locus<'a>(seq: &[u8], windows: &'a [Window], scoring: &Scoring) -> Option<(&'a str, i32)> {
    let mut best: Option<(&str, i32)> = None;
    for w in windows {
        let score = score_at(seq, w, scoring);
        if best.is_none_or(|(_, top)| score > top) {
            best = Some((w.coord.gene.as_str(), score));
        }
    }
    best
}

/// Inserted bases as `pos+SEQ`, where pos is the reference base before them.
fn insertions(aln: &Alignment, query: &[u8], window_start: u64) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut before = aln.ref_start as u64 + window_start - 1;
    let mut open = false;
    for op in &aln.ops {
        match *op {
            Op::Pair { r, .. } | Op::Del { r } => {
                before = window_start + r as u64;
                open = false;
            }
            Op::Ins { q } => {
                if !open {
                    out.push(format!("{before}+"));
                    open = true;
                }
                out.last_mut().unwrap().push(query[q] as char);
            }
        }
    }
    out
}

pub struct PlaceOptions {
    pub min_identity: f64,
    pub paralog_check: bool,
    pub scoring: Scoring,
}

/// Best score of a sequence at one window, over the orientations it allows.
fn score_at(seq: &[u8], window: &Window, scoring: &Scoring) -> i32 {
    window
        .orientations(seq)
        .iter()
        .map(|(_, q)| score_only(q, &window.seq, scoring))
        .max()
        .expect("at least one orientation")
}

/// An allele already compared with every locus, used to spare a similar
/// allele most of that work.
struct Guide<'a> {
    seq: &'a [u8],
    /// Its score at every window; `None` at its own locus.
    scores: &'a [Option<i32>],
}

/// How much higher `seq` can score than `guide`, at any window: every
/// substituted base changes the best score by at most one match-mismatch
/// difference, every inserted or deleted base by at most a match plus a
/// new gap.
fn slack(seq: &[u8], guide: &[u8], s: &Scoring) -> i32 {
    let a = align(seq, guide, s);
    let substituted = a
        .ops
        .iter()
        .filter(|op| matches!(op, Op::Pair { q, r } if seq[*q] != guide[*r]))
        .count();
    let gapped = a.ops.iter().filter(|op| !matches!(op, Op::Pair { .. })).count() + guide.len() - (a.ref_end - a.ref_start);
    let (high, low) = (s.matched.max(s.mismatch).max(0), s.matched.min(s.mismatch).min(0));
    (high - low) * substituted as i32 + (high - s.gap_open - s.gap_extend) * gapped as i32
}

/// The best competing locus of an allele: its score at every window that
/// is not its own locus; of equal scores, the first window wins. With a
/// guide, windows are tried from the most promising down, and the search
/// stops when no remaining window can reach the best score found: the
/// result is the same as trying them all. Returns the winner and, without
/// a guide, all scores.
fn rival(
    seq: &[u8],
    own: &Window,
    windows: &[Window],
    scoring: &Scoring,
    guide: Option<&Guide>,
) -> (Option<(usize, i32)>, Vec<Option<i32>>) {
    let better = |best: Option<(usize, i32)>, i: usize, score: i32| match best {
        Some((at, top)) if score < top || (score == top && at < i) => Some((at, top)),
        _ => Some((i, score)),
    };
    let mut best = None;
    match guide {
        None => {
            let scores: Vec<Option<i32>> = windows
                .iter()
                .map(|w| (!w.same_locus(own)).then(|| score_at(seq, w, scoring)))
                .collect();
            for (i, score) in scores.iter().enumerate() {
                if let Some(score) = score {
                    best = better(best, i, *score);
                }
            }
            (best, scores)
        }
        Some(guide) => {
            let slack = slack(seq, guide.seq, scoring);
            let mut order: Vec<(i32, usize)> = guide
                .scores
                .iter()
                .enumerate()
                .filter_map(|(i, s)| s.map(|s| (s + slack, i)))
                .collect();
            order.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            for (bound, i) in order {
                if best.is_some_and(|(_, top)| bound < top) {
                    break;
                }
                best = better(best, i, score_at(seq, &windows[i], scoring));
            }
            (best, Vec::new())
        }
    }
}

/// Place every allele. The first allele of each gene is compared with all
/// loci; the others use it as their guide.
pub fn place_all(
    alleles: &[(String, Vec<u8>)],
    windows: &[Window],
    aliases: &HashMap<String, String>,
    options: &PlaceOptions,
) -> Result<Vec<Placed>> {
    let window_of = |name: &str| {
        let gene = gene_for(name, aliases);
        windows.iter().position(|w| w.coord.gene == gene)
    };
    let mut first: HashMap<usize, usize> = HashMap::new();
    let (mut leaders, mut followers) = (Vec::new(), Vec::new());
    for (i, (name, _)) in alleles.iter().enumerate() {
        match window_of(name) {
            Some(w) if first.contains_key(&w) => followers.push((i, first[&w])),
            Some(w) => {
                first.insert(w, leaders.len());
                leaders.push(i);
            }
            None => leaders.push(i),
        }
    }
    let led: Vec<(Placed, Vec<Option<i32>>)> = par::map(&leaders, None, |i| {
        place(&alleles[*i].0, &alleles[*i].1, windows, aliases, options, None)
    })
    .into_iter()
    .collect::<Result<_>>()?;
    let followed: Vec<Placed> = par::map(&followers, None, |(i, leader)| {
        let guide = Guide {
            seq: &alleles[leaders[*leader]].1,
            scores: &led[*leader].1,
        };
        place(&alleles[*i].0, &alleles[*i].1, windows, aliases, options, Some(&guide)).map(|(placed, _)| placed)
    })
    .into_iter()
    .collect::<Result<_>>()?;
    Ok(led.into_iter().map(|(placed, _)| placed).chain(followed).collect())
}

/// Align one allele to the window of its gene (in the orientation that
/// scores best), compare with every other locus, and call its variants.
/// Also returns its score at every window when no guide was given.
fn place(
    name: &str,
    seq: &[u8],
    windows: &[Window],
    aliases: &HashMap<String, String>,
    options: &PlaceOptions,
    guide: Option<&Guide>,
) -> Result<(Placed, Vec<Option<i32>>)> {
    let gene = gene_for(name, aliases);
    let Some(index) = windows.iter().position(|w| w.coord.gene == gene) else {
        let placed = Placed {
            allele: name.into(),
            gene: gene.into(),
            status: Status::NoCoordinates,
            placement: None,
            variants: vec![],
        };
        return Ok((placed, Vec::new()));
    };
    let own = &windows[index];
    let (reverse, query, aln): (bool, Vec<u8>, Alignment) = own
        .orientations(seq)
        .into_iter()
        .map(|(rev, q)| {
            let a = align(&q, &own.seq, &options.scoring);
            (rev, q, a)
        })
        .max_by_key(|(_, _, a)| a.score)
        .expect("at least one orientation");

    // Genes annotated at exactly the same interval are one locus, not rivals.
    let shared: Vec<String> = windows
        .iter()
        .filter(|w| w.coord.gene != gene && w.same_locus(own))
        .map(|w| w.coord.gene.clone())
        .collect();

    let (rival, scores) = match options.paralog_check {
        true => rival(seq, own, windows, &options.scoring, guide),
        false => (None, Vec::new()),
    };
    let rival = rival.map(|(i, score)| (windows[i].coord.gene.clone(), score));

    let status = if aln.identity() < options.min_identity {
        Status::LowIdentity
    } else if rival.as_ref().is_some_and(|(_, s)| *s >= aln.score) {
        Status::AmbiguousLocus
    } else {
        Status::Ok
    };
    let c = &own.coord;
    // Gene bases beyond either end of the allele, in gene orientation.
    let span = (own.start + aln.ref_start as u64, own.start + aln.ref_end as u64 - 1);
    let (low, high) = (span.0.saturating_sub(c.start), c.end.saturating_sub(span.1));
    let missing = if reverse { (high, low) } else { (low, high) };

    let mut variants = Vec::new();
    if status != Status::LowIdentity {
        for v in variant::call(&aln, &query, &own.seq).with_context(|| format!("allele {name}"))? {
            variants.push(Called {
                kind: v.kind(),
                allele_pos: if reverse { query.len() - v.query_index } else { v.query_index + 1 },
                strand_pos: v.query_index + 1,
                reverse,
                key: (own.start + v.offset as u64, v.reference, v.alternate),
            });
        }
    }
    let placed = Placed {
        allele: name.into(),
        gene: gene.into(),
        status,
        placement: Some(Placement {
            window: index,
            contig: c.contig.clone(),
            span,
            reverse,
            identity: aln.identity(),
            score: aln.score,
            rival,
            shared,
            missing,
            row: render::row(&aln, &query, &own.seq),
            insertions: insertions(&aln, &query, own.start),
        }),
        variants,
    };
    Ok((placed, scores))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aliases_reach_candidates_and_whole_genes() {
        let aliases: HashMap<String, String> = [("IGHV4-4*02", "IGHV4-4B"), ("IGHV9-9", "IGHV3-30")]
            .iter()
            .map(|(a, g)| (a.to_string(), g.to_string()))
            .collect();
        assert_eq!(gene_for("IGHV4-4*02", &aliases), "IGHV4-4B");
        assert_eq!(gene_for("IGHV4-4*02_S1234", &aliases), "IGHV4-4B");
        assert_eq!(gene_for("IGHV4-4*05", &aliases), "IGHV4-4");
        assert_eq!(gene_for("IGHV9-9*01_S1", &aliases), "IGHV3-30");
        assert_eq!(gene_for("IGHV9-9-1*01", &aliases), "IGHV9-9-1");
    }
}

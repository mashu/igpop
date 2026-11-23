//! Step 1: place every allele on the locus of its gene and call its
//! differences from the reference.

use crate::align::{align, revcomp, score_only, Alignment, Op, Scoring};
use crate::fetch::merge_regions;
use crate::freq::Key;
use crate::io::{Coord, Reference, Strand};
use crate::variant::Kind;
use crate::{render, variant};
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

/// Align one allele to the window of its gene (in the orientation that
/// scores best), compare with every other locus, and call its variants.
pub fn place(name: &str, seq: &[u8], windows: &[Window], aliases: &HashMap<String, String>, options: &PlaceOptions) -> Result<Placed> {
    let gene = aliases.get(name).map(String::as_str).unwrap_or_else(|| gene_of(name));
    let Some(index) = windows.iter().position(|w| w.coord.gene == gene) else {
        return Ok(Placed {
            allele: name.into(),
            gene: gene.into(),
            status: Status::NoCoordinates,
            placement: None,
            variants: vec![],
        });
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

    // Best competing locus: the same allele scored at every other window.
    let mut rival: Option<(String, i32)> = None;
    if options.paralog_check {
        for w in windows.iter().filter(|w| !w.same_locus(own)) {
            for (_, q) in w.orientations(seq) {
                let s = score_only(&q, &w.seq, &options.scoring);
                if rival.as_ref().is_none_or(|(_, best)| s > *best) {
                    rival = Some((w.coord.gene.clone(), s));
                }
            }
        }
    }

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
                key: (own.start + v.offset as u64, v.reference, v.alternate),
            });
        }
    }
    Ok(Placed {
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
    })
}

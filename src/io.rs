//! Input parsing (coordinates table, allele FASTA) and random access to a
//! reference FASTA through its `.fai` index. The reference may be a local
//! file or an http(s) URL, plain or bgzip-compressed (`.gz` with `.gzi`);
//! only the requested bytes are read or downloaded.

use crate::net;
use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    /// Not given: both orientations are aligned and the better one is kept.
    Auto,
}

/// One gene locus on one genome build. Coordinates are 1-based, inclusive.
#[derive(Clone, Debug)]
pub struct Coord {
    pub gene: String,
    pub build: String,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

/// Canonical build label: 37/hg19/b37 -> GRCh37, 38/hg38 -> GRCh38,
/// anything else is kept as written.
pub fn normalise_build(s: &str) -> String {
    match s.trim().to_ascii_lowercase().as_str() {
        "37" | "grch37" | "hg19" | "b37" | "hs37d5" => "GRCh37".into(),
        "38" | "grch38" | "hg38" => "GRCh38".into(),
        _ => s.trim().to_string(),
    }
}

/// Coordinates from tab-separated text: gene, build, contig, start, end,
/// strand (+, - or . for "try both"). `#` lines and a header line starting
/// with `gene` are skipped. `label` names the source in error messages.
pub fn parse_coords(text: &str, label: &str) -> Result<Vec<Coord>> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') || (n == 0 && line.to_lowercase().starts_with("gene")) {
            continue;
        }
        let at = format!("{label}:{}", n + 1);
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            bail!(
                "{at}: expected gene, build, contig, start, end, strand (tab-separated), found {} field(s)",
                f.len()
            );
        }
        let num = |s: &str| s.parse::<u64>().with_context(|| format!("{at}: bad coordinate '{s}'"));
        let coord = Coord {
            gene: f[0].into(),
            build: normalise_build(f[1]),
            contig: f[2].into(),
            start: num(f[3])?,
            end: num(f[4])?,
            strand: match f[5] {
                "+" => Strand::Forward,
                "-" => Strand::Reverse,
                "." | "" => Strand::Auto,
                other => bail!("{at}: strand must be +, - or ., found '{other}'"),
            },
        };
        if coord.start < 1 || coord.end < coord.start {
            bail!("{at}: need 1 <= start <= end");
        }
        out.push(coord);
    }
    Ok(out)
}

/// (name, sequence) per record; names end at the first blank, sequences
/// are upper-cased and stripped of everything but letters (IMGT gap dots).
pub fn parse_fasta(text: &str) -> Vec<(String, Vec<u8>)> {
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1
                .extend(line.bytes().filter(|b| b.is_ascii_alphabetic()).map(|b| b.to_ascii_uppercase()));
        }
    }
    out
}

struct FaiEntry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// A file read by byte range: local, or remote over HTTP.
pub enum Source {
    Local(File),
    Remote(String),
}

impl Source {
    /// A local path or an http(s) URL.
    pub fn open(spec: &str) -> Result<Self> {
        if net::is_url(spec) {
            Ok(Source::Remote(spec.to_string()))
        } else {
            Ok(Source::Local(File::open(spec).with_context(|| format!("opening {spec}"))?))
        }
    }

    /// Bytes `from..=to` of the file; `to = None` reads to the end. A
    /// range that runs past the end of the file is returned short.
    pub fn read(&mut self, from: u64, to: Option<u64>) -> Result<Vec<u8>> {
        match self {
            Source::Local(file) => {
                file.seek(SeekFrom::Start(from))?;
                let mut buf = Vec::new();
                match to {
                    Some(to) => file.take(to - from + 1).read_to_end(&mut buf)?,
                    None => file.read_to_end(&mut buf)?,
                };
                Ok(buf)
            }
            Source::Remote(url) => net::get_range(url, from, to),
        }
    }
}

/// A whole (small) file, e.g. an index: local path or URL.
pub fn read_whole(spec: &str) -> Result<Vec<u8>> {
    if net::is_url(spec) {
        net::get(spec)
    } else {
        std::fs::read(spec).with_context(|| format!("reading {spec}"))
    }
}

pub struct Reference {
    source: Source,
    index: HashMap<String, FaiEntry>,
    /// bgzip block starts as (compressed offset, uncompressed offset).
    blocks: Option<Vec<(u64, u64)>>,
}

impl Reference {
    /// `spec` is a path or an http(s) URL. `<spec>.fai` must exist next to
    /// it, and `<spec>.gzi` as well when the FASTA is bgzip-compressed.
    pub fn open(spec: &str) -> Result<Self> {
        let side = |ext: &str| read_whole(&format!("{spec}.{ext}"));
        let mut index = HashMap::new();
        for line in String::from_utf8(side("fai")?)?.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                bail!("malformed .fai line: {line}");
            }
            index.insert(
                f[0].to_string(),
                FaiEntry {
                    length: f[1].parse()?,
                    offset: f[2].parse()?,
                    line_bases: f[3].parse()?,
                    line_width: f[4].parse()?,
                },
            );
        }
        let blocks = if spec.ends_with(".gz") || spec.ends_with(".bgz") {
            let raw = side("gzi")?;
            let word = |i: usize| u64::from_le_bytes(raw[i..i + 8].try_into().unwrap());
            if raw.len() < 8 || raw.len() as u64 != 8 + 16 * word(0) {
                bail!("{spec}.gzi is not a valid bgzip index");
            }
            let mut blocks = vec![(0u64, 0u64)];
            blocks.extend((0..word(0) as usize).map(|i| (word(8 + 16 * i), word(16 + 16 * i))));
            Some(blocks)
        } else {
            None
        };
        let source = Source::open(spec)?;
        Ok(Reference { source, index, blocks })
    }

    /// The contig name as the reference spells it (`14` vs `chr14`).
    pub fn resolve_contig(&self, name: &str) -> Result<String> {
        let alternative = match name.strip_prefix("chr") {
            Some(bare) => bare.to_string(),
            None => format!("chr{name}"),
        };
        [name.to_string(), alternative]
            .into_iter()
            .find(|c| self.index.contains_key(c))
            .ok_or_else(|| anyhow!("contig '{name}' is not in the reference index"))
    }

    pub fn contig_length(&self, contig: &str) -> Result<u64> {
        Ok(self.index.get(contig).ok_or_else(|| anyhow!("unknown contig '{contig}'"))?.length)
    }

    /// Upper-cased sequence for a 1-based inclusive interval.
    pub fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let e = self.index.get(contig).ok_or_else(|| anyhow!("unknown contig '{contig}'"))?;
        if start < 1 || end > e.length || start > end {
            bail!("{contig}:{start}-{end} is outside the contig (length {})", e.length);
        }
        let byte = |p: u64| e.offset + (p / e.line_bases) * e.line_width + p % e.line_bases;
        let (from, to) = (byte(start - 1), byte(end - 1));
        let mut buf = match &self.blocks {
            None => self.source.read(from, Some(to))?,
            Some(blocks) => {
                let first = blocks.partition_point(|b| b.1 <= from) - 1;
                let last = blocks.partition_point(|b| b.1 <= to) - 1;
                let compressed_end = blocks.get(last + 1).map(|b| b.0 - 1);
                let compressed = self.source.read(blocks[first].0, compressed_end)?;
                let mut plain = Vec::new();
                MultiGzDecoder::new(&compressed[..]).read_to_end(&mut plain)?;
                let skip = (from - blocks[first].1) as usize;
                let len = (to - from + 1) as usize;
                if plain.len() < skip + len {
                    bail!("bgzip index and data disagree for {contig}:{start}-{end}");
                }
                plain[skip..skip + len].to_vec()
            }
        };
        buf.retain(|b| !b.is_ascii_whitespace());
        buf.make_ascii_uppercase();
        Ok(buf)
    }
}

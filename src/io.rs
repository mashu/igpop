//! Input parsing (coordinates table, allele FASTA) and random access to a
//! reference FASTA through its `.fai` index. The reference may be a local
//! file or an http(s) URL, plain or bgzip-compressed (`.gz` with `.gzi`);
//! only the requested bytes are read or downloaded.

use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

pub const DEFAULT_GRCH38: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa";
pub const DEFAULT_GRCH37: &str =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz";

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

/// Tab-separated: gene, build, contig, start, end, strand (+, - or . for
/// "try both"). `#` lines and a header line starting with `gene` are skipped.
pub fn read_coords(path: &Path) -> Result<Vec<Coord>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("opening {}", path.display()))?;
    parse_coords(&text, &path.display().to_string())
}

/// Coordinates from text; `label` names the source in error messages.
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

pub fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let reader = BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1
                .extend(line.bytes().filter(|b| b.is_ascii_alphabetic()).map(|b| b.to_ascii_uppercase()));
        }
    }
    Ok(out)
}

struct FaiEntry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

pub enum Source {
    Local(File),
    Remote { agent: ureq::Agent, url: String },
}

impl Source {
    /// A local path or an http(s) URL.
    pub fn open(spec: &str) -> Result<Self> {
        if is_url(spec) {
            Ok(Source::Remote {
                agent: agent(),
                url: spec.to_string(),
            })
        } else {
            Ok(Source::Local(File::open(spec).with_context(|| format!("opening {spec}"))?))
        }
    }

    /// Bytes `from..=to` of the underlying file; `to = None` reads to the
    /// end. A range that runs past the end of the file is returned short.
    pub fn read(&mut self, from: u64, to: Option<u64>) -> Result<Vec<u8>> {
        match self {
            Source::Local(file) => {
                file.seek(SeekFrom::Start(from))?;
                let mut buf = Vec::new();
                match to {
                    Some(to) => {
                        file.take(to - from + 1).read_to_end(&mut buf)?;
                    }
                    None => {
                        file.read_to_end(&mut buf)?;
                    }
                }
                Ok(buf)
            }
            Source::Remote { agent, url } => {
                let started = std::time::Instant::now();
                let range = match to {
                    Some(to) => format!("bytes={from}-{to}"),
                    None => format!("bytes={from}-"),
                };
                let response = match call_with_retry(|| agent.get(url).set("Range", &range)) {
                    Ok(r) => r,
                    Err(e) if matches!(*e, ureq::Error::Status(416, _)) => return Ok(Vec::new()),
                    Err(e) => return Err(*e).with_context(|| format!("GET {url}")),
                };
                if response.status() != 206 {
                    bail!("{url}: server ignored the byte-range request (HTTP {})", response.status());
                }
                let mut buf = Vec::new();
                response.into_reader().read_to_end(&mut buf)?;
                account(url, buf.len(), started);
                Ok(buf)
            }
        }
    }
}

static REQUESTS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static BYTES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// (number of HTTP requests, bytes downloaded) so far in this run.
pub fn download_totals() -> (u64, u64) {
    use std::sync::atomic::Ordering::Relaxed;
    (REQUESTS.load(Relaxed), BYTES.load(Relaxed))
}

fn account(url: &str, bytes: usize, started: std::time::Instant) {
    use std::sync::atomic::Ordering::Relaxed;
    REQUESTS.fetch_add(1, Relaxed);
    BYTES.fetch_add(bytes as u64, Relaxed);
    if std::env::var_os("IGPOP_TRACE").is_some() {
        let name = url.rsplit('/').next().unwrap_or(url);
        eprintln!("  GET {name}: {bytes} bytes in {:.1}s", started.elapsed().as_secs_f64());
    }
}

pub fn is_url(spec: &str) -> bool {
    spec.starts_with("http://") || spec.starts_with("https://")
}

/// One shared agent, so connections to a host are reused across requests.
fn agent() -> ureq::Agent {
    static AGENT: std::sync::OnceLock<ureq::Agent> = std::sync::OnceLock::new();
    AGENT
        .get_or_init(|| ureq::AgentBuilder::new().try_proxy_from_env(true).build())
        .clone()
}

/// A whole (small) file, e.g. an index: local path or URL.
pub fn read_whole(spec: &str) -> Result<Vec<u8>> {
    if is_url(spec) {
        http_get(&agent(), spec)
    } else {
        std::fs::read(spec).with_context(|| format!("reading {spec}"))
    }
}

/// Send a request, retrying on server overload (429, 5xx) and transport
/// errors with exponential backoff: 1, 2, 4, 8, 16 s.
fn call_with_retry(request: impl Fn() -> ureq::Request) -> std::result::Result<ureq::Response, Box<ureq::Error>> {
    const ATTEMPTS: u32 = 6;
    let mut attempt = 0;
    loop {
        let error = match request().call() {
            Ok(response) => return Ok(response),
            Err(e) => e,
        };
        let transient = match &error {
            ureq::Error::Status(code, _) => *code == 429 || (500..600).contains(code),
            ureq::Error::Transport(_) => true,
        };
        attempt += 1;
        if !transient || attempt == ATTEMPTS {
            return Err(Box::new(error));
        }
        let wait = 1u64 << (attempt - 1);
        if std::env::var_os("IGPOP_TRACE").is_some() {
            eprintln!("  retry {attempt} in {wait}s: {error}");
        }
        std::thread::sleep(std::time::Duration::from_secs(wait));
    }
}

fn http_get(agent: &ureq::Agent, url: &str) -> Result<Vec<u8>> {
    let started = std::time::Instant::now();
    let mut buf = Vec::new();
    call_with_retry(|| agent.get(url))
        .with_context(|| format!("GET {url}"))?
        .into_reader()
        .read_to_end(&mut buf)?;
    account(url, buf.len(), started);
    Ok(buf)
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

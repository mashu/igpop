//! Region queries on a bgzip-compressed, tabix-indexed file (local or
//! remote). Only the small `.tbi` index is read whole; data is read from the
//! first block that can hold an overlapping record up to the first record
//! past the region, in modest byte-range pieces.

use crate::io::{read_whole, Source};
use anyhow::{bail, Context, Result};
use flate2::read::{GzDecoder, MultiGzDecoder};
use std::collections::HashMap;
use std::io::Read;

const PIECE: u64 = 256 * 1024;

struct RefIndex {
    bins: HashMap<u32, Vec<(u64, u64)>>,
    linear: Vec<u64>,
}

pub struct Tabix {
    spec: String,
    names: Vec<String>,
    refs: Vec<RefIndex>,
    col_seq: usize,
    col_beg: usize,
}

struct Cursor<'a> {
    data: &'a [u8],
    at: usize,
}

impl Cursor<'_> {
    fn bytes(&mut self, n: usize) -> Result<&[u8]> {
        if self.at + n > self.data.len() {
            bail!("truncated tabix index");
        }
        self.at += n;
        Ok(&self.data[self.at - n..self.at])
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.bytes(4)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.bytes(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.bytes(8)?.try_into().unwrap()))
    }
}

/// Bins of the UCSC binning scheme that may hold records overlapping the
/// 0-based half-open interval [beg, end).
fn region_bins(beg: u64, end: u64) -> Vec<u32> {
    let last = end - 1;
    let mut bins = vec![0u32];
    for (shift, offset) in [(26u32, 1u64), (23, 9), (20, 73), (17, 585), (14, 4681)] {
        bins.extend((offset + (beg >> shift)..=offset + (last >> shift)).map(|b| b as u32));
    }
    bins
}

/// The `#CHROM` line of a bgzipped VCF, read from the start of the file.
pub fn sample_line(spec: &str) -> Result<String> {
    use std::io::Read;
    let mut source = Source::open(spec)?;
    let mut packed: Vec<u8> = Vec::new();
    // Headers are small; stop after 16 MB of compressed data.
    while (packed.len() as u64) < 64 * PIECE {
        let more = source.read(packed.len() as u64, Some(packed.len() as u64 + PIECE - 1))?;
        if more.is_empty() {
            break;
        }
        packed.extend_from_slice(&more);
        // Decode whatever is complete; a cut-off last block ends the read.
        let mut text = Vec::new();
        let _ = MultiGzDecoder::new(&packed[..]).read_to_end(&mut text);
        let text = String::from_utf8_lossy(&text);
        if let Some(line) = text.lines().find(|l| l.starts_with("#CHROM")) {
            // Complete only if a later line has started.
            if text[text.find("#CHROM").unwrap_or(0)..].contains('\n') {
                return Ok(line.to_string());
            }
        }
    }
    anyhow::bail!("no #CHROM header line found in {spec}")
}

impl Tabix {
    /// `spec` is the data file (path or URL); `<spec>.tbi` must exist.
    pub fn open(spec: &str) -> Result<Self> {
        let packed = read_whole(&format!("{spec}.tbi"))?;
        let mut raw = Vec::new();
        MultiGzDecoder::new(&packed[..])
            .read_to_end(&mut raw)
            .context("decompressing .tbi")?;
        let mut c = Cursor { data: &raw, at: 0 };
        if c.bytes(4)? != b"TBI\x01" {
            bail!("{spec}.tbi is not a tabix index");
        }
        let n_ref = c.i32()? as usize;
        let _format = c.i32()?;
        let col_seq = c.i32()? as usize;
        let col_beg = c.i32()? as usize;
        let (_col_end, _meta, _skip) = (c.i32()?, c.i32()?, c.i32()?);
        let l_nm = c.i32()? as usize;
        let names: Vec<String> = c
            .bytes(l_nm)?
            .split(|b| *b == 0)
            .filter(|s| !s.is_empty())
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect();
        let mut refs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let mut bins = HashMap::new();
            for _ in 0..c.i32()? {
                let bin = c.u32()?;
                let chunks = (0..c.i32()?).map(|_| Ok((c.u64()?, c.u64()?))).collect::<Result<Vec<_>>>()?;
                bins.insert(bin, chunks);
            }
            let linear = (0..c.i32()?).map(|_| c.u64()).collect::<Result<Vec<_>>>()?;
            refs.push(RefIndex { bins, linear });
        }
        Ok(Tabix {
            spec: spec.to_string(),
            names,
            refs,
            col_seq,
            col_beg,
        })
    }

    /// The file's own spelling of a contig: exact, with/without `chr`, or a
    /// RefSeq chromosome accession (NC_000014.9 for chromosome 14).
    pub fn resolve_contig(&self, name: &str) -> Option<String> {
        let bare = name.strip_prefix("chr").unwrap_or(name);
        let number = match bare {
            "X" => Some(23),
            "Y" => Some(24),
            n => n.parse::<u32>().ok(),
        };
        let accession = number.map(|n| format!("NC_{n:06}."));
        self.names
            .iter()
            .find(|n| *n == name)
            .or_else(|| self.names.iter().find(|n| *n == bare || **n == format!("chr{bare}")))
            .or_else(|| self.names.iter().find(|n| accession.as_deref().is_some_and(|a| n.starts_with(a))))
            .cloned()
    }

    /// Lines whose record starts at or before `end` and reaches `start`
    /// (by the length of column 4, VCF REF, or by END for symbolic
    /// alleles). 1-based, inclusive.
    pub fn query(&self, contig: &str, start: u64, end: u64) -> Result<Vec<String>> {
        let mut reader = self.reader(contig, start, end);
        let mut source = Source::open(&self.spec)?;
        while let Some((from, to)) = reader.request() {
            reader.feed(&source.read(from, Some(to))?)?;
        }
        Ok(reader.finish())
    }

    /// The path or URL of the indexed file.
    pub fn spec(&self) -> &str {
        &self.spec
    }

    /// A reader for one region. It does no I/O itself: it says which bytes
    /// of the file it needs next and is fed them, so the same logic serves
    /// blocking reads here and concurrent downloads in the browser.
    pub fn reader(&self, contig: &str, start: u64, end: u64) -> Reader {
        let mut reader = Reader {
            contig: contig.to_string(),
            start,
            end,
            col_seq: self.col_seq,
            col_beg: self.col_beg,
            chunks: Vec::new().into_iter(),
            current: None,
            lines: Vec::new(),
        };
        let Some(index) = self.names.iter().position(|n| n == contig).map(|tid| &self.refs[tid]) else {
            return reader;
        };
        let (beg0, end0) = (start - 1, end);
        let Some(floor) = index.linear.get((beg0 >> 14) as usize).copied() else {
            return reader;
        };
        // Chunks of every bin that overlaps the region. A record spanning
        // far (a structural variant with END) sits in a coarse bin, in a
        // chunk of its own well before the region; chunks are therefore
        // read one by one, never as one stretch from the first to the last.
        let mut chunks: Vec<(u64, u64)> = region_bins(beg0, end0)
            .iter()
            .filter_map(|b| index.bins.get(b))
            .flatten()
            .filter(|chunk| chunk.1 > floor)
            .map(|chunk| (chunk.0.max(floor), chunk.1))
            .collect();
        chunks.sort();
        // Chunks closer than one request are read together.
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for chunk in chunks {
            match merged.last_mut() {
                Some(last) if (chunk.0 >> 16) <= (last.1 >> 16) + PIECE => last.1 = last.1.max(chunk.1),
                _ => merged.push(chunk),
            }
        }
        reader.chunks = merged.into_iter();
        reader.next_chunk();
        reader
    }
}

/// One chunk of the file being read: BGZF blocks from a virtual offset to
/// a virtual offset (file offset of a block << 16 | offset inside it).
struct Chunk {
    /// The block the chunk ends in, and how much of that block belongs to it.
    last_block: u64,
    last_len: usize,
    /// The bytes to ask for next.
    pending: (u64, u64),
    /// File offset of the block at the front of `packed`.
    block_at: u64,
    /// Bytes to skip in the first block.
    skip: usize,
    packed: Vec<u8>,
    text: Vec<u8>,
}

/// See [`Tabix::reader`].
pub struct Reader {
    contig: String,
    start: u64,
    end: u64,
    col_seq: usize,
    col_beg: usize,
    chunks: std::vec::IntoIter<(u64, u64)>,
    current: Option<Chunk>,
    lines: Vec<String>,
}

impl Reader {
    fn next_chunk(&mut self) {
        self.current = self.chunks.next().map(|(first, until)| {
            let (at, last_block) = (first >> 16, until >> 16);
            // First request: to the end of the chunk, capped so a dense
            // file is still read in moderate pieces.
            let to = (last_block + 0xffff).min(at + 4 * PIECE - 1).max(at);
            Chunk {
                last_block,
                last_len: (until & 0xffff) as usize,
                pending: (at, to),
                block_at: at,
                skip: (first & 0xffff) as usize,
                packed: Vec::new(),
                text: Vec::new(),
            }
        });
    }

    /// The bytes `from..=to` of the file the reader needs next; `None`
    /// when the region is read.
    pub fn request(&self) -> Option<(u64, u64)> {
        self.current.as_ref().map(|chunk| chunk.pending)
    }

    /// Give the reader the bytes it asked for (fewer, or none, at the end
    /// of the file).
    pub fn feed(&mut self, bytes: &[u8]) -> Result<()> {
        let Some(chunk) = self.current.as_mut() else {
            return Ok(());
        };
        let from = chunk.pending.0 + bytes.len() as u64;
        chunk.pending = (from, from + PIECE - 1);
        chunk.packed.extend_from_slice(bytes);
        let mut done = bytes.is_empty();
        while !done {
            if chunk.block_at > chunk.last_block || (chunk.block_at == chunk.last_block && chunk.last_len == 0) {
                done = true;
                break;
            }
            // One whole BGZF block must be buffered.
            if chunk.packed.len() < 18 {
                break;
            }
            let block_len = u16::from_le_bytes([chunk.packed[16], chunk.packed[17]]) as usize + 1;
            if chunk.packed.len() < block_len {
                break;
            }
            let mut block = Vec::new();
            GzDecoder::new(&chunk.packed[..block_len])
                .read_to_end(&mut block)
                .context("corrupt BGZF block")?;
            chunk.packed.drain(..block_len);
            if chunk.block_at == chunk.last_block {
                block.truncate(chunk.last_len);
            }
            chunk.block_at += block_len as u64;
            chunk.text.extend_from_slice(&block[chunk.skip.min(block.len())..]);
            chunk.skip = 0;

            let mut consumed = 0;
            while let Some(nl) = chunk.text[consumed..].iter().position(|b| *b == b'\n') {
                let line = std::str::from_utf8(&chunk.text[consumed..consumed + nl])?;
                consumed += nl + 1;
                if line.starts_with('#') {
                    continue;
                }
                let fields: Vec<&str> = line.splitn(9, '\t').collect();
                if fields.len() < self.col_beg.max(5) || fields[self.col_seq - 1] != self.contig {
                    done = true;
                    break;
                }
                let pos: u64 = fields[self.col_beg - 1].parse()?;
                if pos > self.end {
                    done = true;
                    break;
                }
                // Last base of the record: REF, or END for symbolic alleles.
                let stated = fields
                    .get(7)
                    .filter(|_| fields[4].contains('<'))
                    .and_then(|info| info.split(';').find_map(|x| x.strip_prefix("END=")))
                    .and_then(|x| x.parse::<u64>().ok());
                let last_base = stated.unwrap_or(pos + fields[3].len() as u64 - 1).max(pos);
                if last_base >= self.start {
                    self.lines.push(line.to_string());
                }
            }
            chunk.text.drain(..consumed);
        }
        if done {
            self.next_chunk();
        }
        Ok(())
    }

    /// The lines of the region, in file order.
    pub fn finish(mut self) -> Vec<String> {
        // Chunks of different bins can interleave.
        let position = |line: &String| line.split('\t').nth(1).and_then(|p| p.parse::<u64>().ok()).unwrap_or(0);
        self.lines.sort_by_key(position);
        self.lines.dedup();
        self.lines
    }
}

#[cfg(test)]
mod tests {
    use super::region_bins;

    #[test]
    fn bins_cover_all_levels() {
        let bins = region_bins(0, 1);
        assert_eq!(bins, vec![0, 1, 9, 73, 585, 4681]);
        assert!(region_bins(16383, 16385).contains(&4682));
    }
}

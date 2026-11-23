//! Region-wise reading of remote, tabix-indexed VCF files, with a local
//! cache so that a region is downloaded once.

use crate::tabix::Tabix;
use crate::{note, par};
use anyhow::{Context, Result};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// (contig, first position, last position), 1-based and inclusive.
pub type Region = (String, u64, u64);

/// Public servers answer bursts of requests with HTTP 503.
const MAX_CONNECTIONS: usize = 3;

/// Sorted, merged regions; neighbours closer than `gap` are read together.
pub fn merge_regions(mut regions: Vec<Region>, gap: u64) -> Vec<Region> {
    regions.sort();
    let mut merged: Vec<Region> = Vec::new();
    for (contig, start, end) in regions {
        match merged.last_mut() {
            Some(last) if last.0 == contig && start <= last.2 + gap => last.2 = last.2.max(end),
            _ => merged.push((contig, start, end)),
        }
    }
    merged
}

/// FNV-1a, to name cache files after the file they came from.
fn fnv64(text: &str) -> u64 {
    text.bytes()
        .fold(0xcbf29ce484222325, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3))
}

/// Directory of gzipped text files, one per downloaded piece; `None` = off.
#[derive(Clone, Copy)]
pub struct Cache<'a>(pub Option<&'a Path>);

impl Cache<'_> {
    fn path(&self, file: &str, piece: &str, extension: &str) -> Option<PathBuf> {
        self.0.map(|dir| dir.join(format!("{:016x}_{piece}.{extension}.gz", fnv64(file))))
    }

    fn read(path: &Path) -> Option<Vec<String>> {
        let mut text = String::new();
        MultiGzDecoder::new(File::open(path).ok()?).read_to_string(&mut text).ok()?;
        Some(text.lines().map(str::to_string).collect())
    }

    fn write(path: &Path, lines: &[String]) -> Result<()> {
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir)?;
        }
        let partial = path.with_extension("partial");
        let mut gz = GzEncoder::new(File::create(&partial)?, Compression::default());
        for line in lines {
            writeln!(gz, "{line}")?;
        }
        gz.finish()?;
        std::fs::rename(partial, path)?;
        Ok(())
    }

    /// The lines of a small text belonging to `file` (a header, a sample
    /// table): from the cache, or from `fetch` and then cached.
    pub fn text(&self, file: &str, what: &str, fetch: impl Fn() -> Result<Vec<String>>) -> Result<Vec<String>> {
        let path = self.path(file, what, "txt");
        if let Some(lines) = path.as_deref().and_then(Self::read) {
            return Ok(lines);
        }
        let lines = fetch()?;
        if let Some(path) = &path {
            Self::write(path, &lines)?;
        }
        Ok(lines)
    }
}

/// Lines of a tabix-indexed VCF for every region, as (contig, line).
/// `spec` may contain `{contig}` (one file per chromosome). Regions already
/// in the cache are not requested again.
pub fn fetch_vcf(spec: &str, regions: &[Region], cache: Cache) -> Result<Vec<(String, String)>> {
    let file_of = |contig: &str| spec.replace("{contig}", contig);
    let path_of = |(contig, start, end): &Region| cache.path(&file_of(contig), &format!("{contig}_{start}_{end}"), "vcf");
    let mut results: Vec<Option<Vec<String>>> = regions
        .iter()
        .map(|region| path_of(region).as_deref().and_then(Cache::read))
        .collect();
    let missing: Vec<usize> = (0..regions.len()).filter(|i| results[*i].is_none()).collect();
    if !missing.is_empty() {
        let mut indexes: HashMap<&String, (Tabix, String)> = HashMap::new();
        for i in &missing {
            let contig = &regions[*i].0;
            if indexes.contains_key(contig) {
                continue;
            }
            let file = file_of(contig);
            let index = Tabix::open(&file).with_context(|| format!("opening index of {file}"))?;
            match index.resolve_contig(contig) {
                Some(name) => {
                    indexes.insert(contig, (index, name));
                }
                None => note!("warning: contig {contig} is not in {file}"),
            }
        }
        let fetched = par::map(&missing, Some(MAX_CONNECTIONS), |i| -> Result<(usize, Vec<String>)> {
            let (contig, start, end) = &regions[*i];
            let Some((index, name)) = indexes.get(contig) else {
                return Ok((*i, Vec::new()));
            };
            let lines = index
                .query(name, *start, *end)
                .with_context(|| format!("reading {contig}:{start}-{end}"))?;
            if let Some(path) = path_of(&regions[*i]) {
                Cache::write(&path, &lines)?;
            }
            Ok((*i, lines))
        });
        for result in fetched {
            let (i, lines) = result?;
            results[i] = Some(lines);
        }
    }
    note!("  {} regions, {} from cache", regions.len(), regions.len() - missing.len());
    Ok(regions
        .iter()
        .zip(results)
        .flat_map(|(region, lines)| lines.unwrap_or_default().into_iter().map(move |l| (region.0.clone(), l)))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions_merge_within_the_gap() {
        let r = |a, b| ("c".to_string(), a, b);
        assert_eq!(
            merge_regions(vec![r(100, 200), r(50, 120), r(5000, 5100)], 100),
            vec![r(50, 200), r(5000, 5100)]
        );
    }
}

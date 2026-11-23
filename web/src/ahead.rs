//! Downloading ahead: the regions a run is about to read are fetched
//! several at a time and kept, so the analysis finds them without waiting.
//!
//! The reader of a region (`igpop::tabix::Reader`) asks for exactly the
//! same bytes here as it will during the analysis, because what it asks
//! for depends only on the index and on the bytes it was given. That is
//! what makes the pieces fetched here the ones the analysis looks for.
//!
//! At most three requests run against one server at a time: public data
//! servers answer bursts with errors.

use crate::{success, KEPT};
use anyhow::{anyhow, bail, Result};
use futures::stream::{self, StreamExt};
use igpop::evidence::Download;
use igpop::fetch::Region;
use igpop::note;
use igpop::tabix::Tabix;
use js_sys::Uint8Array;
use std::collections::BTreeMap;
use wasm_bindgen::{JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
use web_sys::{Request, RequestInit, RequestMode, Response, WorkerGlobalScope};

const PER_SERVER: usize = 3;

/// One remote file and the regions to read of it.
struct Remote {
    url: String,
    regions: Vec<Region>,
}

fn js(error: JsValue) -> anyhow::Error {
    anyhow!("{}", error.as_string().unwrap_or_else(|| format!("{error:?}")))
}

/// One GET without blocking: status, body and the file's date.
async fn get(url: &str, range: Option<&str>) -> Result<(u16, Vec<u8>, Option<String>)> {
    let init = RequestInit::new();
    init.set_method("GET");
    init.set_mode(RequestMode::Cors);
    let request = Request::new_with_str_and_init(url, &init).map_err(js)?;
    if let Some(range) = range {
        request.headers().set("Range", range).map_err(js)?;
    }
    let scope: WorkerGlobalScope = js_sys::global().unchecked_into();
    let response: Response = JsFuture::from(scope.fetch_with_request(&request))
        .await
        .map_err(js)?
        .dyn_into()
        .map_err(js)?;
    let modified = response.headers().get("Last-Modified").map_err(js)?;
    let body = JsFuture::from(response.array_buffer().map_err(js)?).await.map_err(js)?;
    let bytes = Uint8Array::new(&body).to_vec();
    igpop::net::record(bytes.len());
    Ok((response.status(), bytes, modified))
}

/// Learn the date of a file, once per run, by reading its first byte.
async fn date(url: &str) -> Result<()> {
    if KEPT.with(|kept| kept.borrow().versions.contains_key(url)) {
        return Ok(());
    }
    let (_, _, modified) = get(url, Some("bytes=0-0")).await?;
    KEPT.with(|kept| {
        kept.borrow_mut().versions.entry(url.to_string()).or_insert(modified);
    });
    Ok(())
}

/// Make sure a piece is kept; returns it. `range = None` is the whole file.
async fn piece(url: &str, range: Option<&str>) -> Result<Vec<u8>> {
    if let Some(bytes) = KEPT.with(|kept| kept.borrow_mut().get(url, range)) {
        return Ok(bytes);
    }
    let (status, bytes, _) = get(url, range).await?;
    match status {
        s if s == success(range) => {
            KEPT.with(|kept| kept.borrow_mut().put(url, range, &bytes));
            Ok(bytes)
        }
        // Past the end of the file.
        416 => Ok(Vec::new()),
        other => bail!("{url}: HTTP {other}"),
    }
}

/// Fetch everything the reader of one region will ask for.
async fn region(index: &Tabix, contig: &str, start: u64, end: u64) -> Result<()> {
    let mut reader = index.reader(contig, start, end);
    while let Some((from, to)) = reader.request() {
        // Spelled exactly as `igpop::net::get_range` spells it.
        let range = format!("bytes={from}-{to}");
        reader.feed(&piece(index.spec(), Some(&range)).await?)?;
    }
    Ok(())
}

/// One file with its index: the index first, then its regions.
async fn file(url: String, regions: Vec<Region>) -> Result<()> {
    let index_url = format!("{url}.tbi");
    date(&url).await?;
    date(&index_url).await?;
    piece(&index_url, None).await?;
    // Reads the index through the blocking transport, which now has it.
    let index = Tabix::open(&url)?;
    let jobs = regions.iter().filter_map(|(contig, start, end)| {
        let name = index.resolve_contig(contig)?;
        let index = &index;
        Some(async move { region(index, &name, *start, *end).await })
    });
    let results: Vec<Result<()>> = stream::iter(jobs).buffer_unordered(PER_SERVER).collect().await;
    results.into_iter().collect()
}

/// Download ahead everything in `downloads`. Failures are reported and
/// otherwise ignored: the analysis fetches what is missing itself.
pub async fn download(downloads: Vec<Download>) {
    // Per server, the files to read, each with its regions.
    let mut servers: BTreeMap<String, Vec<Remote>> = BTreeMap::new();
    for download in downloads.into_iter().filter(|d| igpop::net::is_url(&d.spec)) {
        let mut files: BTreeMap<String, Vec<Region>> = BTreeMap::new();
        for region in download.regions {
            files.entry(download.spec.replace("{contig}", &region.0)).or_default().push(region);
        }
        for (url, regions) in files {
            let server = url.split('/').nth(2).unwrap_or("").to_string();
            servers.entry(server).or_default().push(Remote { url, regions });
        }
    }
    let regions: usize = servers.values().flatten().map(|f| f.regions.len()).sum();
    if regions == 0 {
        return;
    }
    note!(
        "downloading {regions} regions from {} servers, up to {PER_SERVER} at a time per server",
        servers.len()
    );
    // Servers side by side; the files of one server one after another.
    let per_server = servers.into_values().map(|files| async move {
        for Remote { url, regions } in files {
            let name = url.rsplit('/').next().unwrap_or(&url).to_string();
            if let Err(error) = file(url, regions).await {
                note!("  {name}: could not be downloaded ahead ({error:#}); it will be read during the analysis");
            }
        }
    });
    futures::future::join_all(per_server).await;
}

//! igpop for the browser. The page runs this module in a Web Worker.
//!
//! A run has three parts. `Session::start` places the alleles. Then the
//! regions the run will need are downloaded ahead, several at a time
//! (`ahead`). Then `Session::finish` does the analysis; it reads its data
//! with blocking requests (allowed in workers), and finds what was
//! downloaded ahead already there. If downloading ahead fails, the
//! analysis still works, one request at a time.
//!
//! Downloaded pieces are kept in a store the page provides (`store.js`),
//! so later runs need little downloading. A piece is stored under the
//! file's `Last-Modified` date, which each run asks the server for once:
//! pieces of a file that has changed are simply no longer found.
//!
//! Only data hosts that allow cross-origin requests can be read from a
//! browser. The reference, 1000 Genomes and dbSNP hosts do; the gnomAD
//! bucket does not, so the web app runs without reading gnomAD directly
//! (dbSNP still relays gnomAD's overall frequencies).

mod ahead;

use anyhow::anyhow;
use igpop::{Choice, Config, Session, Text};
use js_sys::{Array, Function, Object, Reflect, Uint8Array};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use wasm_bindgen::prelude::*;
use web_sys::{XmlHttpRequest, XmlHttpRequestResponseType};

thread_local! {
    /// Where progress messages go during a run.
    static PROGRESS: RefCell<Option<Function>> = const { RefCell::new(None) };
    static KEPT: RefCell<Kept> = RefCell::new(Kept::default());
}

fn progress(message: &str) {
    PROGRESS.with(|slot| {
        if let Some(callback) = slot.borrow().as_ref() {
            let _ = callback.call1(&JsValue::NULL, &JsValue::from_str(message));
        }
    });
}

/// Downloaded pieces: the page's store, which lasts between visits, and
/// what this run learned about the files it read.
#[derive(Default)]
struct Kept {
    /// Object with `get(key) -> Uint8Array | undefined` and `put(key, bytes)`.
    store: Option<Object>,
    /// Per file read in this run: its `Last-Modified` date, or `None` if
    /// the server gave none (such a file is not stored between visits).
    versions: HashMap<String, Option<String>>,
    /// Pieces downloaded ahead that could not go into the store.
    memory: HashMap<String, Vec<u8>>,
    /// Names of the store's pieces this run wrote or has already counted.
    seen: HashSet<String>,
    /// Bytes read from pieces that were in the store before this run.
    reused: u64,
}

impl Kept {
    /// The name of a piece in the store: file, version, byte range.
    fn stored_as(&self, url: &str, range: Option<&str>) -> Option<String> {
        let version = self.versions.get(url)?.as_ref()?;
        Some(format!("{url}\n{version}\n{}", range.unwrap_or("all")))
    }

    fn in_memory_as(url: &str, range: Option<&str>) -> String {
        format!("{url}\n{}", range.unwrap_or("all"))
    }

    fn method(&self, name: &str) -> Option<(&Object, Function)> {
        let store = self.store.as_ref()?;
        Some((store, Reflect::get(store, &name.into()).ok()?.dyn_into().ok()?))
    }

    fn get(&mut self, url: &str, range: Option<&str>) -> Option<Vec<u8>> {
        if let Some(bytes) = self.memory.get(&Self::in_memory_as(url, range)) {
            return Some(bytes.clone());
        }
        let key = self.stored_as(url, range)?;
        let (store, get) = self.method("get")?;
        let found = get.call1(store, &key.as_str().into()).ok()?;
        if found.is_undefined() || found.is_null() {
            return None;
        }
        let bytes = Uint8Array::new(&found).to_vec();
        // Counted once, however often the piece is read in this run.
        if self.seen.insert(key) {
            self.reused += bytes.len() as u64;
        }
        Some(bytes)
    }

    fn put(&mut self, url: &str, range: Option<&str>, bytes: &[u8]) {
        match (self.stored_as(url, range), self.method("put")) {
            (Some(key), Some((store, put))) => {
                let _ = put.call2(store, &key.as_str().into(), &Uint8Array::from(bytes));
                self.seen.insert(key);
            }
            _ => {
                self.memory.insert(Self::in_memory_as(url, range), bytes.to_vec());
            }
        }
    }
}

/// The status a successful request has: 206 for a range, 200 for a file.
fn success(range: Option<&str>) -> u16 {
    if range.is_some() {
        206
    } else {
        200
    }
}

/// One blocking GET; `Err` is a network failure (including a blocked
/// cross-origin request), not an HTTP error status. Also returns the
/// file's `Last-Modified` date.
fn get_once(url: &str, range: Option<&str>) -> Result<(igpop::net::Response, Option<String>), JsValue> {
    let request = XmlHttpRequest::new()?;
    request.open_with_async("GET", url, false)?;
    request.set_response_type(XmlHttpRequestResponseType::Arraybuffer);
    if let Some(range) = range {
        request.set_request_header("Range", range)?;
    }
    request.send()?;
    let status = request.status()?;
    let body = request.response()?;
    let bytes = if body.is_null() || body.is_undefined() {
        Vec::new()
    } else {
        Uint8Array::new(&body).to_vec()
    };
    igpop::net::record(bytes.len());
    Ok(((status, bytes), request.get_response_header("Last-Modified")?))
}

/// The transport igpop uses in the browser: a kept piece if there is one,
/// otherwise from the server (and then kept). A failed request is tried
/// three times; servers occasionally drop a connection.
fn transport(url: &str, range: Option<&str>) -> anyhow::Result<igpop::net::Response> {
    if let Some(bytes) = KEPT.with(|kept| kept.borrow_mut().get(url, range)) {
        return Ok((success(range), bytes));
    }
    let mut last = JsValue::UNDEFINED;
    for _ in 0..3 {
        match get_once(url, range) {
            // Status 0: the browser refused or lost the request.
            Ok(((0, _), _)) => last = JsValue::from_str("no response"),
            Ok(((status, _), _)) if status == 429 || status >= 500 => last = JsValue::from_str(&format!("HTTP {status}")),
            Ok((response, modified)) => {
                KEPT.with(|kept| {
                    let mut kept = kept.borrow_mut();
                    kept.versions.entry(url.to_string()).or_insert(modified);
                    if response.0 == success(range) {
                        kept.put(url, range, &response.1);
                    }
                });
                return Ok(response);
            }
            Err(error) => last = error,
        }
    }
    let host = url.split('/').nth(2).unwrap_or(url);
    Err(anyhow!(
        "could not read from {host} ({}). Check the internet connection; if it persists, the server may be down or may not allow requests from web pages.",
        last.as_string().unwrap_or_else(|| format!("{last:?}"))
    ))
}

fn set(object: &Object, key: &str, value: &JsValue) {
    let _ = Reflect::set(object, &JsValue::from_str(key), value);
}

fn strings(items: &[String]) -> Array {
    items.iter().map(|s| JsValue::from_str(s)).collect()
}

/// Give the module the page's store of downloaded pieces: an object with
/// `get(key) -> Uint8Array | undefined` and `put(key, bytes)`, or null.
#[wasm_bindgen]
pub fn use_store(store: JsValue) {
    KEPT.with(|kept| kept.borrow_mut().store = store.dyn_into::<Object>().ok());
}

/// The version of igpop this module was built from.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// Analyse the alleles of one FASTA text.
///
/// Resolves to `{ files: { "<name>": "<content>", ... }, summary: [...],
/// warnings: [...] }`, the same files the command-line program writes.
/// `on_progress` is called with one line of text at a time.
#[wasm_bindgen]
pub async fn run(
    fasta: String,
    fasta_name: String,
    build: String,
    candidate_marker: String,
    haplotypes: bool,
    paralog_check: bool,
    on_progress: Function,
) -> Result<JsValue, JsValue> {
    igpop::net::set_transport(transport);
    igpop::log::set_sink(progress);
    PROGRESS.with(|slot| *slot.borrow_mut() = Some(on_progress));
    // Every run asks each file's server once whether the file has changed.
    KEPT.with(|kept| {
        let mut kept = kept.borrow_mut();
        kept.versions.clear();
        kept.memory.clear();
        kept.seen.clear();
        kept.reused = 0;
    });
    let (_, before) = igpop::net::totals();

    let mut config = Config::new(Text {
        name: fasta_name,
        content: fasta,
    });
    config.build = build;
    config.candidate_marker = candidate_marker;
    config.paralog_check = paralog_check;
    config.gnomad = Choice::Off;
    if !haplotypes {
        config.phased = Choice::Off;
    }
    config.prefix = "igpop".to_string();

    let failed = |e: anyhow::Error| JsValue::from_str(&format!("{e:#}"));
    let session = Session::start(&config).map_err(failed)?;
    ahead::download(session.downloads()).await;
    let result = session.finish();
    PROGRESS.with(|slot| *slot.borrow_mut() = None);
    let reused = KEPT.with(|kept| {
        let mut kept = kept.borrow_mut();
        kept.memory.clear();
        kept.reused
    });
    let mut output = result.map_err(failed)?;
    let downloaded = igpop::net::totals().1 - before;
    output.summary.retain(|line| !line.starts_with("downloaded"));
    output.summary.push(format!(
        "{:.1} MB downloaded, {:.1} MB taken from this browser's store",
        downloaded as f64 / 1e6,
        reused as f64 / 1e6
    ));

    let files = Object::new();
    for (name, content) in &output.files {
        set(&files, name, &JsValue::from_str(content));
    }
    let out = Object::new();
    set(&out, "files", &files);
    set(&out, "summary", &strings(&output.summary));
    set(&out, "warnings", &strings(&output.warnings));
    Ok(out.into())
}

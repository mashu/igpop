//! igpop for the browser. The page runs this module in a Web Worker:
//! `run` blocks until the analysis is done, reads remote data with
//! synchronous requests (allowed in workers) and reports progress through
//! a JavaScript callback.
//!
//! Downloaded pieces are kept in a store the page provides (`store.js`),
//! so later runs need little downloading. A piece is stored under the
//! file's `Last-Modified` date: the first request to each file in a run
//! always goes to the server and tells whether the file has changed, in
//! which case its old pieces are simply no longer found.
//!
//! Only data hosts that allow cross-origin requests can be read from a
//! browser. The reference, 1000 Genomes and dbSNP hosts do; the gnomAD
//! bucket does not, so the web app runs without reading gnomAD directly
//! (dbSNP still relays gnomAD's overall frequencies).

use anyhow::anyhow;
use igpop::{Choice, Config, Text};
use js_sys::{Array, Function, Object, Reflect, Uint8Array};
use std::cell::RefCell;
use std::collections::HashMap;
use wasm_bindgen::prelude::*;
use web_sys::{XmlHttpRequest, XmlHttpRequestResponseType};

thread_local! {
    /// Where progress messages go during a run.
    static PROGRESS: RefCell<Option<Function>> = const { RefCell::new(None) };
}

/// The page's store of downloaded pieces, and what this run learned about
/// the files it read.
#[derive(Default)]
struct Kept {
    /// Object with `get(key) -> Uint8Array | undefined` and `put(key, bytes)`.
    store: Option<Object>,
    /// Per file read in this run: its `Last-Modified` date, or `None` if
    /// the server gave none (such a file is not stored).
    versions: HashMap<String, Option<String>>,
    /// Requests answered from the store, and their bytes.
    hits: (u64, u64),
}

thread_local! {
    static KEPT: RefCell<Kept> = RefCell::new(Kept::default());
}

impl Kept {
    fn key(&self, url: &str, range: Option<&str>) -> Option<String> {
        let version = self.versions.get(url)?.as_ref()?;
        Some(format!("{url}\n{version}\n{}", range.unwrap_or("all")))
    }

    fn get(&mut self, url: &str, range: Option<&str>) -> Option<Vec<u8>> {
        let key = self.key(url, range)?;
        let get: Function = Reflect::get(self.store.as_ref()?, &"get".into()).ok()?.dyn_into().ok()?;
        let found = get.call1(self.store.as_ref()?, &key.into()).ok()?;
        if found.is_undefined() || found.is_null() {
            return None;
        }
        let bytes = Uint8Array::new(&found).to_vec();
        self.hits = (self.hits.0 + 1, self.hits.1 + bytes.len() as u64);
        Some(bytes)
    }

    fn put(&self, url: &str, range: Option<&str>, bytes: &[u8]) {
        let (Some(key), Some(store)) = (self.key(url, range), self.store.as_ref()) else {
            return;
        };
        if let Some(put) = Reflect::get(store, &"put".into()).ok().and_then(|f| f.dyn_into::<Function>().ok()) {
            let _ = put.call2(store, &key.into(), &Uint8Array::from(bytes));
        }
    }
}

fn progress(message: &str) {
    PROGRESS.with(|slot| {
        if let Some(callback) = slot.borrow().as_ref() {
            let _ = callback.call1(&JsValue::NULL, &JsValue::from_str(message));
        }
    });
}

/// One synchronous GET; `Err` is a network failure (including a blocked
/// cross-origin request), not an HTTP error status.
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
    Ok(((status, bytes), request.get_response_header("Last-Modified")?))
}

/// The transport igpop uses in the browser: from the store if the piece
/// is there, otherwise from the server (and then into the store). A failed
/// request is tried three times; servers occasionally drop a connection.
fn transport(url: &str, range: Option<&str>) -> anyhow::Result<igpop::net::Response> {
    // A successful response: 206 for a range, 200 for a whole file.
    let wanted = if range.is_some() { 206 } else { 200 };
    if let Some(bytes) = KEPT.with(|kept| kept.borrow_mut().get(url, range)) {
        return Ok((wanted, bytes));
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
                    if response.0 == wanted {
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
/// Returns `{ files: { "<name>": "<content>", ... }, summary: [...],
/// warnings: [...] }`, the same files the command-line program writes.
/// `on_progress` is called with one line of text at a time.
#[wasm_bindgen]
pub fn run(
    fasta: &str,
    fasta_name: &str,
    build: &str,
    candidate_marker: &str,
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
        kept.hits = (0, 0);
    });

    let mut config = Config::new(Text {
        name: fasta_name.to_string(),
        content: fasta.to_string(),
    });
    config.build = build.to_string();
    config.candidate_marker = candidate_marker.to_string();
    config.paralog_check = paralog_check;
    config.gnomad = Choice::Off;
    if !haplotypes {
        config.phased = Choice::Off;
    }
    config.prefix = "igpop".to_string();

    let result = igpop::run(&config);
    PROGRESS.with(|slot| *slot.borrow_mut() = None);
    let mut output = result.map_err(|e| JsValue::from_str(&format!("{e:#}")))?;
    let (hits, bytes) = KEPT.with(|kept| kept.borrow().hits);
    if hits > 0 {
        output.summary.push(format!(
            "{hits} requests ({:.1} MB) were answered from this browser's store",
            bytes as f64 / 1e6
        ));
    }

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

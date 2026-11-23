//! igpop for the browser. The page runs this module in a Web Worker:
//! `run` blocks until the analysis is done, reads remote data with
//! synchronous requests (allowed in workers) and reports progress through
//! a JavaScript callback.
//!
//! Only data hosts that allow cross-origin requests can be read from a
//! browser. The reference, 1000 Genomes and dbSNP hosts do; the gnomAD
//! bucket does not, so the web app runs without reading gnomAD directly
//! (dbSNP still relays gnomAD's overall frequencies).

use anyhow::anyhow;
use igpop::{Choice, Config, Text};
use js_sys::{Array, Function, Object, Reflect, Uint8Array};
use std::cell::RefCell;
use wasm_bindgen::prelude::*;
use web_sys::{XmlHttpRequest, XmlHttpRequestResponseType};

thread_local! {
    /// Where progress messages go during a run.
    static PROGRESS: RefCell<Option<Function>> = const { RefCell::new(None) };
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
fn get_once(url: &str, range: Option<&str>) -> Result<igpop::net::Response, JsValue> {
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
    Ok((status, bytes))
}

/// The transport igpop uses in the browser. A failed request is tried
/// three times; servers occasionally drop a connection.
fn transport(url: &str, range: Option<&str>) -> anyhow::Result<igpop::net::Response> {
    let mut last = JsValue::UNDEFINED;
    for _ in 0..3 {
        match get_once(url, range) {
            // Status 0: the browser refused or lost the request.
            Ok((0, _)) => last = JsValue::from_str("no response"),
            Ok((status, _)) if status == 429 || status >= 500 => last = JsValue::from_str(&format!("HTTP {status}")),
            Ok(response) => return Ok(response),
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
    let output = result.map_err(|e| JsValue::from_str(&format!("{e:#}")))?;

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

//! Reading remote files: whole, or a byte range. The transport is
//! pluggable so the same pipeline runs as a command-line program (blocking
//! HTTP, below) and inside a browser (the embedding crate registers its own).

use anyhow::{bail, Result};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::OnceLock;

/// A response: the HTTP status and the body.
pub type Response = (u16, Vec<u8>);

/// `GET url`, optionally with a `Range` header value such as `bytes=0-99`.
/// Errors are transport failures; HTTP error statuses are returned as such.
pub type Transport = fn(url: &str, range: Option<&str>) -> Result<Response>;

static TRANSPORT: OnceLock<Transport> = OnceLock::new();
static REQUESTS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);

/// Register the transport; the first registration wins.
pub fn set_transport(transport: Transport) {
    let _ = TRANSPORT.set(transport);
}

fn request(url: &str, range: Option<&str>) -> Result<Response> {
    #[cfg(feature = "cli")]
    let transport = TRANSPORT.get_or_init(|| native::transport);
    #[cfg(not(feature = "cli"))]
    let Some(transport) = TRANSPORT.get() else {
        bail!("no transport registered (net::set_transport)");
    };
    transport(url, range)
}

/// A transport reports each response it actually received over the
/// network (and not what it answered from a store of its own).
pub fn record(bytes: usize) {
    REQUESTS.fetch_add(1, Relaxed);
    BYTES.fetch_add(bytes as u64, Relaxed);
}

/// (number of requests, bytes downloaded) so far.
pub fn totals() -> (u64, u64) {
    (REQUESTS.load(Relaxed), BYTES.load(Relaxed))
}

pub fn is_url(spec: &str) -> bool {
    spec.starts_with("http://") || spec.starts_with("https://")
}

/// A whole (small) remote file, e.g. an index.
pub fn get(url: &str) -> Result<Vec<u8>> {
    match request(url, None)? {
        (200, body) => Ok(body),
        (status, _) => bail!("GET {url}: HTTP {status}"),
    }
}

/// Bytes `from..=to` of a remote file; `to = None` reads to the end. A
/// range that starts past the end of the file gives an empty result.
pub fn get_range(url: &str, from: u64, to: Option<u64>) -> Result<Vec<u8>> {
    let range = match to {
        Some(to) => format!("bytes={from}-{to}"),
        None => format!("bytes={from}-"),
    };
    match request(url, Some(&range))? {
        (206, body) => Ok(body),
        (416, _) => Ok(Vec::new()),
        (200, _) => bail!("{url}: the server ignored the byte-range request"),
        (status, _) => bail!("GET {url}: HTTP {status}"),
    }
}

#[cfg(feature = "cli")]
mod native {
    use super::Response;
    use anyhow::{Context, Result};
    use std::io::Read;

    /// One shared agent, so connections to a host are reused.
    fn agent() -> ureq::Agent {
        static AGENT: std::sync::OnceLock<ureq::Agent> = std::sync::OnceLock::new();
        AGENT
            .get_or_init(|| ureq::AgentBuilder::new().try_proxy_from_env(true).build())
            .clone()
    }

    /// Blocking HTTP. Server overload (429, 5xx) and dropped connections
    /// are retried with exponential backoff: 1, 2, 4, 8, 16 s. Set
    /// `IGPOP_TRACE` to print every request.
    pub fn transport(url: &str, range: Option<&str>) -> Result<Response> {
        const ATTEMPTS: u32 = 6;
        let trace = std::env::var_os("IGPOP_TRACE").is_some();
        let started = std::time::Instant::now();
        let mut attempt = 0;
        loop {
            let mut call = agent().get(url);
            if let Some(range) = range {
                call = call.set("Range", range);
            }
            let (transient, error) = match call.call() {
                Ok(response) => {
                    let status = response.status();
                    let mut body = Vec::new();
                    response
                        .into_reader()
                        .read_to_end(&mut body)
                        .with_context(|| format!("GET {url}"))?;
                    if trace {
                        let name = url.rsplit('/').next().unwrap_or(url);
                        eprintln!("  GET {name}: {} bytes in {:.1}s", body.len(), started.elapsed().as_secs_f64());
                    }
                    super::record(body.len());
                    return Ok((status, body));
                }
                Err(ureq::Error::Status(code, _)) if code != 429 && !(500..600).contains(&code) => return Ok((code, Vec::new())),
                Err(e @ ureq::Error::Status(..)) => (true, e),
                Err(e @ ureq::Error::Transport(_)) => (true, e),
            };
            attempt += 1;
            if !transient || attempt == ATTEMPTS {
                return Err(error).with_context(|| format!("GET {url}"));
            }
            let wait = 1u64 << (attempt - 1);
            if trace {
                eprintln!("  retry {attempt} in {wait}s: {error}");
            }
            std::thread::sleep(std::time::Duration::from_secs(wait));
        }
    }
}

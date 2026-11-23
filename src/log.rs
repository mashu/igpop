//! Progress messages. The command-line program prints them to stderr; an
//! embedding program (the web app) registers its own sink.

use std::sync::OnceLock;

static SINK: OnceLock<fn(&str)> = OnceLock::new();

/// Register where messages go; the first registration wins.
pub fn set_sink(sink: fn(&str)) {
    let _ = SINK.set(sink);
}

pub fn emit(message: &str) {
    match SINK.get() {
        Some(sink) => sink(message),
        None => eprintln!("{message}"),
    }
}

/// `note!("...", args)`: one progress message.
#[macro_export]
macro_rules! note {
    ($($arg:tt)*) => { $crate::log::emit(&format!($($arg)*)) };
}

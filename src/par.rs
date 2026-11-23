//! `map` over a slice: on several threads in the command-line program, one
//! after another where threads do not exist (the browser).

/// Results in input order. `threads`: `None` = all cores.
#[cfg(feature = "cli")]
pub fn map<T: Sync, R: Send>(items: &[T], threads: Option<usize>, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    use rayon::prelude::*;
    match threads.and_then(|n| rayon::ThreadPoolBuilder::new().num_threads(n).build().ok()) {
        Some(pool) => pool.install(|| items.par_iter().map(&f).collect()),
        None => items.par_iter().map(&f).collect(),
    }
}

#[cfg(not(feature = "cli"))]
pub fn map<T, R>(items: &[T], _threads: Option<usize>, f: impl Fn(&T) -> R) -> Vec<R> {
    items.iter().map(f).collect()
}

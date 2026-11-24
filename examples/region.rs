//! Print the records of a tabix-indexed file in one region:
//! `cargo run --release --example region -- URL contig start end`
use igpop::tabix::Tabix;

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let index = Tabix::open(&a[1])?;
    let name = index
        .resolve_contig(&a[2])
        .ok_or_else(|| anyhow::anyhow!("contig {} not in the file", a[2]))?;
    for line in index.query(&name, a[3].parse()?, a[4].parse()?)? {
        println!("{line}");
    }
    Ok(())
}

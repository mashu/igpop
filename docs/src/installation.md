# Installation

## Pre-built binary

Download the archive for your platform from the
[releases page](https://github.com/mashu/igpop/releases), unpack it, and put `igpop` on your `PATH`. The archive also contains the
`data/` directory with ready-made coordinates.

## From source

You need a Rust toolchain (<https://rustup.rs>).

```bash
git clone https://github.com/mashu/igpop.git
cd igpop
cargo build --release
# binary: target/release/igpop
```

No other software is required: samtools, bcftools and bwa are not used.

## Network access

By default igpop reads from these hosts over HTTPS:

| Host | What |
|---|---|
| `ftp.1000genomes.ebi.ac.uk` | reference FASTA, 1000 Genomes 30x calls, phased panel and its sample table |
| `storage.googleapis.com` | gnomAD v4.1 genomes |
| `ftp.ncbi.nih.gov` | dbSNP |

Every source can be replaced by a local file with the same indexes
(`--reference`, `--vcf`, `--gnomad`, `--dbsnp`, `--phased-vcf`), or
switched off (`--no-gnomad`, `--no-dbsnp`, `--no-haplotypes`,
`--no-frequencies`).

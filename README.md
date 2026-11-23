# igpop

[![CI](https://github.com/mashu/igpop/actions/workflows/ci.yml/badge.svg)](https://github.com/mashu/igpop/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-mashu.github.io%2Figpop-blue)](https://mashu.github.io/igpop/)
[![Release](https://img.shields.io/github/v/release/mashu/igpop?include_prereleases&sort=semver)](https://github.com/mashu/igpop/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

Are the nucleotide differences that define a germline allele seen in human
populations, and how often?

igpop aligns each allele (IGH out of the box; IGK, IGL and TR with your own
coordinates) to the locus of its gene, calls its differences from the
reference genome, and judges each against 1000 Genomes 30x, gnomAD and the
datasets relayed by dbSNP. Every allele gets a verdict in words, the reasons
when it does not pass, and a check whether its whole combination of variants
occurs on one chromosome in phased genomes. Only the regions needed are
downloaded.

## Quick start (human IGH)

1. Install Rust with [rustup](https://rustup.rs/) (Linux and macOS; on
   Windows run the installer from that page):

   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. Install igpop (open a new terminal first, so that `cargo` is found):

   ```bash
   cargo install --git https://github.com/mashu/igpop
   ```

3. Run it on your alleles:

   ```bash
   igpop variants --alleles alleles.fasta -o run
   ```

4. Open `run.report.md`.

Nothing else is needed for IGH: gene coordinates (GRCh38) are built in and
population data are read over the internet.

### Allele names

FASTA headers must follow the usual `GENE*ALLELE` convention, for example
`IGHV1-69*01`, `IGHD3-10*01`, `IGHJ4*02`. The text before `*` is the gene and
decides where the sequence is aligned. Names containing `_S` (as in
`IGHV1-69*01_S1234`) are treated as novel candidates to be judged; all others
are treated as trusted alleles that calibrate the score.

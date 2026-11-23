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

```bash
cargo install --git https://github.com/mashu/igpop
igpop variants --alleles alleles.fasta -o run
```

Then open `run.report.md`.

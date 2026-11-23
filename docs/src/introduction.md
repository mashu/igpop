# igpop

igpop answers one question about a germline allele (for example a candidate
IGHV allele): **are the nucleotide differences that define it seen in human
populations, and how often?**

For every allele it

1. aligns the sequence to the locus *you* specify on the reference genome,
2. writes each difference from the reference as a normalised variant,
3. looks each variant up in the 1000 Genomes 30x callset and in gnomAD,
4. reports per-population frequencies with confidence bounds, a frequency
   class, and the quality of the site the variant sits on,
5. scores each variant against trusted alleles, calibrating itself, and
   gives every allele a verdict in words: `supported`, `inconclusive` or
   `needs_confirmation`,
6. writes a report that lists the alleles to look at and checks the run.

New here? Go to the [Quick start](quickstart.md), or try the
[web app](app/index.html): nothing to install.

The locus is supplied by the user because short-read aligners cannot place
sequences reliably in immunoglobulin and T-cell receptor loci: the genes are
too similar to each other.

igpop never downloads whole genome or VCF files. It reads only the byte
ranges that cover the requested loci, using the `.fai`, `.gzi` and `.tbi`
indexes that sit next to the remote files.

It is the successor of [IgSNPer.jl](https://gitlab.com/gkhlab/igsnper.jl),
which took the first dbSNP sub-dataset above a fixed frequency threshold.
igpop replaces that with allele counts, confidence bounds and explicit site
quality, and keeps IgSNPer's column names where the meaning is the same.

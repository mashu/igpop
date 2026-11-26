# igpop

igpop answers one question about a germline allele (for example a candidate
IGHV allele): **are the nucleotide differences that define it seen in human
populations, and how often?**

For every allele it

1. aligns the sequence to the gene *you* name, at coordinates you can
   see and change,
2. writes each difference from the reference as a normalised variant,
3. looks each variant up in seven collections of people (1000 Genomes,
   gnomAD, ALFA, GenomeAsia, ToMMo, the Westlake BioBank and the studies
   collected by dbSNP), group by group,
4. scores each variant by how much more typical its frequency is of real
   variants than of arbitrary base changes, measured in the same data,
5. looks for the whole sequence on the phased chromosomes of 1000 Genomes,
6. gives every allele a verdict in words (`supported`, `inconclusive` or
   `needs_confirmation`) with the reason and the numbers behind it, and
   writes a report that lists the alleles to look at and checks the run.

What happens in which order, and why, is on the page
[Method, step by step](method.md).

New here? Go to the [Quick start](quickstart.md), or try the
[web app](app/index.html): nothing to install.

The locus is supplied by the user because short-read aligners cannot place
sequences reliably in immunoglobulin and T-cell receptor loci: the genes are
too similar to each other.

igpop never downloads whole genome or VCF files. It reads only the byte
ranges that cover the requested loci, using the `.fai`, `.gzi` and `.tbi`
indexes that sit next to the remote files.

It is the successor of [IgSNPer.jl](https://gitlab.com/gkhlab/igsnper.jl),
which took the first dbSNP sub-dataset above a fixed frequency threshold
and counted the differences below it. igpop replaces that with allele
counts per population, confidence bounds, a calibrated score and a check
of the whole sequence, and reports site quality next to every verdict.

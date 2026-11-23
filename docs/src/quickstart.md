# Quick start

You need one file: your alleles in FASTA format, named like `IGHV1-69*01`
(the gene is the text before `*`).

```bash
igpop variants --alleles my_alleles.fasta -o my_run
```

That is the whole command for human IGH alleles. igpop ships the IGH gene
coordinates, reads the reference genome, 1000 Genomes, gnomAD and dbSNP over
the internet (only the few regions it needs, roughly 150 MB for a full IGHV set the
first time, almost nothing on later runs), and calibrates its score by itself.

When it finishes, the last lines on screen look like this:

```
566 alleles: 445 ok, 68 ambiguous_locus, 0 low_identity, 53 no_coordinates
verdicts: 371 supported, 87 inconclusive, 12 needs_confirmation, 43 matches_reference, 53 not_placed
read my_run.report.md first
```

## What to open, in this order

1. **`my_run.report.md`**: the answer in words. Which alleles are fine,
   which need a second look and why, and whether the run itself is healthy.
   It is plain text (Markdown); any text editor shows it, and GitHub, VS Code
   or RStudio render the tables.
2. **`my_run.final.tsv`**: the alleles that passed.
   **`my_run.alleles.tsv`**: all alleles, as a table for Excel or R. The
   first columns are `verdict`, `filters` (why it did not pass), `reason`,
   `notes` and `hap_status` (was the whole combination of variants found
   on one chromosome?).
3. **`my_run.variants.tsv`**: only when you want the numbers behind one
   variant: frequencies in every population, site quality, rsID.
4. **`my_run.alignments.txt`**: only when a placement looks odd: the
   alignments of all alleles of a gene, stacked.

[Reading the results](reading-results.md) explains the verdicts and what to
do with each.

## Which alleles are judged?

igpop needs examples of real variation to calibrate itself. By default,
alleles whose name contains `_S` (the suffix IgDiscover-style tools give to
novel candidates, `IGHV1-18*01_S2706`) are **candidates**; all other alleles
are **trusted**. Both get a verdict; only trusted ones shape the score.

If your names follow another convention, either change the marker
(`--candidate-marker _novel`) or list the trusted names in a file, one per
line (`--trusted trusted.txt`).

## Genes that exist only on the older genome build

A few IGHV genes are absent from the GRCh38 reference. They are reported as
`not_placed`, with the reason "its gene has coordinates on GRCh37 only". Run
them with a second command:

```bash
igpop variants --alleles my_alleles.fasta --build 37 -o my_run_37
```

On that build the evidence is limited to what dbSNP relays, and there is
no haplotype check ([Limitations](limitations.md)).

## Other loci (IGK, IGL, TCR) or other coordinates

Supply your own coordinates table ([Input files](inputs.md)):

```bash
igpop variants --alleles trb.fasta --coords trb_coordinates.tsv -o trb
```

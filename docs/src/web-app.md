# The web app

<https://mashu.github.io/igpop/app/> runs igpop in the browser. Drop a FASTA
of human IGH alleles on the page; nothing has to be installed.

## What it does

The same program as the command line, compiled to WebAssembly. It runs
inside your browser tab:

- **Your sequences stay on your computer.** Nothing is uploaded. The page
  only downloads public data: the regions of the reference genome, of 1000
  Genomes and of dbSNP that cover your genes.
- **The same result files.** The *Files* tab offers exactly the tables and
  the report the command line writes ([Output files](outputs.md)).

## How results are shown

| Part of the page | Shows |
|---|---|
| Tiles | Number of alleles per verdict. Click one to filter the list. |
| Allele list | One line per allele, grouped by gene: verdict, the reasons it did not pass in plain words, the haplotype check. Search and filters above. |
| An opened allele | Every difference from the reference with its verdict, score and a frequency line: one mark per population, placed by the share of chromosomes that carry the base, from one in a million to all, with the 1% line. The `?` beside a column title explains the column. A difference at a position where the 1000 Genomes genotypes are less dependable carries the label "less reliable position"; pointing at it says which check failed. Below: the haplotype check as bars per superpopulation, and the alignment of the gene. In the alignment, pointing at a column highlights it and shows the position, the reference base, the base of the allele and, for a difference, its rsID, verdict and frequency. With the keyboard, the arrow keys jump between the differences of the allele. |
| Score and checks | The calibration table as bars, and how the trusted alleles fared. |
| Report | The run report. |

The page can also display results of the command-line program: drop the
`.alleles.tsv` and `.variants.tsv` files (and any of the others) on it.

## Differences from the command line

- **No gnomAD ancestry groups.** Browsers may only read from servers that
  permit it. The 1000 Genomes and dbSNP servers do; the gnomAD server does
  not. gnomAD's overall frequencies still arrive through dbSNP, but the
  nine ancestry groups do not, so support is a little lower for variants
  that are frequent in one gnomAD group only. In the IGHV example, 15 of
  566 alleles get a stricter verdict in the web app (14 supported become
  inconclusive, 1 inconclusive becomes needs_confirmation); none gets a
  more lenient one. For final results use the command line.
- **Slower the first time.** One processor core, and downloads limited
  to three at a time per data server (the servers refuse bursts). The
  downloads are most of the waiting. Later runs on the same genes use the
  stored data. D and J sets
  take well under a minute; a full V set downloads roughly 100 MB and takes
  several minutes.
- **Fewer options.** Genome build, candidate marker, haplotype check on or
  off. Own coordinates, trusted lists and the other options need the
  command line.

## Downloaded data is kept in the browser

The pieces of the public data files that a run downloads are stored in the
browser (its private file storage for this site), so the next run on the
same genes downloads almost nothing, also after closing the browser. Only
public data is stored, never your sequences or results.

- The start page shows how much is stored, with a link to delete it.
- At each run the page asks every data file's server once whether the
  file has changed (its date). Pieces of a changed file are not reused,
  so a new dbSNP release is picked up by itself.
- The store starts afresh when it passes 1.5 GB or six months.
- A browser may still remove stored data when the disk runs short; the
  page asks it not to. Private windows keep nothing.
- With the page open in two tabs, only the first uses the store.

## Requirements

A current Chrome, Edge, Firefox or Safari, and a connection that reaches
`ftp.1000genomes.ebi.ac.uk` and `ftp.ncbi.nih.gov`.

## While it runs

The page shows four steps: placing the alleles on their genes, fetching
the population data (with a bar that counts the regions), reading
frequencies and haplotypes, and judging the alleles. "Show every step"
lists each file as it is read.

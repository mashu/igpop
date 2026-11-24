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

## Gene names are checked before the run

An allele is matched to a gene by the part of its name before `*`, and the
gene must be in the gene list, which says where each gene lies on the
reference genome. After you choose a file, and before anything is
downloaded, the page shows how the names of the file connect to the genes
of the list, one line per connection:

| Arrow | Meaning |
|---|---|
| solid, "same name" | The gene list has a gene of this name. These lines are folded away. |
| solid, "built-in rule" | The list knows these alleles under another gene name. What the databases call IGHV4-4 is two genes in the list: `IGHV4-4*01`, `*02`, `*03` and candidates named after them belong to `IGHV4-4B`, and `*07`, `*08`, `*09` to `IGHV4-4A`. |
| solid, "your match" | You said which gene these alleles belong to. |
| dashed, "by sequence" | Each allele will be placed on the gene its sequence fits best; the results say which. |
| broken, "not matched" | No gene of the list has this name. Without a choice the alleles are reported as `not_placed`. |
| broken, "left out" | You took these alleles out of the analysis. |

Above the lines, counters say how many alleles are in each state, and the
run button repeats how many are still not matched.

Both sides can be changed. On the side of the file, the menu of a line
sets what to do with its alleles, and "edit names" opens the alleles one
by one: correct a name (your file is not changed; the corrected name is
used for this run) or decide for a single allele. On the side of the list,
**Gene list** lets you look through the genes, remove one, add one (gene,
build, contig, first and last base, strand), load your own list from a
tab-separated file, or save the current one. Your list and your matches
are kept in this browser; "Back to the built-in list" discards the list.

## Where the evidence comes from

The same step shows, before the run, which data will be read for which
genes: one card per genome build with its sources (frequencies or whole
sequences, and how many people), a card for genes that lie outside the
chromosomes and get no population data, and for every gene of the file the
sources that will be used for it.

After the run, the tab **Data used** shows what each source had to say:
how many of the distinct differences of the run it has seen, for how many
it showed a frequency of 1% or more, how many alleles the phased panel
could look up completely and how many it found. The studies that dbSNP
passes on are listed underneath, with their estimated sizes.

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
| An opened allele | Every difference from the reference with its verdict, score and a frequency strip: one mark per population on a scale from one in a million to 100%, with the 1% line. Below: the haplotype check as bars per superpopulation, and the alignment of the gene with this allele highlighted. |
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
- **Slower.** One download at a time and one processor core; the
  downloads are most of the waiting. D and J sets
  take well under a minute; a full V set downloads roughly 100 MB and takes
  several minutes.
- **No cache.** Each run downloads its data again.
- **Fewer options.** Genome build, candidate marker, haplotype check on or
  off. Own coordinates, trusted lists and the other options need the
  command line.

## Requirements

A current Chrome, Edge, Firefox or Safari, and a connection that reaches
`ftp.1000genomes.ebi.ac.uk` and `ftp.ncbi.nih.gov`.

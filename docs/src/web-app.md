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

## While it runs

The page shows four steps. During the second, "Fetching population data",
each source has its own line: how many of its regions are ready, how many
were downloaded now and how many were already in this browser from an
earlier run. A source that is new to the browser (gnomAD after an update,
or the data of genes not analysed before) downloads in full while the
others are taken from the browser. "Steps so far" lists what has been
done.

## Data kept in the browser

Downloaded pieces of the population files are kept in the browser, so a
later run on the same genes needs little or no downloading. They belong
to the site, not to a version of the app: an update does not remove them.
The start page says how much is kept and offers to delete it.

Several tabs can use the kept data at once. One tab adds to them; a
second tab reads them, sees what the first one adds, and keeps its own
downloads only until it is closed (the page says so). Data are downloaded
again when the genes are different (neighbouring genes are fetched as one
region, so another set of genes gives other regions), when the kept data
are older than 180 days or larger than 1.5 GB, or when the browser's site
data are cleared.

## Differences from the command line

- **gnomAD comes from its copy at Amazon.** Browsers may only read from
  servers that permit it. gnomAD's server at Google does not, but the
  project's copy at Amazon holds the same files and does, so the web app
  reads it there and gives the same evidence as the command line. gnomAD
  is the largest download (about 60 MB for a full V set); it can be
  switched off under Options, in which case its overall frequencies still
  arrive through dbSNP but the ancestry groups do not, and a few alleles
  get a stricter verdict.
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

## An opened allele

An opened allele has four parts, each with a heading and the question it
answers.

1. **Verdict.** The verdict and the reason for it in one sentence.
2. **Differences from the reference genome.** Is each difference, taken
   by itself, carried by people? One row per difference.
3. **Haplotype check.** Does the whole sequence, with all its differences
   together, exist on single chromosomes of real people? The bars say
   what they count in their title. When every difference is in the phased
   panel they count chromosomes with exactly this sequence. When only some
   are, they count chromosomes that carry those differences, which is not
   the frequency of the allele, and the panel says so.
4. **Alignment.** Where the differences lie and what the other alleles of
   the gene have there.

## The evidence behind one difference

Pointing at a difference in the table (its description or its frequency
line) shows a short summary: how many chromosomes carry it in each source,
the three populations where it is most frequent, the guaranteed frequency
and the score. The full evidence opens with "show".

In an opened allele, "show" in the Details column of a difference (or a click
on it in the alignment) shows:

- **Where**: its position in your sequence and on the genome, and the two
  bases.
- **Who carries it**: every population that looked at the position, with
  the number of chromosomes that carry the base out of those examined: the
  five continental groups of 1000 Genomes, the ancestry groups of gnomAD, the populations of ALFA, GenomeAsia, ToMMo (Japan) and WBBC (China), and the studies dbSNP passes on. Populations without
  a single copy are listed too. For dbSNP studies the counts are worked
  out from the estimated study size and marked "about".
- **Why this verdict and score**: which population gives the guaranteed
  frequency (the support), the score that support has on the curve of the
  run and how well that score is known, how many of all possible changes
  are supported at least as well, and whether a test showed a frequency of
  1% or more.
- **What it was compared with**: a picture on one axis of support. The
  line is the score, with a band for its uncertainty and dashed lines at
  the two limits (2 and 0.1 unless you move the first). Below it, every real variant of the run is
  a dot (point at it for its rsID and gene) and all possible changes are
  bars. The difference is marked on the line. Because the line has no
  steps, a difference just below a limit is visibly close to one just
  above it. See [The score](score.md).

## The score first, the limit as a setting

The list of alleles shows search, gene and order; the remaining filters
(trusted or candidate, haplotype, placement) are under **More filters**,
and the limit described below under **Verdict settings**. Both fold away and
say in their title what is set.

Every allele shows its weakest score as a number with a dot on a
continuous colour scale (from the colour of "needs confirmation" at 0.01
through grey at 1 to the colour of "supported" at 100), and **Order** sorts
the list by it. An allele is only as well supported by frequency as its
weakest difference.

A yes/no verdict needs a cut, and no cut is given by the data. The field
"Count a difference as supported from a score of ..." sets it, in steps of
0.1 from 0 to 100; the slider under it covers 0 to 20. The analysis runs
with 2. Under the slider the app says what the chosen limit means:

- the support at which the score reaches it (0.28% for 2 in the IGHV
  example, 1.7% for 10),
- the share of all possible single-base changes that score this high,
  which is how often an arbitrary change would pass,
- how many alleles are supported, and how many more or fewer than at 2,
- how such a score is usually read: 1 to 3 weak, 3 to 10 positive, 10 or
  more strong evidence. At 1 or below a difference passes although its
  support is no more typical of real variants than of arbitrary changes,
  and at 0 everything with a score passes.

Moving the slider does not run the analysis again: scores, classes and
haplotype checks stay as measured, and only their reading changes, exactly
as `--supported-score` would change it on the command line. The verdict
labels, the counts at the top and the tables under **Files**
(`final.tsv`, `alleles.tsv`, `variants.tsv`) follow the limit; the written
report describes the run at 2. The limit for "needs confirmation" (0.1)
is not moved.

The site files give counts for continental groups and ancestry groups, not
for the single populations within them.

## The alignment of a gene

"Alignment of this gene" in an opened allele shows every allele of the
gene under the reference, 5' to 3'.

- **Positions** are counted from the first base of the gene; the few bases
  shown before and after it are not numbered.
- **Differences that are not supported are marked**, in every allele of
  the gene: underlined when the difference is inconclusive or has no
  population data, framed when it needs confirmation. A supported
  difference is a plain coloured letter. The marks follow the limit set
  under Verdict settings.
- **Regions.** For heavy-chain V genes a line shows FR1, CDR1, FR2, CDR2,
  FR3 and the start of CDR3, and the popup of a column names its region.
  Pointing at a region highlights its columns in all alleles; a click
  keeps the highlight until the region is clicked again.
  CDR3 is only a few bases long here: a V gene ends just after the second
  cysteine, and the rest of CDR3 comes from the D and J genes.
  The borders are placed from the conserved residues of a V gene (the two
  cysteines and the tryptophan) and standard framework lengths, not from
  an IMGT alignment, so they can be off by a codon or two. For the 53
  IGHV genes of the built-in list the resulting CDR1 and CDR2 lengths are
  the usual ones (8, 9 or 10 and 7, 8, 9 or 10 codons).
- Pointing at a column shows what is there; a click on a difference of the
  opened allele shows its evidence.
- The header gives the gene, its place on the genome and a key to the
  marks; the "?" beside it holds the full explanation. The dot before
  each allele's name is the verdict of the whole allele.

## Alleles that were not placed cleanly

A sequence is "not placed cleanly" when it fits another gene of the list
at least as well as the gene of its name (`ambiguous_locus`), is too
different from its gene (`low_identity`), or names a gene the list does
not have. Such alleles carry a label in the allele list ("Fits another
gene too"), the **Placement** filter shows only them, and the warning
about them has a link that sets this filter. Opening one explains why,
with the gene it fits and the two alignment scores. An allele that fits
another gene is still judged, but its differences, and the frequencies
found for them, may belong to the other copy of the gene.

## Trying it without a file

"Try the example" on the first page loads 566 human IGHV alleles of the
KIARVA set (174 database alleles, 392 candidates) as if you had chosen the
file. The first run downloads about 100 MB.

# Questions and answers

**My candidate has an rsID. Is it real?**
Not for that reason. dbSNP assigns an rsID to any submitted change,
including errors. Look at the verdict and `support`.

**My candidate is `needs_confirmation`. Is it an error?**
Not necessarily. It means population data do not back it. Real alleles
restricted to a family or to a population that is poorly represented in
1000 Genomes and gnomAD look the same. Confirm it another way.

**My candidate is `inconclusive`. What now?**
Population data cannot decide. Use what you know from your own data: number
of unrelated carriers, expression level, linkage to known alleles, genomic
validation.

**A trusted database allele is not `supported`.**
Either one of its bases is doubtful, or it is real but rare, or its
variants are written differently in the population data (a cluster of
neighbouring changes can be recorded as one insertion or deletion there;
see [Limitations](limitations.md)). Check its variants in `variants.tsv`
and the alignment.

**Why is `supported` not the same as "all variants common"?**
A variant is `supported` when it is classed `common`, or when its
guaranteed frequency is at a level where variants of real alleles are at
least ten times more frequent than arbitrary changes. And a whole allele is
`supported` when chromosomes with exactly its sequence exist in phased
genomes, however rare its single differences are.

**The score is above 1 but the allele is inconclusive.**
A score near 1 means the frequency says nothing either way. IGHV is so
variable that about 1 in 17 arbitrary base changes coincides with a real
variant at 0.1% or more, so a single rare difference cannot score high.
Such an allele passes only through the haplotype check, and only if all
its variants are in the phased panel.

**Many alleles say "fits another locus at least as well".**
Duplicated genes (IGHV1-69/1-69D, IGHV3-23/3-23D, IGHV3-30 and its
relatives) cannot be told apart by sequence. Short reads in the population
datasets have the same problem there, so frequencies in these genes are
less certain. The verdict stands, with that caveat.

**Half my alleles are `not_placed`.**
Read the reason. "coordinates on GRCh37 only": run again with `--build 37`.
"missing from the coordinates table": the gene name before `*` does not
match a row of the table; add it or use `--gene-aliases`.

**The run stopped with a download error.**
Servers occasionally refuse connections. Run the same command again: every
region already fetched is in `.igpop_cache` and is not downloaded twice.

**Do results change between runs?**
With the same inputs and the same dataset releases, no. dbSNP is read from
its "latest release" address, so rsIDs and COMMON flags follow dbSNP
updates. Delete `.igpop_cache` to refresh.

**Can I use it for mouse or macaque?**
No. The population datasets are human.

**How do I cite the numbers?**
State the igpop version (top of the report), the datasets (section 4), and
the calibration source (section 5).

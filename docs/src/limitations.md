# Limitations

## The evidence

- **Short reads in duplicated genes.** Every dataset igpop reads was made
  from short reads placed on a reference genome. Where a read fits two
  near-identical genes it can be placed on the wrong one, and a base that
  differs between the two genes then appears as a variant of one of them.
  All datasets share this error, so their agreement does not rule it out.
  igpop reports the symptoms ([site quality](site-qc.md)) and cannot
  recover the truth. In the IGHV example, 390 of 537 distinct differences
  lie in sequence that gnomAD flags as duplicated, and 276 of the 439
  supported alleles have at least one difference at a site that 1000
  Genomes or gnomAD flag. Most of these alleles are real; the flags say
  that the population data alone cannot prove it.
- **No long reads.** Assemblies from long reads and pangenome graphs
  resolve these genes much better. igpop does not read them yet: they
  cover hundreds of people, not tens of thousands, and have no common
  format per population.
- **Overlapping people.** The 1000 Genomes samples are part of gnomAD, and
  several studies relayed by dbSNP overlap. igpop never adds counts
  across datasets, so nobody is counted twice within a result. A variant
  "seen in three datasets" may still have been seen in the same people
  three times: the number of datasets is not a number of replications.
- **Shared loci.** Genes annotated at the same interval (`shares_locus_with`)
  cannot be told apart in population data.
- **Genes missing from the reference.** Four IGHV genes are absent from
  GRCh38 and are judged on GRCh37, with older and smaller data; one more
  (IGHV4-38-2) is on an unplaced scaffold and gets no population evidence.
  See [Datasets](datasets.md#genes-that-the-reference-genome-lacks).
- **Absent is not disproof.** A variant with no record may sit where
  genotyping fails. igpop does not yet read per-site coverage, so it cannot
  separate "confidently not there" from "not callable".
- **By frequency, the verdict is per variant.** Two common variants that
  never occur on the same chromosome are both `common`, so a chimera of
  two real alleles passes by frequency. The
  [haplotype check](haplotypes.md) can confirm an allele as a whole, but
  it cannot refute one where the phased panel lacks some of its variants,
  which in IGHV is over half of the alleles.
- **A single rare difference is the hardest case.** Its frequency cannot
  separate it from an error that happens to hit a real variant, and the
  haplotype check adds little for one difference.
- **dbSNP-relayed datasets have no counts.** Their sizes are estimated and
  per-site coverage is unknown ([Datasets](datasets.md)).
- **Superpopulations only.** The 1000 Genomes source provides counts for
  five superpopulations, not the 26 populations.
- **Representation.** A cluster of neighbouring changes can be recorded in
  the population data as one insertion or deletion. igpop compares
  normalised variants one by one and then misses the match (known case:
  IGHV4-28*02). Such alleles come out weaker than they are.
- **Verdicts are about population support.** A real private allele gets
  `needs_confirmation`; an artefact shared by all short-read datasets can
  get `supported`. igpop reports evidence, not the probability that an allele is real.

## The score

- **"Real" means "like the accepted alleles".** The curve is fitted to the
  differences of alleles that the run takes as real: by default the
  alleles without the candidate marker in their name. These are mostly
  alleles known for a long time, found in many people. The score
  therefore says how typical a support is *of the variants of accepted
  alleles*, which favours common variants: a real allele private to one
  family looks like an error to it, and should.
- **The yardstick comes from the same file.** The trusted alleles
  calibrate the curve and are then judged by it. In the example this
  circle does no harm: fitting the curve 42 times, each time without the
  alleles of one gene, and judging that gene with the reduced curve gives
  the same verdicts (99 of 104 trusted alleles pass either way; for 93 of
  them every difference is `common`, which does not involve the curve).
  An allele set that is wrong in a consistent way would still calibrate
  itself wrong.
- **All possible changes are treated alike.** The background does not
  distinguish kinds of base change. Transitions are more frequent in
  people than transversions, so an arbitrary transition is found more
  often: of artificial single-base changes of real alleles, 5.3% of the
  transitions and 1.8% of the transversions come out `supported`. An
  error that is a transition is somewhat more likely to pass.
- **One difference at a time.** The class and the score are computed per
  difference and are not corrected for the number of differences or
  alleles in a run. A run of thousands of candidates will contain some
  that pass by coincidence; the share is the one in the last column of
  the score table.
- **Best case in the example.** The alleles of the worked example (KIARVA)
  were inferred from people of the 1000 Genomes collection, the same
  people igpop then consults, though with another technology. Alleles
  inferred from other people will be supported
  less often for the same quality.
- **Insertions and deletions borrow their score.** It is read from the
  curve of single-base changes as a lower bound, and only where that bound
  is above 1; a rarely seen insertion or deletion has no score
  ([The score](score.md)).
- **Built-in calibration is human IGHV.** For other loci it is a stand-in
  until a run of that locus can calibrate itself.

# Limitations

- **Short reads in repetitive loci.** Both frequency datasets are short-read
  callsets. Where reads cannot be placed, igpop can report the symptoms
  (site filters, `segdup`) but cannot recover the truth.
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
- **Insertions and deletions borrow their score.** It is read from the
  curve of single-base changes as a lower bound, and only where that bound
  is above 1; a rarely seen insertion or deletion has no score
  ([The score](score.md)).
- **Built-in calibration is human IGHV.** For other loci it is a stand-in
  until a run of that locus can calibrate itself.

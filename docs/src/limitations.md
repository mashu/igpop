# Limitations

- **Short reads in repetitive loci.** Both frequency datasets are short-read
  callsets. Where reads cannot be placed, igpop can report the symptoms
  (site filters, `segdup`) but cannot recover the truth.
- **Shared loci.** Genes annotated at the same interval (`shares_locus_with`)
  cannot be told apart in population data.
- **Genes missing from the reference.** A gene absent from the build cannot
  be assessed on that build. Five IGHV genes in the shipped coordinates
  exist only on GRCh37.
- **Absent is not disproof.** A variant with no record may sit where
  genotyping fails. igpop does not yet read per-site coverage, so it cannot
  separate "confidently not there" from "not callable".
- **The verdict is per variant.** Two common variants that never occur on
  the same chromosome are both `common`. The
  [haplotype check](haplotypes.md) addresses this, but only where the
  phased panel contains all of an allele's variants, which in IGHV is
  under half of the alleles.
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
  get `supported`. A probability needs a prior that only you can supply.
- **Insertions and deletions have no score.** Their verdict comes from the
  class alone: `common` is `supported`, `absent` is `needs_confirmation`.
- **Built-in calibration is human IGHV.** For other loci it is a stand-in
  until a run of that locus can calibrate itself.

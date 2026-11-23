# Site quality filters

A frequency is only as good as the genotypes behind it. In immunoglobulin
and T-cell receptor loci short-read genotypes are wrong more often than
elsewhere, for three reasons:

- **Paralogs.** Genes and gene copies are nearly identical. Reads from one
  copy are placed on another, which creates variants that do not exist at
  that position or distorts the frequency of ones that do.
- **Copy-number variation.** Many genes are deleted or duplicated on some
  haplotypes. A person with one copy looks homozygous; a person with three
  looks like a heterozygote with odd proportions.
- **Cell-line artefacts.** 1000 Genomes DNA comes from B-cell lines. B cells
  delete part of the IGH locus on at least one chromosome when they
  rearrange V(D)J, so one haplotype can be missing from the data.

The filters below detect the traces these problems leave. They do not say
"the variant is false"; they say "the frequency at this site is biased, and
here is the symptom". That is why igpop keeps `class` (how frequent) and
site quality (how trustworthy) in separate columns.

## 1000 Genomes: `kgp_site_qc`

These are the four criteria the 1000 Genomes consortium applied before
building its phased panel. `PASS` means all four are met.

| Value | Criterion | What it detects |
|---|---|---|
| `FILTER=<name>` | The record's FILTER is not `PASS`. | GATK's variant quality model (VQSR) judged the read evidence untypical of true variants: poor mapping quality, strand bias, variants bunched at read ends. Typical of mismapped reads. |
| `missingness>5%` | Fewer than 95% of chromosomes have a genotype. | Many samples had too little usable coverage here, usually because the region is deleted in them or reads cannot be placed. The frequency then describes a biased subset of people. |
| `HWE` | Genotype proportions depart from Hardy–Weinberg expectation (p ≤ 1e-10) in **every** superpopulation. | Too many heterozygotes: two gene copies collapsed onto one position, so everyone appears to carry both bases. Too few: one haplotype is dropping out (deletion, cell-line artefact). |
| `mendelian_errors>5%` | In more than 5% of the trios the child's genotype cannot come from the parents. | Genotypes are unreliable at the site. In IGH the usual cause is a haplotype missing in the cell line of one family member. |

Supporting numbers are in `kgp_filter`, `kgp_vqslod`, `kgp_AN`,
`kgp_excess_het_p` and `kgp_mendelian_error_rate`.

### What to expect in IGH

Well-known, genuinely common polymorphisms often fail
`mendelian_errors>5%`, most of all in the V genes closest to the D–J
region (IGHV1-2, IGHV1-3, IGHV2-5, IGHV3-7, IGHV3-11), where cell-line
deletions are most frequent. A flagged `common` variant is still strong
evidence that the variant exists in people: it was seen in hundreds of
individuals. What the flag puts in doubt is the exact frequency.

## gnomAD: `gnomad_filter` and `gnomad_flags`

| Value | Meaning | What it detects |
|---|---|---|
| `PASS` | No filter. | |
| `AC0` | No alternate allele remains after removing low-confidence genotypes (GQ < 20, DP < 10, het allele balance < 0.2). | The variant was only ever seen in poor genotypes. igpop classes it `absent`. |
| `AS_VQSR` | Failed gnomAD's allele-specific variant quality model. | Same idea as the VQSR filter above. |
| `InbreedingCoeff` | Inbreeding coefficient below −0.3. | A strong excess of heterozygotes: the collapsed-paralog signature. |
| flag `segdup` | The site lies in an annotated segmental duplication. | Mapping is ambiguous by construction. Most of IGH carries this flag. |
| flag `lcr` | The site lies in a low-complexity region. | Indel and alignment errors are frequent. |

## How to read score and site quality together

| Score | Site quality | Reading |
|---|---|---|
| high | passes | Population-supported variant with a usable frequency. |
| high | flagged | The variant exists in the population; its frequency is uncertain. |
| near 1 | either | The data do not decide. |
| low | either | Looks like an arbitrary change: a sequencing or inference error, a private variant, or a site short reads cannot genotype. |

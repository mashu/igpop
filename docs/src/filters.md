# Why an allele did not pass: the `filters` column

igpop writes two allele tables with the same columns:

| File | Content |
|---|---|
| `<prefix>.final.tsv` | Alleles that passed: verdict `supported` (by frequency or by haplotype) or `matches_reference`. |
| `<prefix>.alleles.tsv` | Every allele, passed or not, with the reasons in `filters`. |

Nothing is ever dropped silently: an allele missing from `final.tsv` is in
`alleles.tsv` with a `filters` entry, and each of its variants is in
`variants.tsv` with its own `filters` entry and all the numbers.

## Reasons for a variant

A variant that is not `supported` by frequency gets two reasons, separated by a comma:
first **what the data show**, then **how that scores**. A `supported`
variant has `.`.

What the data show (exactly one, checked in this order):

| Reason | Meaning | Where to look |
|---|---|---|
| `no_population_data` | No dataset could be consulted for this position. | Report, section 4. |
| `never_observed` | Consulted, and not seen in any dataset. | `kgp_AN`, `gnomad_AN`: how many chromosomes were looked at. |
| `single_observation` | Seen, but never more than once in any population. One copy can be a genotyping error. | `kgp_AC`, `gnomad_AC`, `dbsnp_freq`. |
| `rare_in_every_population` | Seen repeatedly, and every population large enough to tell puts the frequency below `--common` (1%). The variant is real-looking but rare. | `support`, `kgp_any_pop_AF_upper95`, `gnomad_any_group_AF_upper95`. |
| `common_not_demonstrated` | Seen repeatedly; the data neither show that it reaches 1% somewhere nor rule it out. Usually a variant near 1% in one population, or data too thin. | `class_p`, `kgp_max_pop`, `gnomad_max_group`. |

How that scores (exactly one):

| Reason | Meaning |
|---|---|
| `low_score` | The score is at most `--confirm-score` (0.1): at this level of support, arbitrary base changes outnumber variants of trusted alleles at least ten to one. Verdict `needs_confirmation`. |
| `intermediate_score` | The score lies between the two limits. Verdict `inconclusive`. |
| `not_scored` | An insertion or deletion; these have no score, the verdict follows the first reason. |

Typical pairs:

| `filters` | Reading |
|---|---|
| `never_observed,low_score` | Looks like an error or a private variant. |
| `rare_in_every_population,low_score` | Exists in people, at a frequency where errors are far more numerous than real allele-defining variants. |
| `rare_in_every_population,intermediate_score` | Rare but plausible. |
| `common_not_demonstrated,intermediate_score` | Borderline; often around 0.5 to 1% in one population. |

## Reasons for an allele

The allele's `filters` lists the reasons of its failing variants with
counts, in the same order:

```
rare_in_every_population:1,common_not_demonstrated:1,low_score:1,intermediate_score:1
```

reads: two variants failed; one is rare everywhere and scores low, one is
borderline. The `reason` column names them. Alleles that were not placed
have `no_coordinates` or `low_identity` instead.

## What is not a filter

These are reported, but do not make an allele fail, because they describe
the reliability of the data, not the allele:

- site quality in 1000 Genomes (`kgp_site_qc`, `site_qc_flagged`; see
  [Site quality filters](site-qc.md));
- gnomAD's own flags (`gnomad_filter`, `gnomad_flags`);
- an ambiguous or shared locus (`notes`);
- a haplotype check without an answer (`incomplete`, `untestable`). The
  check can only help an allele: found on enough chromosomes, it makes
  the allele `supported` ([The haplotype check](haplotypes.md)).

To be stricter, filter `alleles.tsv` yourself, e.g. keep only
`supported_by == frequency+haplotype`, or drop alleles with a `notes`
entry.

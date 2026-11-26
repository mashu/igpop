# Output files

Start with `<prefix>.report.md`; it is described in
[The run report and its checks](report.md). The tables below hold the
same results in full.

All tables are tab-separated with a header line. A `.` means "not
applicable or not available"; `-` in `evidence` means "no record".

| File | One row per | Holds |
|---|---|---|
| `variants.tsv` | difference of an allele | verdict, score, class, and where the evidence comes from |
| `populations.tsv` | difference and group of people | copies, chromosomes, frequency with its bounds |
| `site_quality.tsv` | difference and dataset | what the dataset says about the site |
| `datasets.tsv` | group of people of the run | size, whether it was used and why not |
| `alleles.tsv`, `final.tsv` | allele | verdict, reason, haplotype check |
| `score_model.tsv`, `score_levels.tsv` | — | the score curve of the run |
| `haplotype_calibration.tsv` | number of chromosomes | what a found haplotype is worth |
| `alignments.tsv`, `alignments.txt` | allele | where it was placed |

A difference is identified everywhere by the same five columns: `build`,
`contig`, `genome_pos`, `reference_bases`, `query_variant`. Use them to
join `variants.tsv` with `populations.tsv` and `site_quality.tsv`.

## `<prefix>.variants.tsv`

One row per difference between an allele and the reference. An allele
identical to the reference has no rows here (it still appears in
`alleles.tsv`).

### Identity of the variant

| Column | Meaning |
|---|---|
| `allele_name` | Allele name from the FASTA. |
| `gene` | Locus the allele was placed on (after `--gene-aliases`). |
| `build` | Genome build the allele was placed and judged on: `--build`, or another build for a gene the first one lacks. |
| `rsid` | dbSNP identifier of exactly this position and base change, or `.`. A label only: having an rsID is not evidence that a variant is real. |
| `evidence` | Datasets in which the variant was observed at least once, e.g. `1kGP_30x;gnomAD_v4.1;ToMMo_61k;dbSNP` (`dbSNP` = a study relayed by dbSNP reports a non-zero frequency), or `-` if none. The counts are in `populations.tsv`. |
| `contig` | Reference contig. |
| `genome_pos` | 1-based reference position. For insertions and deletions this is the base before the event, as in VCF. |
| `aln_pos` | 1-based position in the allele sequence as supplied. |
| `aln_pos_genome_strand` | The same base counted along the forward strand of the genome. Equal to `aln_pos` for a gene on the forward strand; for a gene on the reverse strand it counts from the other end of the allele (allele length - `aln_pos` + 1). IgSNPer reported this number. |
| `strand` | `+` or `-`: the strand of the gene. |
| `type` | `SNV`, `INS` or `DEL`. |
| `reference_bases` | The bases of the reference genome that the difference replaces, on the reference strand. For an insertion or deletion they begin with the base before it, as in VCF. |
| `query_variant` | The allele's base(s), written on the **reference strand** (for a minus-strand gene this is the complement of what the allele FASTA shows). |

### Verdict, score and class

| Column | Meaning |
|---|---|
| `verdict` | The verdict of this variant by frequency: `supported`, `inconclusive`, `needs_confirmation` or `not_assessed`. `supported`: `class` is `common`, or `score` is at least `--supported-score` (2). `needs_confirmation`: `score` at most `--confirm-score` (0.1), or an unscored variant that is `absent`. See [Reading the results](reading-results.md). |
| `score` | The main result: how many times more typical this support is of trusted variants than of arbitrary base changes. Above 1 supports the variant, below 1 speaks against it. See [The score](score.md). |
| `data_show` | What the data show about a variant that is not `common`, whatever its verdict: `never_observed`, `single_observation`, `rare_in_every_population`, `common_not_demonstrated` or `no_population_data`. The first entry of `filters` when the variant is not supported. |
| `score_basis` | `curve` for a single-base change; `borrowed` for an insertion or deletion, whose score is read from the same curve as a lower bound; `.` without a score. |
| `score_low`, `score_high` | A 95% interval for the score: how well the calibration knows it at this support. |
| `chance` | Share of all possible single-base changes in the genes of the run that are supported at least as well as this variant. |
| `filters` | Why the variant is not `supported`: what the data show, then how that scores, e.g. `rare_in_every_population,low_score`. `.` for a supported variant. See [Why an allele did not pass](filters.md). |
| `support` | The number the score is based on: the largest lower 95% bound of the allele frequency over all groups of all datasets. 0 if never seen. It is the largest `AF_lower95` of the variant in `populations.tsv`. |
| `support_from` | The group that gives it, as `dataset:group` (`gnomAD_v4.1:afr`, `ToMMo_61k:Japanese`, `dbSNP:TOPMED`). |
| `class` | `common`, `uncommon`, `absent` or `not_assessed` (below). |
| `class_basis` | What made the variant `common`: the dataset and population that passed the test (`gnomAD_v4.1:afr`), or `dbSNP:COMMON`. |
| `class_p` | p-value of the best group, multiplied by the number of groups of the run (Bonferroni; the number of rows with `used = yes` for this build in `datasets.tsv`). `.` if the allele was nowhere seen twice. See [Method](method.md#7-the-class-is-the-difference-demonstrably-common). |
| `site_quality` | What the datasets say against the site, as `dataset:criterion,...;dataset:...`, or `.` if nothing. Reported, never used in a verdict. Details in `site_quality.tsv`. |
| `dbsnp_common` | `yes` if dbSNP's COMMON flag applies to exactly this allele. |

| Class | Meaning |
|---|---|
| `common` | In at least one group the allele was seen at least twice **and** an exact binomial test rejects "frequency below `--common`" (default 1%) at 5% after Bonferroni correction for the number of groups of the run. Or: dbSNP flags exactly this allele as COMMON and no dataset read directly shows it to be rare everywhere. |
| `uncommon` | Seen in at least one dataset, but not shown to reach `--common` anywhere. |
| `absent` | Looked for by at least one dataset that examines every position, and seen in none. |
| `not_assessed` | No dataset could be consulted, and no COMMON flag. |

The test adapts to dataset size: three copies in 200 chromosomes do not
pass, 300 in 20,000 do. "Seen at least twice" guards against a single
miscalled genotype. The dbSNP COMMON flag is a curated external label, not
a test result; it applies to one allele per site, the most frequent
alternate allele.

### Locus context

| Column | Meaning |
|---|---|
| `locus_status` | `ok`, or `ambiguous_locus` if the allele aligns at least as well to another locus. |
| `shares_locus_with` | Other genes annotated at exactly the same interval. Population data at such a locus mixes reads from all of them. |

## `<prefix>.populations.tsv`

The evidence itself: one row per difference and per group of people that
looked at its position. Everything igpop concludes about a difference is
computed from its rows here.

| Column | Meaning |
|---|---|
| `build`, `contig`, `genome_pos`, `reference_bases`, `query_variant` | The difference, as in `variants.tsv`. |
| `rsid` | Its dbSNP identifier, or `.`. |
| `dataset` | `1kGP_30x`, `gnomAD_v4.1`, `ALFA`, `GenomeAsia`, `ToMMo_61k`, `WBBC`, or `dbSNP` for a study relayed by dbSNP (on GRCh37: `1kGP_phase3`, `gnomAD_v2.1`, `dbSNP`). |
| `group` | The group of people within the dataset: `AFR`, `nfe`, `Japanese`, `South China`, `TOPMED`, ... Described in `datasets.tsv`. |
| `counts` | `reported`: the dataset gives copies and chromosomes. `estimated`: a study relayed by dbSNP gives a frequency only; the counts are worked out from it and from the estimated size of the study. |
| `AC`, `AN` | Copies of the allele, and chromosomes examined at this position. `AC` is 0 where a dataset that examines every position has no record of the difference. |
| `AF` | `AC / AN`. |
| `AF_lower95`, `AF_upper95` | The lowest and the highest frequency the counts are compatible with (exact one-sided 95% bounds). |

A difference has no row for a group that says nothing about it: a group
of ALFA or a relayed study without a record at the position.

## `<prefix>.site_quality.tsv`

One row per difference and per dataset that has a record of it: what the
dataset itself says about the site. See
[Site quality filters](site-qc.md).

| Column | Meaning |
|---|---|
| `build` ... `rsid` | The difference, as above. |
| `dataset` | The dataset the row is about. |
| `failed` | `PASS`, or the criteria the site fails in this dataset: `FILTER=...` (the dataset's own filter), `missingness>5%`, `HWE`, `mendelian_errors>5%`. Only measures the dataset publishes are judged. |
| `filter` | The FILTER field of the record. |
| `region_flags` | Flags on the record, e.g. `segdup` (segmental duplication), `lcr` (low-complexity sequence). |
| `AC`, `AN` | Copies and chromosomes over all people of the dataset. |
| `people_one_copy`, `people_two_copies` | Heterozygous and homozygous carriers, where published. |
| `hwe_p` | Hardy-Weinberg p-value per group, where published. |
| `excess_het_p` | p-value for an excess of heterozygotes (small = excess). |
| `mendelian_error_rate` | Share of parent-child trios whose genotypes are impossible by inheritance (1000 Genomes). |
| `vqslod` | The call set's own log-odds that the site is a true variant. |

## `<prefix>.datasets.tsv`

One row per group of people the run could consult: what the evidence was
drawn from. The rows with `used = yes` for a build are the strata of that
build; their number is the factor of the Bonferroni correction.

| Column | Meaning |
|---|---|
| `build` | Genome build the row belongs to. |
| `dataset`, `dataset_title` | Short name as in the other tables, and the name in words. |
| `group`, `group_title` | The group as in the other tables, and in words. |
| `ancestry` | Continental group, for colour and orientation only: `african`, `american`, `east_asian`, `european`, `middle_eastern`, `south_asian`, `other`. |
| `counts` | `reported` or `estimated`, as in `populations.tsv`. |
| `chromosomes` | The most chromosomes the group has at any position read (for an estimated study: its estimated size). |
| `used`, `reason` | Whether the group is evidence in this run, and if not, why: left out by name, a copy of a dataset read directly, too small to show a frequency of `--common`, or no chromosomes in the genes read. |
| `every_position_examined` | `yes` if a difference without a record counts as "looked for and not seen" (0 copies). |
| `records` | Records of the dataset in the regions read. |
| `smallest_nonzero_AF` | For an estimated study: the smallest frequency it reports, from which its size is estimated. |
| `people` | People in the dataset, where it has a fixed number. |
| `read_from` | File or address the dataset was read from. |

## `<prefix>.alleles.tsv` and `<prefix>.final.tsv`

One row per allele. `alleles.tsv` holds every allele; `final.tsv` holds
the rows with verdict `supported` (by frequency or by haplotype) or
`matches_reference`, same columns.

| Column | Meaning |
|---|---|
| `allele_name`, `gene`, `build` | As above. |
| `verdict` | By frequency, the verdict of the allele's weakest variant; raised to `supported` when the haplotype check supports the allele. Or `matches_reference` (no variants), or `not_placed`. |
| `supported_by` | `frequency`, `haplotype` or `frequency+haplotype` for a supported allele; `.` otherwise. |
| `filters` | Reasons of the variants that fail by frequency, with counts, e.g. `never_observed:1,low_score:1`; `no_coordinates` or `low_identity` if not placed; `.` if none. An allele supported by haplotype keeps them: they say why frequency alone did not suffice. |
| `hap_status` | Outcome of the [haplotype check](haplotypes.md): `found_exact`, `found_with_extra`, `not_found`, `incomplete`, `untestable`, `not_checked`. |
| `reason` | The variants behind the verdict, or why the allele was not placed. |
| `notes` | Circumstances to keep in mind: ambiguous or shared locus, truncated allele, variants at sites failing site QC. |
| `locus_status` | `ok`, `ambiguous_locus`, `low_identity`, or `no_coordinates`. |
| `shares_locus_with` | As above. |
| `trusted` | Whether the allele is named as a trusted example of real variation (by name or `--trusted`) and was placed cleanly. |
| `yardstick` | Why the allele calibrates the score: `named` (trusted), `haplotype` (two or more differences confirmed together on phased chromosomes), or `.`. |
| `missing_5prime`, `missing_3prime` | Bases of the gene interval that the allele does not reach at its 5' and 3' ends. A truncated allele is reported here, not as extra variants. |
| `n_variants` | Number of differences from the reference. |
| `min_score`, `min_score_variant` | The lowest score among the allele's single-base variants, and which variant it is. An allele is only as well supported as its weakest variant. |
| `common`, `common_str` | Number and list of variants with `class = common`. |
| `uncommon`, `uncommon_str` | Number and list of all variants **not** classed common (`uncommon` + `absent` + `not_assessed`). Same meaning as in IgSNPer. |
| `absent`, `absent_str` | The subset never seen in any dataset. |
| `not_assessed` | The subset for which no dataset could be consulted. |
| `site_qc_flagged`, `site_qc_flagged_str` | Variants at a site that fails a criterion in some dataset (`site_quality` in `variants.tsv` is not `.`). These overlap the classes above. |
| `hap_tested` | Variants that are sites of the phased panel / all variants. |
| `hap_untested_str` | The variants that are not. |
| `hap_total` | Founder chromosomes in the panel. |
| `hap_exact` | Chromosomes carrying all variants and no other panel variant inside the allele. `.` unless every variant was tested. |
| `hap_exact_donors` | The people carrying those chromosomes. |
| `hap_coincidence` | Share of arbitrary single-base changes of real chromosomes whose sequence is found on at least `hap_exact` chromosomes. Below `--coincidence` (0.05), with at least two chromosomes, the allele is supported by haplotype. |
| `hap_exact_AF`, `hap_exact_AF_lower95`, `hap_exact_AF_upper95` | Its frequency with exact 95% bounds. |
| `hap_exact_AFR` ... `hap_exact_SAS` | Carriers / founder chromosomes per superpopulation. |
| `hap_exact_max_pop` | Superpopulation with the highest lower bound. |
| `hap_exact_pop_heterogeneity_p` | Test of equal frequency in the five superpopulations; small = population-specific. |
| `hap_carrying`, `hap_carrying_AF`, `hap_carrying_AFR` ... | Chromosomes carrying all tested variants, with or without others. |

Variant lists use `rsid(REF>ALT:aln_pos,genome_pos)`, separated by `;`.

`uncommon` and `min_score` tell the same story at different resolution:
more uncommon variants go with a lower `min_score`, but two alleles with
the same count can differ several hundred-fold in score
([The score](score.md)).

## `<prefix>.score_model.tsv`

The curve behind the score, as `key`, `value` rows. Pass it to a later run
with `--calibration`.

| Key | Meaning |
|---|---|
| `never` | Score of a change no dataset has seen. |
| `b0`, `b1`, `b2` | log(score) = b0 + b1 x + b2 x² for a change that was seen, x = log10(support). `b2` is 0 when a straight line fits as well. |
| `cov` | Covariance of the three coefficients (six numbers: the upper triangle), from which `score_low` and `score_high` are computed. |
| `arbitrary`, `arbitrary_seen` | Number of possible single-base changes in the loci, and how many of them some dataset has seen. |
| `trusted`, `trusted_seen` | The same for the variants of yardstick alleles. |
| `tail` | 49 rows: a support, and the share of all possible changes supported at least that well (column `chance` is read from these). |
| `real` | One row per real variant that was seen: its support. The points the curve was drawn through. |

## `<prefix>.score_levels.tsv`

The two groups counted per level of support: a ruler for reading the
curve. The levels do not enter the score.

| Column | Meaning |
|---|---|
| `support` | Level of support: the range of the lower frequency bound. |
| `arbitrary_changes`, `arbitrary_fraction` | How many of all possible single-base changes in the loci fall in this level, and their share. |
| `trusted_variants`, `trusted_fraction` | The same for the variants of trusted alleles. |
| `ratio_of_fractions` | The raw ratio of the two shares in this level. |
| `score_from`, `score_to` | The score at the lower and at the upper edge of the level. |
| `score_source` | `run`: the curve was fitted to this run. `built-in`: the counts are those of this run, the curve is the built-in one. `file`: curve from `--calibration`. See [Automatic calibration](calibration.md). |

## `<prefix>.haplotype_calibration.tsv`

How much a found haplotype is worth in this run
([The haplotype check](haplotypes.md)).

| Column | Meaning |
|---|---|
| `chromosomes` | A number of founder chromosomes. |
| `coincidence` | Share of arbitrary single-base changes of real chromosomes whose sequence is found on at least that many chromosomes. |
| `supports` | Whether an allele found on that many chromosomes is supported by its haplotype. |

## `<prefix>.alignments.tsv`

One row per allele: where it was placed.

| Column | Meaning |
|---|---|
| `allele`, `gene`, `build`, `contig` | Identity and location. |
| `ref_start`, `ref_end` | Reference interval the allele covers. Positions outside it are unknown, not reference. |
| `strand` | Strand of the reference the allele aligned to. |
| `identity` | Matching columns divided by alignment columns. |
| `score` | Alignment score at its own locus. |
| `best_other_gene`, `best_other_score` | Best-scoring different locus and its score. |
| `shares_locus_with` | As above. |
| `missing_5prime`, `missing_3prime` | Gene bases not reached by the allele at each end. |
| `status` | `ok`; `ambiguous_locus` (another locus scores at least as well); `low_identity` (below `--min-identity`, no variants emitted); `no_coordinates` (gene not in the coordinates table for this build). |

## `<prefix>.alignments.txt`

Alignments for reading by eye, one block per gene, in gene (5'→3')
orientation. The reference is on top, and each allele is a row under it:

```
>IGHV1-3  GRCh38 chr14:106005095-106005390  strand -  2 allele(s)
               106005390 .. 106005291
reference      CAGGTCCAGCTTGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCC…
IGHV1-3*01     .....T.....G....................................…
IGHV1-3*02     ................................................…

  IGHV1-3*01  identity 0.9797  status ok  variants 6  insertions none
```

| Symbol | Meaning |
|---|---|
| `.` | Same base as the reference. |
| letter | Mismatch: the allele's base. |
| `-` | Reference base missing in the allele (deletion). |
| blank | Not covered by the allele. |

The two numbers above each chunk are the reference positions of its first
and last column. Inserted bases cannot be drawn in a reference-anchored
stack, so they are listed under the block as
`<reference position before the insertion>+<bases on the reference strand>`.

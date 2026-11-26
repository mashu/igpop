// The fixed words of the tables in plain language, and how numbers are
// written.

import { el } from './dom.js?v=BUILD';
import { given } from './tsv.js?v=BUILD';

export const VERDICTS = {
  supported: { label: 'Supported', icon: '✓', help: 'Every difference from the reference is an established population variant, or the exact sequence of the allele is established on phased chromosomes.' },
  inconclusive: { label: 'Inconclusive', icon: '?', help: 'At least one difference is seen in people, but too rarely to tell a real variant from an error.' },
  needs_confirmation: { label: 'Needs confirmation', icon: '!', help: 'At least one difference looks like an arbitrary base change. Confirm it independently.' },
  not_assessed: { label: 'Not assessed', icon: '–', help: 'No population data could be read for at least one difference.' },
  matches_reference: { label: 'Matches reference', icon: '=', help: 'No difference from the reference genome.' },
  not_placed: { label: 'Not placed', icon: '×', help: 'The allele could not be aligned to a locus.' },
};

const FILTERS = {
  no_population_data: 'no population data',
  never_observed: 'never observed',
  single_observation: 'seen only once',
  rare_in_every_population: 'rare in every population',
  common_not_demonstrated: 'not shown to reach 1%',
  low_score: 'low score',
  intermediate_score: 'intermediate score',
  not_scored: 'not scored (insertion or deletion)',
  sequence_not_in_panel: 'each difference is known, but no phased chromosome has this exact sequence',
  no_coordinates: 'no coordinates for its gene',
  low_identity: 'too different from its gene',
};

// How well a sequence could be placed on the gene its name gives.
export const PLACEMENTS = {
  ok: { label: 'Placed on its gene', help: 'The sequence fits the gene of its name better than any other gene of the list.' },
  ambiguous_locus: { label: 'Fits another gene too', help: 'The sequence fits another gene of the list at least as well as the gene of its name. The genes are near-identical copies, so its differences, and their frequencies, may belong to the other gene.' },
  low_identity: { label: 'Too different from its gene', help: 'The sequence is too different from the gene of its name to compare; no differences are called.' },
  no_coordinates: { label: 'Gene not in the list', help: 'The gene of its name is not in the gene list, so the sequence could not be placed.' },
};

export const HAPLOTYPES = {
  found_exact: { label: 'Haplotype found', help: 'Chromosomes with exactly this sequence exist in the phased 1000 Genomes panel.' },
  found_with_extra: { label: 'Found with extra variants', help: 'The combination exists, but always together with further variants inside the allele.' },
  not_found: { label: 'Haplotype not found', help: 'Every variant is in the phased panel, but no chromosome carries them all together.' },
  incomplete: { label: 'Partly testable', help: 'Some variants are not sites of the phased panel, so the full combination cannot be looked up.' },
  untestable: { label: 'Not testable', help: 'None of the variants is a site of the phased panel.' },
  not_checked: { label: 'Not checked', help: 'Not placed, or no phased data for this build.' },
};

export const OUTPUTS = ['report.md', 'final.tsv', 'alleles.tsv', 'variants.tsv', 'alignments.tsv', 'alignments.txt', 'score_model.tsv', 'score_levels.tsv', 'haplotype_calibration.tsv', 'populations.tsv', 'site_quality.tsv', 'datasets.tsv'];

export const ROUTES = {
  frequency: { label: 'by frequency', help: 'Each difference is an established population variant.' },
  haplotype: { label: 'by haplotype', help: 'The exact sequence is established on phased chromosomes; by frequency alone it would not pass.' },
  'frequency+haplotype': { label: 'by frequency and haplotype', help: 'Each difference is an established variant, and the exact sequence is established on phased chromosomes.' },
};

export const FILE_NOTES = {
  'report.md': 'Results and checks in words',
  'final.tsv': 'Alleles that passed',
  'alleles.tsv': 'Every allele: verdict, reasons, haplotype check',
  'variants.tsv': 'Every difference with all population evidence',
  'alignments.tsv': 'Where each allele was placed',
  'alignments.txt': 'Alignments to read by eye',
  'score_model.tsv': 'The curve behind the score',
  'score_levels.tsv': 'Real variants and possible changes counted per level of support',
  'haplotype_calibration.tsv': 'How often an arbitrary change is found on phased chromosomes',
  'populations.tsv': 'Every count behind the verdicts: one row per difference and population',
  'site_quality.tsv': 'How far each dataset\u2019s counts can be trusted at each difference',
  'datasets.tsv': 'The populations of the run: who they are and how many chromosomes',
};

export const sortGenes = (a, b) => a.localeCompare(b, 'en', { numeric: true });

export function badge(verdict) {
  const v = VERDICTS[verdict] || { label: verdict, icon: '' };
  return el('span', { class: `badge v-${verdict}`, title: v.help }, el('i', { 'aria-hidden': 'true', text: v.icon }), v.label);
}

export function chips(filters) {
  if (!given(filters)) return [];
  return filters.split(',').map((part) => {
    const [code, count] = part.split(':');
    const label = FILTERS[code] || code;
    return el('span', { class: 'chip', text: count && Number(count) > 1 ? `${label} ×${count}` : label });
  });
}

// Frequencies are always percentages: one decimal from 1% up, two
// significant digits below (0.70%, 0.0074%).
export const percent = (af) => (!(af > 0) ? '0%' : af >= 0.01 ? `${(100 * af).toFixed(1)}%` : `${Number((100 * af).toPrecision(2))}%`);

export const scoreText = (x) => (x >= 10 ? x.toFixed(0) : x >= 1 ? x.toFixed(1) : x >= 0.1 ? x.toFixed(2) : x.toPrecision(1));

export const pct = (x) => (x >= 0.995 ? '100%' : `${Number((100 * x).toPrecision(2))}%`);

// Why the genotypes of a dataset at a position are less reliable: the
// site criteria of site_quality.tsv, in plain words.
export const SITE_QC = [
  [/^FILTER=/, 'the dataset\u2019s own variant filter did not pass this call: the reads that support it look less like a true variant than usual'],
  [/^missingness/, 'more than 5% of the people have no genotype here, usually because too few reads could be placed'],
  [/^HWE/, 'the numbers of carriers of one and of two copies do not fit together in any population, as happens when reads of a second copy of the gene are mixed in'],
  [/^mendelian/, 'in more than 5% of the families the child\u2019s genotype does not fit the parents\u2019'],
];

/** Variant ids of the tables, such as rs1(A>T:163,106805343), in words. */
export const idWords = (text) => text.replace(/([^\s;,(]+)\(([A-Z-]+)>([A-Z-]+):(\d+),\d+\)/g, (_, id, from, to, at) => `${id === '.' ? 'a difference without rsID' : id} (${from} to ${to} at position ${at})`).replace(/;/g, ', ');

export const BASE_WORDS = { '.': 'same as the reference', '-': 'base missing (deletion)', ' ': 'not covered by this allele' };


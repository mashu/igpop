// The limit from which a score counts as support: a setting of the
// page, applied to the tables of the run without running again.

import { $, el, explain } from './dom.js?v=BUILD';
import { scoreModel } from './score.js?v=BUILD';
import { emit, state } from './state.js?v=BUILD';
import { given } from './tsv.js?v=BUILD';
import { pct, percent, scoreText } from './words.js?v=BUILD';

const RUN_LIMIT = 2; // the limit the analysis itself is run with

export const CONFIRM = 0.1;

const FILTER_ORDER = ['no_population_data', 'never_observed', 'single_observation', 'rare_in_every_population', 'common_not_demonstrated', 'low_score', 'intermediate_score', 'not_scored'];

const SCORE_FILTERS = ['low_score', 'intermediate_score', 'not_scored'];

/** Keep what the analysis decided, so that the limit can be moved and
 *  moved back. */
export function rememberVerdicts(files) {
  state.limit = RUN_LIMIT;
  state.runFiles = { 'alleles.tsv': files['alleles.tsv'], 'variants.tsv': files['variants.tsv'], 'final.tsv': files['final.tsv'] };
  for (const a of state.alleles) a.run = { verdict: a.verdict, supported_by: a.supported_by, filters: a.filters, reason: a.reason };
  for (const rows of state.variants.values()) for (const v of rows) v.run = { verdict: v.verdict, filters: v.filters };
}

/** The verdict of one difference by frequency, as the program decides it. */
function variantVerdict(v, limit) {
  if (v.class === 'common') return 'supported';
  if (v.class === 'not_assessed') return 'not_assessed';
  if (given(v.score)) return Number(v.score) >= limit ? 'supported' : Number(v.score) <= CONFIRM ? 'needs_confirmation' : 'inconclusive';
  return v.class === 'absent' ? 'needs_confirmation' : 'inconclusive';
}

/** Judge every difference and allele again with another limit for
 *  "supported". Scores, haplotype checks and everything measured stay as
 *  they are; only the reading changes. */
function applyLimit(limit) {
  state.limit = limit;
  const rerun = limit === RUN_LIMIT;
  for (const a of state.alleles) {
    const variants = state.variants.get(a.allele_name) || [];
    for (const v of variants) {
      v.verdict = rerun ? v.run.verdict : variantVerdict(v, limit);
      if (rerun || v.verdict === v.run.verdict) v.filters = v.run.filters;
      else if (v.verdict === 'supported') v.filters = '.';
      else {
        // What the data show stays; how that scores is read again.
        const shown = given(v.data_show) ? [v.data_show] : given(v.run.filters) ? v.run.filters.split(',').filter((f) => !SCORE_FILTERS.includes(f)) : ['common_not_demonstrated'];
        v.filters = [...shown, given(v.score) ? (Number(v.score) <= CONFIRM ? 'low_score' : 'intermediate_score') : 'not_scored'].join(',');
      }
    }
    if (rerun || ['not_placed', 'matches_reference'].includes(a.run.verdict)) { Object.assign(a, a.run); continue; }
    const verdicts = variants.map((v) => v.verdict);
    const byFrequency = ['needs_confirmation', 'inconclusive', 'not_assessed'].find((x) => verdicts.includes(x)) || 'supported';
    // Known from the run: whether the haplotype check supports the allele,
    // and whether a complete lookup found its sequence on no chromosome.
    const haplotype = (a.run.supported_by || '').includes('haplotype');
    const absent = given(a.hap_exact) && Number(a.hap_exact) === 0;
    const words = { needs_confirmation: 'without population support', inconclusive: 'with weak population support', not_assessed: 'without population data' };
    const weakest = variants.filter((v) => v.verdict === byFrequency);
    const frequencyReason = byFrequency === 'supported'
      ? `all variants (${variants.length}) have population support`
      : `${weakest.length} of ${variants.length} variant(s) ${words[byFrequency]}: ${weakest.map((v) => (given(v.rsid) ? v.rsid : `${v.contig}:${v.genome_pos}`)).join(';')}`;
    const counted = FILTER_ORDER.map((f) => [f, variants.filter((v) => given(v.filters) && v.filters.split(',').includes(f)).length]).filter(([, k]) => k > 0);
    if (byFrequency === 'supported' && absent) {
      Object.assign(a, { verdict: 'inconclusive', supported_by: '.', filters: 'sequence_not_in_panel', reason: `every variant has population support, but none of the ${given(a.hap_total) ? a.hap_total : ''} phased chromosomes has this exact sequence` });
    } else if (byFrequency !== 'supported' && haplotype) {
      Object.assign(a, { verdict: 'supported', supported_by: 'haplotype', filters: counted.map(([f, k]) => `${f}:${k}`).join(',') || '.',
        reason: a.run.supported_by === 'haplotype' ? a.run.reason : `exact sequence on phased chromosomes; by frequency alone: ${frequencyReason}` });
    } else {
      Object.assign(a, { verdict: byFrequency, supported_by: byFrequency === 'supported' ? (haplotype ? 'frequency+haplotype' : 'frequency') : '.',
        filters: counted.map(([f, k]) => `${f}:${k}`).join(',') || '.', reason: frequencyReason });
    }
  }
  // The tables to download follow the limit.
  if (rerun) Object.assign(state.files, state.runFiles);
  else {
    const table = (text, rows) => {
      const header = text.split('\n')[0].split('\t');
      return `${[header.join('\t'), ...rows.map((row) => header.map((column) => row[column] ?? '.').join('\t'))].join('\n')}\n`;
    };
    state.files['alleles.tsv'] = table(state.runFiles['alleles.tsv'], state.alleles);
    state.files['final.tsv'] = table(state.runFiles['final.tsv'] || state.runFiles['alleles.tsv'], state.alleles.filter((a) => ['supported', 'matches_reference'].includes(a.verdict)));
    state.files['variants.tsv'] = table(state.runFiles['variants.tsv'], state.alleles.flatMap((a) => state.variants.get(a.allele_name) || []));
  }
}

/** The setting itself, and what the chosen limit costs and gains. */
export function renderLimit() {
  const model = scoreModel();
  $('limit-setting').hidden = !model;
  if (!model) return;
  $('limit').value = String(Math.min(20, state.limit));
  $('limit-value').value = String(state.limit);
  $('limit-reset').hidden = state.limit === RUN_LIMIT;
  // The smallest support that passes; a limit at or below the score of a
  // base never seen lets everything pass.
  const everything = model.score(0) >= state.limit;
  const at = everything ? 0 : model.supportAt(state.limit);
  const supported = state.alleles.filter((a) => a.verdict === 'supported').length;
  const before = state.alleles.filter((a) => a.run.verdict === 'supported').length;
  const change = supported === before ? '' : ` (${supported > before ? `${supported - before} more` : `${before - supported} fewer`} than at ${RUN_LIMIT})`;
  const reading = state.limit >= 10 ? 'From 10 on, a score is usually read as strong evidence.'
    : state.limit >= 3 ? 'Between 3 and 10, a score is usually read as positive but not strong evidence.'
      : state.limit > 1 ? 'Between 1 and 3, a score is weak evidence: such a difference is more typical of real variants than of arbitrary changes, but not by much.'
        : 'At 1 or below, a difference passes although its support is no more typical of real variants than of arbitrary changes.';
  // A figure, what it counts, and on request what exactly it means.
  const fact = (number, words, meaning) => el('span', {},
    explain(el('b', { class: 'hint-value', tabindex: '0', text: number }), el('span', {}, meaning.map((line) => el('span', { text: line })))), words);
  const judged = state.alleles.filter((a) => !['matches_reference', 'not_placed'].includes(a.verdict)).length;
  const counts = model.counts;
  $('limit-note').replaceChildren(
    everything
      ? fact('any', 'even a base never seen passes', [`A base that no dataset has seen scores ${scoreText(model.score(0))}, which is at or above this limit.`])
      : Number.isFinite(at)
        ? fact(percent(at), 'frequency a difference must be shown to reach', [
          `A difference scores ${state.limit} or more when the data guarantee, with 95% confidence, a frequency of at least ${percent(at)} in some population of some source.`,
          'This guaranteed frequency is lower than the observed one, the more so the fewer people were examined: 3 carriers among 200 chromosomes guarantee less than 300 among 20,000.',
          'A difference shown to reach 1% in a population is supported whatever this limit.'])
        : fact('none', 'no frequency reaches this score in this run', ['The score curve of this run stays below this limit even for a base that everyone carries.']),
    everything
      ? fact('100%', 'of arbitrary base changes would pass', ['With this limit, every single-base change that is possible in these genes would count as supported.'])
      : Number.isFinite(at)
        ? fact(pct(model.chance(at)), 'of arbitrary base changes would pass', [
          `Take every single-base change that is possible in the genes of this run: ${counts.arbitrary.toLocaleString('en')} of them, three at each position. ${pct(model.chance(at))} of them score ${state.limit} or more.`,
          'This is how often a random error in one base would pass by its score alone. It is the price of the limit: a lower limit lets more real differences through, and more errors with them.',
          'It counts single differences. An allele with several differences passes only if each of them does.'])
        : null,
    fact(`${supported} of ${judged}`, `alleles supported${change}`, [
      `${judged} of the ${state.alleles.length} alleles differ from the reference genome and could be placed on a gene; the others are identical to the reference or were not placed.`,
      `${supported} of them are supported with this limit: every difference passes by frequency, or the whole sequence was found on phased chromosomes.`]),
    el('p', { text: `${reading}${state.limit === RUN_LIMIT ? '' : ` Verdicts, counts and the tables under Files follow this limit; the report describes the run at ${RUN_LIMIT}.`}` }));
}

export function wireLimit() {
  const move = (value) => {
    const limit = Math.min(100, Math.max(0, Math.round(10 * Number(value)) / 10));
    if (!Number.isFinite(limit)) return renderLimit();
    applyLimit(limit);
    renderLimit();
    emit('verdicts');
  };
  $('limit').addEventListener('input', () => {
    // While dragging, only the number moves; the list follows on release.
    $('limit-value').value = $('limit').value;
  });
  $('limit').addEventListener('change', () => move($('limit').value));
  $('limit-value').addEventListener('change', () => move($('limit-value').value));
  $('limit-reset').addEventListener('click', () => move(RUN_LIMIT));
}


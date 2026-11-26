// The score of a difference: the curve of the run, and the picture that
// places a difference on it.

import { el, hint, plural, svg } from './dom.js?v=BUILD';
import { state } from './state.js?v=BUILD';
import { given, number, parseTSV } from './tsv.js?v=BUILD';
import { pct, percent, scoreText } from './words.js?v=BUILD';

/** The curve the scores of this run come from (score_model.tsv), with the
 *  real variants it was drawn through. Worked out once per result. */
export function scoreModel() {
  if (state.model !== undefined) return state.model;
  state.model = null;
  const text = state.files['score_model.tsv'];
  if (!text) return null;
  const lines = text.split('\n').map((line) => line.split('\t'));
  const value = (key) => (lines.find((f) => f[0] === key) || [])[1];
  const beta = ['b0', 'b1', 'b2'].map((key) => Number(value(key)));
  const never = Number(value('never'));
  const c = (value('cov') || '').split(',').map(Number);
  if (!Number.isFinite(never) || beta.some((b) => !Number.isFinite(b)) || c.length !== 6) return null;
  const cov = [[c[0], c[1], c[2]], [c[1], c[3], c[4]], [c[2], c[4], c[5]]];
  const counts = Object.fromEntries(['arbitrary', 'arbitrary_seen', 'trusted', 'trusted_seen'].map((key) => [key, Number(value(key))]));
  const curve = (support) => {
    const x = Math.log10(Math.min(1, Math.max(1e-6, support)));
    const basis = [1, x, x * x];
    const log = basis.reduce((sum, b, i) => sum + b * beta[i], 0);
    let variance = 0;
    for (let i = 0; i < 3; i += 1) for (let j = 0; j < 3; j += 1) variance += basis[i] * cov[i][j] * basis[j];
    return [log, Math.sqrt(Math.max(0, variance))];
  };
  const score = (support) => (support > 0 ? Math.max(never, Math.exp(curve(support)[0])) : never);
  const interval = (support) => {
    const [log, se] = curve(support);
    return [Math.max(Math.min(never, score(support)), Math.exp(log - 1.96 * se)), Math.max(never, Math.exp(log + 1.96 * se))];
  };
  // The real variants: with their names when this run made the curve
  // itself, otherwise only their supports as kept with the curve.
  const levels = state.files['score_levels.tsv'] ? parseTSV(state.files['score_levels.tsv']) : [];
  const own = levels.some((r) => r.score_source === 'run');
  const accepted = new Set(state.alleles.filter((a) => given(a.yardstick)).map((a) => a.allele_name));
  const named = new Map();
  if (own) for (const name of accepted) for (const v of state.variants.get(name) || []) {
    if (v.type === 'SNV' && v.build === state.build) named.set(`${v.contig}|${v.genome_pos}|${v.query_variant}`, v);
  }
  const real = own && named.size
    ? [...named.values()].map((v) => ({ support: number(v.support), label: `${given(v.rsid) ? v.rsid : 'no rsID'} in ${v.gene}` }))
    : lines.filter((f) => f[0] === 'real').map((f) => ({ support: Number(f[1]), label: 'a variant of a real allele' }));
  const tail = lines.filter((f) => f[0] === 'tail').map((f) => ({ support: Number(f[1]), share: Number(f[2]) }));
  // Share of all possible changes supported at least this well.
  const chance = (support) => (support > 0 ? (tail.filter((t) => t.support <= support * (1 + 1e-9)).pop() || { share: 1 }).share : 1);
  // The smallest support whose score reaches a limit.
  const supportAt = (limit) => {
    if (score(1) < limit) return NaN;
    let [low, high] = [-6, 0];
    for (let i = 0; i < 50; i += 1) {
      const middle = (low + high) / 2;
      if (score(10 ** middle) >= limit) high = middle; else low = middle;
    }
    return 10 ** high;
  };
  state.model = {
    never, beta, counts, score, interval, real, own, levels, chance, supportAt,
    tail,
    seenShare: counts.arbitrary_seen / counts.arbitrary,
    alleles: accepted.size,
    genes: new Set(state.alleles.filter((a) => accepted.has(a.allele_name)).map((a) => a.gene)).size,
  };
  return state.model;
}

/** A picture of what a difference was compared with. On one axis of
 *  support: the score curve with its uncertainty and the two limits, the
 *  real variants as dots, and all possible changes as bars. `v` marks one
 *  difference; without it the picture shows the curve alone. */
export function comparisonChart(v) {
  const model = scoreModel();
  if (!model) return null;
  const support = v ? number(v.support) : NaN;
  const W = 700;
  const [slot0, slot1, x0, x1] = [52, 104, 122, 688];
  // Three plots under each other: the same axis of support, an axis of
  // its own for each.
  const SCORE = { top: 40, bottom: 170 };
  const REAL = { top: 212, bottom: 268 };
  const ALL = { top: 310, bottom: 366 };
  const H = ALL.bottom + 38;
  const xOf = (s) => (s > 0 ? x0 + ((Math.log10(Math.min(1, Math.max(1e-6, s))) + 6) / 6) * (x1 - x0) : (slot0 + slot1) / 2);
  const [yLow, yHigh] = [Math.log10(0.002), Math.log10(300)];
  const yOf = (score) => SCORE.bottom - ((Math.log10(Math.min(300, Math.max(0.002, score))) - yLow) / (yHigh - yLow)) * (SCORE.bottom - SCORE.top);
  const picture = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'compare', role: 'img',
    'aria-label': v
      ? `Three plots over the same axis of support: the score, the real variants and all possible changes. This difference has a support of ${support > 0 ? percent(support) : 'zero (never seen)'} and a score of ${scoreText(Number(v.score))}.`
      : 'Three plots over the same axis of support: the score as a smooth curve, the real variants and all possible changes it was made from.' });
  const text = (x, y, cls, words, anchor = 'middle') => svg('text', { x, y, class: cls, 'text-anchor': anchor }, words);
  const frame = (plot, title) => {
    picture.append(svg('rect', { x: slot0, y: plot.top, width: slot1 - slot0, height: plot.bottom - plot.top, class: 'panel' }));
    picture.append(svg('rect', { x: x0, y: plot.top, width: x1 - x0, height: plot.bottom - plot.top, class: 'panel' }));
    picture.append(text(slot0, plot.top - 8, 'plot-title', title, 'start'));
    for (let decade = -6; decade <= 0; decade += 1) {
      picture.append(svg('line', { x1: xOf(10 ** decade), x2: xOf(10 ** decade), y1: plot.top, y2: plot.bottom, class: 'grid' }));
    }    // The frequency a variant must demonstrably reach to be called common.
    picture.append(svg('line', { x1: xOf(0.01), x2: xOf(0.01), y1: plot.top, y2: plot.bottom, class: 'one-percent' }));
  };

  // 1. The score.
  frame(SCORE, 'Score');
  for (const tick of [0.01, 0.1, 1, 10, 100]) {
    picture.append(svg('line', { x1: slot0, x2: x1, y1: yOf(tick), y2: yOf(tick), class: tick === 0.1 ? 'limit low' : 'grid' }));
    picture.append(text(slot0 - 5, yOf(tick) + 3, 'tick', String(tick), 'end'));
  }
  picture.append(svg('line', { x1: slot0, x2: x1, y1: yOf(state.limit), y2: yOf(state.limit), class: 'limit high' }));
  picture.append(text(x0 + 6, yOf(state.limit) - 4, 'limit-name high', `supported from ${state.limit}`, 'start'));
  picture.append(text(x1 - 4, yOf(0.1) + 11, 'limit-name low', 'needs confirmation up to 0.1', 'end'));
  const samples = Array.from({ length: 121 }, (_, i) => 10 ** (-6 + i / 20));
  const band = [...samples.map((s) => `${xOf(s).toFixed(1)},${yOf(model.interval(s)[1]).toFixed(1)}`),
    ...samples.reverse().map((s) => `${xOf(s).toFixed(1)},${yOf(model.interval(s)[0]).toFixed(1)}`)];
  samples.reverse();
  picture.append(svg('polygon', { points: band.join(' '), class: 'band' }));
  picture.append(svg('polyline', { points: samples.map((s) => `${xOf(s).toFixed(1)},${yOf(model.score(s)).toFixed(1)}`).join(' '), class: 'curve' }));
  picture.append(svg('line', { x1: slot0 + 8, x2: slot1 - 8, y1: yOf(model.never), y2: yOf(model.never), class: 'curve' }));

  // 2 and 3. How the two groups spread over support: the share of the
  // group in each eighth of a step of ten, and the share never seen.
  const step = (x1 - x0) / 48;
  const histogram = (plot, title, shares, never, cls, noun) => {
    frame(plot, title);
    const tall = plot.bottom - plot.top;
    // The axis ends at a round share just above the tallest bar.
    const most = Math.max(...shares, 1e-9);
    const top = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1].find((edge) => edge >= most) || 1;
    for (const [share, y] of [[0, plot.bottom], [top, plot.top]]) picture.append(text(slot0 - 5, y + (share ? 8 : 0), 'tick', `${Number((100 * share).toPrecision(1))}%`, 'end'));
    shares.forEach((share, k) => {
      if (!(share > 0)) return;
      const height = Math.max(1, (share / top) * tall);
      const bar = svg('rect', { x: x0 + k * step + 0.5, y: plot.bottom - height, width: step - 1, height, class: cls });
      bar.append(svg('title', {}, `${pct(share)} of ${noun} have a support from ${percent(10 ** (-6 + k / 8))}`));
      picture.append(bar);
    });
    // Never seen: a bar of its own, with its share written on it.
    const height = never > 0 ? Math.max(1, Math.min(1, never / top) * (tall - 12)) : 0;
    if (height) picture.append(svg('rect', { x: slot0 + 10, y: plot.bottom - height, width: slot1 - slot0 - 20, height, class: cls }));
    picture.append(text((slot0 + slot1) / 2, plot.bottom - height - 3, 'count', never > 0 ? pct(never) : 'none'));
    picture.append(svg('line', { x1: slot0, x2: x1, y1: plot.bottom, y2: plot.bottom, class: 'base' }));
  };
  const realShares = Array(48).fill(0);
  let realNever = 0;
  for (const variant of model.real) {
    if (!(variant.support > 0)) realNever += 1;
    else realShares[Math.min(47, Math.max(0, Math.floor((Math.log10(variant.support) + 6) * 8 + 1e-9)))] += 1;
  }
  const reals = Math.max(1, model.real.length);
  histogram(REAL, `Real variants (${model.real.length.toLocaleString('en')})`, realShares.map((n) => n / reals), realNever / reals, 'real', 'the real variants');
  const allShares = model.tail.slice(0, 48).map((t, k) => t.share - (model.tail[k + 1] ? model.tail[k + 1].share : 0));
  while (allShares.length < 48) allShares.push(0);
  histogram(ALL, `All possible changes (${model.counts.arbitrary.toLocaleString('en')})`, allShares, 1 - model.seenShare, 'arbitrary', 'all possible changes');

  // The axis they share.
  for (let decade = -5; decade <= 0; decade += 1) {
    picture.append(text(xOf(10 ** decade), ALL.bottom + 14, 'tick', `${Number((100 * 10 ** decade).toPrecision(1))}%`, decade === 0 ? 'end' : 'middle'));
  }
  picture.append(text((slot0 + slot1) / 2, ALL.bottom + 14, 'tick', 'never seen'));
  picture.append(text((x0 + x1) / 2, ALL.bottom + 30, 'axis-name', 'Support: the lowest plausible frequency in the best population'));

  if (v) {
    const at = xOf(support);
    for (const plot of [SCORE, REAL, ALL]) picture.append(svg('line', { x1: at, x2: at, y1: plot.top, y2: plot.bottom, class: 'this' }));
    picture.append(svg('circle', { cx: at, cy: yOf(Number(v.score)), r: 4.5, class: 'this-dot' }));
    picture.append(text(Math.min(Math.max(at, 150), W - 70), SCORE.top - 8, 'this-label', `this difference: score ${scoreText(Number(v.score))}`));
  }
  return picture;
}

/** What the marks of the comparison mean. */
function comparisonLegend(v) {
  return el('p', { class: 'legend compare-key' },
    el('span', { class: 'key' }, el('i', { class: 'swatch line-score' }), ' score'),
    el('span', { class: 'key' }, el('i', { class: 'swatch band-score' }), ' its 95% range'),
    el('span', { class: 'key' }, el('i', { class: 'swatch real' }), ' real variants'),
    el('span', { class: 'key' }, el('i', { class: 'swatch' }), ' all possible changes'),
    el('span', { class: 'key' }, el('i', { class: 'line' }), ' 1%'),
    v ? el('span', { class: 'key' }, el('i', { class: 'swatch this-mark' }), ' this difference') : null);
}

/** The part of the evidence panel that shows what a single-base difference
 *  was compared with, and what its score is worth. */
export function comparison(v) {
  if (!given(v.score)) return null;
  const model = scoreModel();
  if (!model) return null;
  const score = Number(v.score);
  const [low, high] = given(v.score_low) ? [Number(v.score_low), Number(v.score_high)] : model.interval(number(v.support));
  const reading = score >= 10 ? 'That is strong evidence for a real variant.' : score >= state.limit && score > 1 ? 'That speaks for a real variant and is above the limit for a supported difference, without being strong evidence.'
    : score > 1 ? 'That speaks for a real variant, the more so the higher the score, but is not strong evidence by itself. "Inconclusive" means exactly that, not "ruled out".'
      : score > 0.1 ? 'That speaks mildly against a real variant without ruling it out.'
        : 'That is what an arbitrary change looks like: the difference needs confirmation from other data.';
  const help = el('span', { class: 'help-text' },
    el('span', { text: model.own
      ? `The real variants are the ${model.real.length} distinct single-base variants of the ${model.alleles} alleles (${plural(model.genes, 'gene')}) that this run accepts as real, namely the database alleles and the candidates whose whole sequence was found on several phased chromosomes. All possible changes are the ${model.counts.arbitrary.toLocaleString('en')} single-base changes that are possible in the genes of the run, which is what an arbitrary error looks like. The line is the score: how much more typical a support is of real variants than of possible changes. It is fitted to all of them at once, so it has no steps.`
      : `The scores of this run come from igpop's built-in curve, because the run has too few established alleles to make its own. The curve was made from ${model.counts.trusted} real variants and ${model.counts.arbitrary.toLocaleString('en')} possible changes in human IGHV genes.` }),
    el('span', { text: 'The three plots share the axis of support. The upper one is the score; the two below show how the real variants and all possible changes spread over support, each as a share of its own group.' }),
    el('span', { text: `The limits for the verdict (${state.limit} for "supported", which you can move above the list of alleles, and 0.1 for "needs confirmation") are conventions, not measurements. Where likelihood ratios are used, 1 to 3 is usually read as weak, 3 to 10 as positive and 10 or more as strong evidence. Because the score is smooth, a difference just below a limit is nearly as well supported as one just above it: read the number, not only the verdict.` }));
  return el('div', { class: 'evidence-part' },
    el('h4', {}, 'Score against support ', hint('How to read this figure', help)),
    el('p', {}, `Score ${scoreText(score)} (95% range ${scoreText(low)} to ${scoreText(high)}). ${reading}`),
    comparisonLegend(v),
    el('div', { class: 'scroll' }, comparisonChart(v)),
    v.score_basis === 'borrowed' ? el('p', { class: 'muted', text: 'An insertion or deletion: its score is read from the same line as a lower bound.' }) : null);
}


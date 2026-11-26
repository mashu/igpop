// The differences of an allele: the table, and the evidence of one
// difference when it is opened.

import { el, explain, hint, plural, svg } from './dom.js?v=BUILD';
import { ancestryLegend, failures, frequencies, populations, stratumWords } from './populations.js?v=BUILD';
import { comparison, scoreModel } from './score.js?v=BUILD';
import { state } from './state.js?v=BUILD';
import { given, number } from './tsv.js?v=BUILD';
import { SITE_QC, badge, chips, pct, percent, scoreText } from './words.js?v=BUILD';

/** The key to the frequency lines, shown with each table of differences. */
function frequencyLegend() {
  return el('div', { class: 'legend' },
    el('span', { class: 'key' }, 'One dot per population that has seen the base, coloured like the bars. '),
    el('span', { class: 'key' }, el('i', { class: 'line' }), ' 1% line'));
}

/** What a score of this size means, in a few words. */
export function scoreWords(score) {
  if (score >= 10) return 'typical of real variants';
  if (score >= 3) return 'speaks for a real variant';
  if (score > 1 / 3) return 'says little either way';
  if (score > 0.1) return 'speaks against';
  return 'typical of arbitrary changes';
}

/** The score explained, with a few points of this run's curve. */
function scoreHint() {
  const model = scoreModel();
  const points = model ? [[0, 'never seen in any dataset'], [1e-4, 'shown to reach 0.01% in a population'], [1e-3, 'shown to reach 0.1%'], [1e-2, 'shown to reach 1%'], [1e-1, 'shown to reach 10%']] : [];
  return el('span', {},
    el('b', { text: 'How to read the score' }),
    el('span', { text: 'The score says how much more typical the frequency of this base is of the variants of real alleles than of an arbitrary base change (what an error would give). It rises smoothly with the lowest frequency the counts are compatible with:' }),
    points.length ? el('span', { class: 'hint-table' }, points.map(([support, words]) => el('span', {},
      el('b', { text: scoreText(model.score(support)) }),
      el('span', { text: `${words}: ${scoreWords(model.score(support))}` })))) : null,
    el('span', { text: `${state.limit} or more counts as supported (the limit you can set above the list), 0.1 or less as needing confirmation. These limits are conventions; around 1 the frequency cannot tell, and 10 or more is usually read as strong evidence.` }));
}

const LOG_MIN = -6; // strip runs from 1 in a million to 100%

const STRIP = { width: 260, height: 38, pad: 14, line: 12 };

// Labelled points of the scale, as shares of the chromosomes examined.
const SCALE = [[1e-6, '0.0001%'], [1e-4, '0.01%'], [1e-2, '1%'], [1, '100%']];

/** How to read a frequency strip. */
function stripHint() {
  return el('span', {},
    el('b', { text: 'How to read the line' }),
    el('span', { text: 'Each mark is one population or dataset in which this exact base was seen. Its place on the line shows how many of the chromosomes examined there carry it.' }),
    el('span', { text: 'The line runs from 0.0001% (1 in a million) at the left end to 100% at the right end. Each tick is a step of a hundred, because real variants range from very rare to very frequent.' }),
    el('span', { text: 'The dashed line marks 1%. A variant to the right of it in some population is called common, which is strong evidence that it is real. No mark at all means the base was never seen.' }));
}

export const xOf = (af) => STRIP.pad + ((Math.max(Math.log10(af), LOG_MIN) - LOG_MIN) / -LOG_MIN) * (STRIP.width - 2 * STRIP.pad);

function strip(v) {
  const marks = frequencies(v);
  const best = marks.reduce((top, m) => (m.af > top.af ? m : top), { af: 0 });
  const label = marks.length
    ? `Seen in ${marks.length} population${marks.length > 1 ? 's' : ''}; highest ${percent(best.af)} in ${best.name}.`
    : 'Not seen in any population.';
  const mid = STRIP.line;
  const picture = svg('svg', { viewBox: `0 0 ${STRIP.width} ${STRIP.height}`, class: 'strip', role: 'img', 'aria-label': label });
  for (const [share, name] of SCALE) {
    const anchor = share === 1e-6 ? 'start' : share === 1 ? 'end' : 'middle';
    const x = share === 1e-6 ? 0 : share === 1 ? STRIP.width : xOf(share);
    picture.append(svg('text', { x, y: STRIP.height - 3, 'text-anchor': anchor, class: 'scale' }, name));
  }
  picture.append(svg('line', { x1: STRIP.pad, x2: STRIP.width - STRIP.pad, y1: mid, y2: mid, class: 'axis' }));
  for (let power = LOG_MIN; power <= 0; power += 2) {
    picture.append(svg('line', { x1: xOf(10 ** power), x2: xOf(10 ** power), y1: mid - 3, y2: mid + 3, class: 'axis' }));
  }
  picture.append(svg('line', { x1: xOf(0.01), x2: xOf(0.01), y1: 1, y2: 2 * mid - 1, class: 'one-percent' }));
  for (const m of marks) {
    const x = xOf(m.af);
    const shape = svg('circle', { cx: x, cy: mid, r: 3.5, class: `seen a-${m.ancestry}` });
    shape.append(svg('title', {}, `${m.name}: ${percent(m.af)}`));
    picture.append(shape);
  }
  return el('div', { class: 'strip-cell' }, picture, el('small', { text: label }));
}

/** Why a difference has its verdict, as a few facts. */
function whyBrief(v) {
  const support = number(v.support);
  const fact = (name, value) => el('div', {}, el('dt', { text: name }), el('dd', {}, value));
  return el('dl', { class: 'why' },
    fact('Best evidence', support > 0 && given(v.support_from) ? stratumWords(v.build, v.support_from).replace(/^the /, '') : 'never seen'),
    support > 0 ? fact('Lowest plausible frequency', percent(support)) : null,
    fact('Score', given(v.score) ? `${v.score_basis === 'borrowed' ? 'at least ' : ''}${scoreText(Number(v.score))}` : 'none'),
    fact('Common (1% or more)', v.class === 'common' ? 'yes' : 'no'),
    fact('Verdict', badge(v.verdict)));
}

/** Why a difference has its verdict and score, step by step. */
function whyScore(v) {
  const steps = [];
  const support = number(v.support);
  if (!(support > 0)) steps.push('No population of any dataset has seen this base.');
  else if (given(v.support_from)) steps.push(`The strongest evidence comes from ${stratumWords(v.build, v.support_from)}. With its counts, the frequency there is at least ${percent(support)}: the lower end of the 95% confidence range. This lowest plausible frequency is what igpop calls the support. It is below the observed frequency, the more so the fewer people were examined.`);
  else steps.push(`In its best population the frequency is at least ${percent(support)}: the lower end of the 95% confidence range.`);
  if (given(v.score)) {
    const model = scoreModel();
    const score = Number(v.score);
    const known = given(v.score_low) ? ` How well the real variants pin the curve down at this point: between ${scoreText(Number(v.score_low))} and ${scoreText(Number(v.score_high))}.` : '';
    steps.push(support > 0
      ? `On the curve of this run, a base with this support is ${score >= 1 ? `${scoreText(score)} times more typical of the variants of real alleles than of arbitrary changes` : `${scoreText(1 / score)} times more typical of arbitrary changes than of the variants of real alleles`}: a score of ${scoreText(score)}.${known}`
      : `${model ? `${pct(1 - model.seenShare)} of all possible single-base changes in these genes have never been seen either, but ${model.counts.trusted - model.counts.trusted_seen === 0 ? 'none' : model.counts.trusted - model.counts.trusted_seen} of the ${model.counts.trusted} variants of real alleles. ` : ''}That gives a score of ${scoreText(score)}.${known}`);
    if (given(v.chance) && support > 0) steps.push(`How often does chance do as well? Of all ${model ? model.counts.arbitrary.toLocaleString('en') : ''} single-base changes that are possible in the genes of the run, ${pct(Number(v.chance))} are supported at least as well as this one.`);
    if (v.score_basis === 'borrowed') steps.push('This is an insertion or deletion. The curve is made from single-base changes, so the score is borrowed from it and is a lower bound: random insertions and deletions turn up in population data far more rarely than random base changes, so the same support counts at least as much.');
  } else steps.push(v.type === 'SNV' ? 'No score could be given.' : 'No score: the curve is made from single-base changes and is borrowed for insertions and deletions only as a lower bound. For this support the bound is below 1 and says nothing, so the verdict rests on whether the change was seen at all.');
  if (v.class === 'common') steps.push(given(v.class_basis) && v.class_basis !== 'dbSNP:COMMON'
    ? `Apart from the score, a test shows that the frequency reaches 1% in ${stratumWords(v.build, v.class_basis)}${given(v.class_p) ? ` (p = ${Number(v.class_p).toPrecision(2)} after correction for the number of populations)` : ''}. A common variant is supported whatever its score.`
    : 'dbSNP marks this base as common. A common variant is supported whatever its score.');
  else if (support > 0) steps.push(`No population shows a frequency of 1% or more${given(v.class_p) ? ` (best test: p = ${Number(v.class_p).toPrecision(2)})` : ''}, so the verdict follows the score: ${state.limit} or more is supported, 0.1 or less needs confirmation, in between inconclusive.`);
  return steps;
}

/** Everything behind one difference: where it is, who carries it, and why
 *  it has its score. */
function evidencePanel(v) {
  const groups = populations(v);
  const [carried, without] = bySeen(groups);
  const top = Math.max(1e-12, ...carried.flatMap((g) => g.rows.map(share)));
  const bar = (r) => el('tr', { class: share(r) > 0 ? '' : 'zero' },
    el('th', { scope: 'row', text: r.name }),
    el('td', {}, el('span', { class: `bar-fill a-${r.ancestry}`, style: `--share:${(100 * share(r)) / top}%` })),
    el('td', { class: 'num', text: Number.isFinite(r.ac) && r.an ? `${r.estimated ? 'about ' : ''}${r.ac.toLocaleString('en')} of ${r.an.toLocaleString('en')}` : '' }),
    el('td', { class: 'num', text: percent(share(r)) }));
  const kind = { SNV: 'single-base change', INS: 'insertion', DEL: 'deletion' }[v.type] || v.type;
  return el('div', { class: 'evidence' },
    el('div', { class: 'evidence-part' },
      el('h4', { text: 'Where' }),
      el('p', {}, `Base ${v.aln_pos} of your sequence (counted from its first base${given(v.aln_pos_genome_strand) && v.aln_pos_genome_strand !== v.aln_pos ? `; base ${v.aln_pos_genome_strand} when the allele is read along the forward strand of the genome, because this gene lies on the reverse strand` : ''}), ${v.build || ''} ${v.contig}:${Number(v.genome_pos).toLocaleString('en')} on the genome. `,
        `A ${kind}: the reference genome has ${v.reference_bases}, this allele has ${v.query_variant}`,
        el('span', { class: 'muted', text: ' (written on the forward strand of the genome, as in the population files; the alignment below shows the gene in its own direction).' }))),
    el('div', { class: 'evidence-part' },
      el('h4', {}, 'Who carries it ', hint('How to read the bars', el('span', { class: 'help-text' }, el('span', { text: 'Each person has two chromosomes, so the counts are of chromosomes.' }), el('span', { text: 'Bars are drawn relative to the population where the base is most frequent. Their colour is the continental group of the population.' })))),
      carried.map((g) => el('table', { class: 'people' },
        el('caption', {}, el('b', { text: g.source }), g.total && g.total.an > 0 ? `: ${g.total.ac.toLocaleString('en')} of ${g.total.an.toLocaleString('en')} chromosomes carry it` : '',
          g.note ? el('small', { class: 'muted', text: ` ${g.note}` }) : null),
        el('thead', { class: 'sr' }, el('tr', {}, ...['Population', 'Share', 'Chromosomes with it', 'Frequency'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, g.rows.sort((a, b) => share(b) - share(a)).map(bar)))),
      without.length ? el('p', { class: 'not-seen' }, el('b', { text: carried.length ? 'Looked for and not seen in: ' : 'Looked for and seen nowhere: ' }), without.map(examined).join('; '), '.') : null,
      groups.length ? null : el('p', { class: 'muted', text: 'No dataset has a record of this base at this position.' }),
      carried.length ? ancestryLegend() : null),
    el('div', { class: 'evidence-part' },
      el('h4', {}, 'Why this verdict ', hint('The reasoning step by step', el('span', { class: 'help-text' }, whyScore(v).map((step, i) => el('span', { text: `${i + 1}. ${step}` }))))),
      whyBrief(v)),
    comparison(v));
}

/** A difference in words: from what to what. */
export function changeWords(v) {
  const [from, to] = [v.reference_bases, v.query_variant];
  if (v.type === 'SNV') return from ? `base change ${from} to ${to}` : `base change to ${to}`;
  // Insertions and deletions are written with the base before them.
  if (v.type === 'DEL') return from ? `deletion of ${from.slice(to.length)}` : `deletion, leaving ${to}`;
  if (v.type === 'INS') return from ? `insertion of ${to.slice(from.length)}` : `insertion, ${to}`;
  return `${v.type} ${to}`;
}

/** The frequency of one row of population data. */
const share = (r) => r.af ?? (r.an ? r.ac / r.an : 0);

/** The datasets of a difference in two lists: those that have seen it, and
 *  those that looked for it and have not. */
function bySeen(groups) {
  const seen = (g) => g.rows.some((r) => share(r) > 0);
  return [groups.filter(seen), groups.filter((g) => !seen(g))];
}

/** A dataset that has not seen a difference, with the chromosomes it examined. */
const examined = (g) => (g.total && g.total.an > 0 ? `${g.source} (${g.total.an.toLocaleString('en')} chromosomes)` : g.source);

/** The population data of one difference, briefly, for a popup: per
 *  source the chromosomes that carry it, and the populations where it is
 *  most frequent as bars on the scale of the frequency line. */
function peoplePopup(v) {
  const [carried, without] = bySeen(populations(v));
  const support = number(v.support);
  // The same scale as the frequency line: one in a million to 100%.
  const width = (af) => `${Math.max(2, Math.min(100, ((Math.log10(af) - LOG_MIN) / -LOG_MIN) * 100))}%`;
  const blocks = carried.map((g) => {
    const seen = g.rows.filter((r) => share(r) > 0).sort((x, y) => share(y) - share(x));
    const top = seen.slice(0, 3);
    return el('span', { class: 'tip-source' },
      el('span', { class: 'tip-head' }, el('b', { text: g.source }),
        g.total && g.total.an > 0 ? el('span', { text: `${g.total.ac.toLocaleString('en')} of ${g.total.an.toLocaleString('en')}` }) : null),
      el('span', { class: 'tip-rows' }, top.flatMap((r) => [
        el('span', { class: 'tip-name', text: r.name }),
        el('span', { class: 'tip-track' }, el('i', { class: `a-${r.ancestry}`, style: `width:${width(share(r))}` })),
        el('span', { class: 'tip-count', text: Number.isFinite(r.ac) && r.an ? `${r.estimated ? '~' : ''}${r.ac.toLocaleString('en')} of ${r.an.toLocaleString('en')}` : '' }),
        el('span', { class: 'tip-share', text: percent(share(r)) }),
      ])),
      seen.length > top.length ? el('span', { class: 'tip-none', text: `and ${plural(seen.length - top.length, 'more population')}` }) : null);
  });
  if (without.length) blocks.push(el('span', { class: 'tip-none', text: `${carried.length ? 'Not seen in' : 'Looked for and seen nowhere'}: ${without.map(examined).join('; ')}.` }));
  return el('span', { class: 'people-tip' },
    blocks.length ? blocks : el('span', { text: 'No dataset has a record of this base at this position.' }),
    carried.length ? el('span', { class: 'tip-none', text: 'Counts are chromosomes. Bars run from 0.0001% to 100% in steps of ten; the line marks 1%.' }) : null,
    support > 0 || given(v.score) ? el('span', { class: 'tip-foot' },
      support > 0 ? el('span', {}, 'Lowest plausible frequency ', el('b', { text: percent(support) })) : el('span', { text: 'Never seen' }),
      given(v.score) ? el('span', {}, 'Score ', el('b', { text: scoreText(Number(v.score)) })) : null) : null);
}

export function variantTable(variants, opened) {
  const rsLink = (v) => (given(v.rsid)
    ? el('a', { href: `https://www.ncbi.nlm.nih.gov/snp/${v.rsid}`, target: '_blank', rel: 'noopener', text: v.rsid })
    : el('span', { class: 'muted', text: 'no rsID' }));
  const change = changeWords;
  const rows = variants.flatMap((v) => {
    const detail = el('tr', { class: 'evidence-row', hidden: true }, el('td', { colspan: '7' }));
    const button = el('button', { type: 'button', class: 'link', 'aria-expanded': 'false', text: 'show', 'aria-label': 'Show populations and reasons' });
    const toggle = (open) => {
      detail.hidden = !open;
      button.setAttribute('aria-expanded', String(open));
      button.textContent = open ? 'hide' : 'show';
      if (open && !detail.firstElementChild.childElementCount) detail.firstElementChild.append(evidencePanel(v));
    };
    button.addEventListener('click', () => toggle(detail.hidden));
    // The alignment opens the evidence of the difference under the pointer.
    if (opened) opened.set(Number(v.genome_pos), () => { toggle(true); detail.scrollIntoView({ block: 'center', behavior: 'smooth' }); });
    return [el('tr', {},
      el('th', { scope: 'row' }, rsLink(v),
        explain(el('span', { class: 'hint-value plain', tabindex: '0' }, el('small', { text: change(v) }), el('small', { text: `${v.contig}:${Number(v.genome_pos).toLocaleString('en')}` })), peoplePopup(v))),
      el('td', {}, v.aln_pos, given(v.aln_pos_genome_strand) && v.aln_pos_genome_strand !== v.aln_pos
        ? el('small', { class: 'other-end', title: 'The same base counted along the forward strand of the genome, that is, from the other end of the allele. Tools that work on the genome strand, such as IgSNPer, report this number.', text: `${v.aln_pos_genome_strand} on genome strand` }) : null),
      el('td', {}, badge(v.verdict)),
      el('td', {}, chips(v.filters), siteQuality(v)),
      scoreCell(v),
      el('td', {}, explain(el('span', { class: 'hint-value plain', tabindex: '0' }, strip(v)), peoplePopup(v))),
      el('td', {}, button)), detail];
  });
  return el('div', { class: 'scroll' }, el('table', { class: 'variants' },
    el('caption', { class: 'sr', text: 'Differences from the reference genome' }),
    el('thead', {}, el('tr', {},
      ...['Variant', 'Position in allele', 'Verdict by frequency', 'Why'].map((h) => el('th', { scope: 'col', text: h })),
      el('th', { scope: 'col' }, 'Score ', hint('What the score means', scoreHint())),
      el('th', { scope: 'col' }, 'How common in people ', hint('How to read the frequency line', stripHint())),
      el('th', { scope: 'col', text: 'Details' }))),
    el('tbody', {}, rows)), frequencyLegend());
}

/** A note under the reasons when the position fails a dataset's site checks. */
function siteQuality(v) {
  const failed = failures(v);
  if (!failed.length) return null;
  const words = (item) => (SITE_QC.find(([pattern]) => pattern.test(item)) || [null, item])[1];
  return el('span', { class: 'qc' }, explain(
    el('span', { class: 'chip chip-qc hint-value', tabindex: '0', text: 'less reliable position' }),
    el('span', {},
      el('b', { text: 'The genotypes at this position are less reliable' }),
      ...failed.map(({ source, failed: items }) => el('span', { class: 'hint-list' }, el('span', { text: `${source}:` }), items.map((item) => el('span', { text: `\u2022 ${words(item)}` })))),
      el('span', { text: 'This is common in antibody genes, which exist in several similar copies. It does not mean the variant is wrong. It means its frequency in that dataset is less certain.' }))));
}

function scoreCell(v) {
  if (!given(v.score)) return el('td', { class: 'num', text: '–', title: v.type === 'SNV' ? 'No score.' : 'No score: for an insertion or deletion seen this rarely the borrowed score says nothing.' });
  const score = Number(v.score);
  const borrowed = v.score_basis === 'borrowed';
  return el('td', { class: 'num' }, explain(el('span', { class: 'hint-value', tabindex: '0', text: `${borrowed ? '\u2265 ' : ''}${scoreText(score)}` }), el('span', {},
    el('b', { text: `${borrowed ? 'At least ' : ''}${scoreText(score)}: ${scoreWords(score)}` }),
    borrowed ? el('span', { text: 'Borrowed from the curve of single-base changes, as a lower bound: random insertions and deletions are found in population data far more rarely than random base changes.' }) : null,
    el('span', { text: `A base seen this often is ${score >= 1 ? `${scoreText(score)} times more typical of real variants than of arbitrary changes` : `${scoreText(1 / score)} times more typical of arbitrary changes than of real variants`}.` }),
    given(v.score_low) ? el('span', { text: `The calibration knows this score to between ${scoreText(Number(v.score_low))} and ${scoreText(Number(v.score_high))}.` }) : null,
    given(v.chance) && number(v.support) > 0 ? el('span', { text: `${pct(Number(v.chance))} of all possible changes are supported at least as well.` }) : null)));
}

/** A part of an opened allele: a heading that says what it is about, a
 *  line that says what question it answers, and its content. */
export function detailSection(title, question, ...content) {
  return el('section', { class: 'detail-section' },
    el('header', {}, el('h4', { text: title }), question ? el('p', { class: 'muted', text: question }) : null),
    ...content.filter(Boolean));
}


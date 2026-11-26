// The tabs beside the allele list: data used, score and checks, report,
// files.

import { $, el, plural } from './dom.js?v=BUILD';
import { datasetsOf, populations } from './populations.js?v=BUILD';
import { comparisonChart, scoreModel } from './score.js?v=BUILD';
import { isCandidate, state } from './state.js?v=BUILD';
import { given, number, parseTSV } from './tsv.js?v=BUILD';
import { FILE_NOTES, OUTPUTS, badge, pct, percent, scoreText } from './words.js?v=BUILD';

/** What each data source contributed to this run, from the result tables:
 *  how many of the differences it has seen and for how many it settled
 *  that they are common. */
export function renderData() {
  const view = $('view-data');
  view.replaceChildren();
  // One entry per distinct difference; alleles often share them.
  const distinct = new Map();
  for (const rows of state.variants.values()) for (const v of rows) distinct.set(`${v.build}|${v.contig}|${v.genome_pos}|${v.query_variant}`, v);
  const all = [...distinct.values()];
  if (!all.length) {
    view.append(el('p', { class: 'muted', text: 'No differences from the reference in these results, so no population data were needed.' }));
    return;
  }
  const builds = [...new Set(all.map((v) => v.build).filter(given))];
  const sections = builds.map((build) => {
    const mine = all.filter((v) => v.build === build);
    const assessed = mine.filter((v) => v.class !== 'not_assessed');
    // Per dataset: how many chromosomes it holds, how many of the
    // differences it has seen, and how many it shows to be common.
    const seenBy = mine.map((v) => populations(v));
    const rows = datasetsOf(build).map((from) => ({
      from,
      chromosomes: [...from.groups.values()].reduce((sum, group) => sum + (group.used ? group.chromosomes : 0), 0),
      seen: seenBy.filter((groups) => groups.some((g) => g.source === from.title && g.rows.some((r) => r.ac > 0))).length,
      common: mine.filter((v) => v.class === 'common' && given(v.class_basis) && v.class_basis.split(':')[0] === from.name).length,
    })).filter((r) => r.chromosomes > 0);
    const never = assessed.filter((v) => !given(v.evidence)).length;
    const alleles = state.alleles.filter((a) => a.build === build && Number(a.n_variants) > 0);
    const looked = alleles.filter((a) => ['found_exact', 'found_with_extra', 'not_found'].includes(a.hap_status));
    const found = looked.filter((a) => a.hap_status === 'found_exact');
    const panel = Math.max(0, ...alleles.map((a) => number(a.hap_total)).filter((x) => x > 0));
    return el('section', { class: 'data-build' },
      el('h3', { text: builds.length > 1 ? `Genes on ${build}: ${plural(mine.length, 'distinct difference')} in ${plural(alleles.length, 'allele')}` : `${plural(mine.length, 'distinct difference')} from the reference, in ${plural(alleles.length, 'allele')}` }),
      build !== state.build ? el('p', { class: 'muted', text: `These genes are missing from ${state.build} and were judged with the older, smaller data of ${build}.` }) : null,
      el('div', { class: 'scroll' }, el('table', { class: 'data-table' },
        el('caption', { text: 'Frequencies: has anyone been seen with the same base?' }),
        el('thead', {}, el('tr', {}, ...['Dataset', 'Size', `Differences it has seen (of ${mine.length})`, 'Shown by it to be common'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, rows.map((r) => el('tr', {},
          el('th', { scope: 'row', text: r.from.title }),
          el('td', { class: 'muted', text: r.from.estimated ? `${r.from.groups.size} studies, sizes estimated` : `${plural([...r.from.groups.values()].filter((g) => g.used).length, 'group')}, up to ${r.chromosomes.toLocaleString('en')} chromosomes` }),
          el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * r.seen) / mine.length}%` }), ` ${r.seen}`),
          el('td', { class: 'num', text: String(r.common) })))),
        el('tfoot', {}, el('tr', {}, el('td', { colspan: '4', class: 'muted', text: `${never} of the differences ${never === 1 ? 'was' : 'were'} seen in no dataset${mine.length - assessed.length ? `; ${mine.length - assessed.length} could not be looked up at all` : ''}. A difference counts once per source; "shown to be common" names the source whose population proved a frequency of 1% or more.` }))))),
      panel ? el('p', {}, el('b', { text: 'Whole sequences: ' }), `the phased panel has ${panel.toLocaleString('en')} chromosomes. ${looked.length} of ${alleles.length} alleles could be looked up completely, and ${found.length} of those were found with their exact sequence.`) : null);
  });
  view.append(el('p', { text: 'Which data the verdicts rest on, and how much each dataset had to say about the alleles of this run.' }), ...sections);

  // Every population of the run, by how many of the differences it has seen.
  for (const build of builds) {
    const mine = all.filter((v) => v.build === build).map((v) => populations(v));
    const groups = datasetsOf(build).flatMap((from) => [...from.groups.values()].filter((group) => group.used).map((group) => ({
      from,
      group,
      seen: mine.filter((sources) => sources.some((g) => g.source === from.title && g.rows.some((r) => r.name === group.title && r.ac > 0))).length,
    })));
    if (!groups.length) continue;
    const top = Math.max(1, ...groups.map((g) => g.seen));
    view.append(el('details', { class: 'fold' }, el('summary', { text: `The ${groups.length} populations${builds.length > 1 ? ` of ${build}` : ''}, one by one` }),
      el('p', { class: 'muted', text: 'Every group of people the verdicts can rest on. A larger group can demonstrate rarer variants. dbSNP passes on frequencies of its studies without the number of people; igpop estimates each size from the smallest frequency the study reports.' }),
      el('div', { class: 'scroll' }, el('table', { class: 'data-table' },
        el('thead', {}, el('tr', {}, ...['Dataset', 'Population', 'Chromosomes', 'Differences it has seen'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, groups.map(({ from, group, seen }) => el('tr', {},
          el('td', { class: 'muted', text: from.title }),
          el('th', { scope: 'row' }, el('i', { class: `swatch a-${group.ancestry}` }), ` ${group.title}`),
          el('td', { class: 'num', text: `${from.estimated ? 'about ' : ''}${group.chromosomes.toLocaleString('en')}` }),
          el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * seen) / top}%` }), ` ${seen}`))))))));
  }
}

export function renderMethod() {
  const view = $('view-method');
  view.replaceChildren();
  const model = scoreModel();
  if (!model) {
    view.append(el('p', { class: 'muted', text: 'No score in these results.' }));
    return;
  }
  const source = model.own ? 'calibrated on this run' : model.levels.some((r) => r.score_source === 'file') ? 'a curve from an earlier run' : 'the built-in human IGHV curve (this run has too few trusted variants of its own)';
  const points = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1];
  view.append(
    el('h3', { text: 'How the score was made' }),
    el('p', { text: 'The score of a difference is a likelihood ratio: how many times more typical its population support is of the variants of real alleles than of arbitrary single-base changes in the same genes. Above 1 the data speak for a real variant, below 1 against; near 1 they say nothing either way. Real alleles are the named trusted ones and those confirmed as haplotypes (two or more differences found together on phased chromosomes).' }),
    el('p', { text: `The score is a smooth curve of support, fitted to all ${model.counts.trusted} real variants and all ${model.counts.arbitrary.toLocaleString('en')} possible changes at once. Two differences with nearly the same support therefore get nearly the same score. A change that no dataset has seen is a case of its own, with a score of ${scoreText(model.never)}.` }),
    el('p', {}, 'Scores in use: ', el('b', { text: source }), '.'),
    el('div', { class: 'scroll' }, comparisonChart(null)),
    el('div', { class: 'scroll' }, el('table', { class: 'data-table' },
      el('caption', { text: 'Points on the curve' }),
      el('thead', {}, el('tr', {}, ...['Support: lowest plausible frequency in the best population', 'Score', 'How well the score is known', 'Possible changes supported at least as well'].map((h) => el('th', { scope: 'col', text: h })))),
      el('tbody', {}, points.map((support) => el('tr', {},
        el('th', { scope: 'row', text: support > 0 ? percent(support) : 'never seen' }),
        el('td', { class: 'num', text: scoreText(model.score(support)) }),
        el('td', { class: 'num', text: support > 0 ? `${scoreText(model.interval(support)[0])} to ${scoreText(model.interval(support)[1])}` : '' }),
        el('td', { class: 'num', text: pct(model.chance(support)) })))))),
    el('p', { class: 'muted' }, 'By frequency, a difference is ', badge('supported'), ` from a score of ${state.limit} or when it is shown to reach 1% in a population, `, badge('needs_confirmation'), ' up to 0.1, and ', badge('inconclusive'), ' in between. The two limits are conventions for strong evidence, not measured.'));
  const hapText = state.files['haplotype_calibration.tsv'];
  if (hapText) {
    const hapRows = parseTSV(hapText).filter((r) => ['1', '2', '3', '5', '10', '20', '50'].includes(r.chromosomes));
    view.append(
      el('h3', { text: 'How much is a found haplotype worth?' }),
      el('p', { text: 'Every founder chromosome of the phased panel was changed one base at a time, and the changed sequence looked up among the chromosomes. This is what an error on a real allele would look like. An allele is supported by its haplotype when its exact sequence is found more often than all but 5% of such changes, and on at least two chromosomes.' }),
      el('div', { class: 'scroll' }, el('table', { class: 'calibration' },
        el('thead', {}, el('tr', {}, ...['Found on at least', 'Share of arbitrary changes', 'Counts as support'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, hapRows.map((r) => el('tr', {},
          el('th', { scope: 'row', text: `${r.chromosomes} chromosome${r.chromosomes === '1' ? '' : 's'}` }),
          shareCell(r.coincidence, null, 'arbitrary'),
          el('td', { text: r.supports })))))));
  }
  const trusted = state.alleles.filter((a) => !isCandidate(a) && Number(a.n_variants) > 0);
  if (trusted.length) {
    const supported = trusted.filter((a) => a.verdict === 'supported').length;
    view.append(el('h3', { text: 'Check: trusted alleles, judged like candidates' }),
      el('p', { text: `${supported} of ${trusted.length} trusted alleles that differ from the reference come out as supported. Most should; if they do not, the coordinates or the choice of trusted alleles need a look before the candidates are read.` }));
  }
}

function shareCell(fraction, count, kind) {
  const share = given(fraction) ? Number(fraction) : 0;
  const text = count == null ? ` ${(100 * share).toFixed(2)}%` : ` ${(100 * share).toFixed(1)}% (${count})`;
  return el('td', {}, el('span', { class: `bar-fill ${kind}`, style: `--share:${100 * share}%` }), el('span', { class: 'num', text }));
}

export function inline(text) {
  const out = document.createDocumentFragment();
  for (const piece of text.split(/(`[^`]+`|\*\*[^*]+\*\*)/)) {
    if (piece.startsWith('`')) out.append(el('code', { text: piece.slice(1, -1) }));
    else if (piece.startsWith('**')) out.append(el('b', { text: piece.slice(2, -2) }));
    else if (piece) out.append(piece);
  }
  return out;
}

export function renderReport() {
  const view = $('view-report');
  view.replaceChildren();
  const text = state.files['report.md'];
  if (!text) {
    view.append(el('p', { class: 'muted', text: 'No report in these results.' }));
    return;
  }
  const lines = text.split('\n');
  for (let i = 0; i < lines.length; i += 1) {
    const line = lines[i];
    const heading = line.match(/^(#{1,3}) (.*)$/);
    if (heading) {
      view.append(el(`h${heading[1].length + 2}`, {}, inline(heading[2])));
    } else if (line.startsWith('|')) {
      const table = [];
      for (; i < lines.length && lines[i].startsWith('|'); i += 1) table.push(lines[i].slice(1, -1).split('|').map((c) => c.trim()));
      i -= 1;
      const [head, , ...body] = table;
      view.append(el('div', { class: 'scroll' }, el('table', {},
        el('thead', {}, el('tr', {}, head.map((c) => el('th', { scope: 'col' }, inline(c))))),
        el('tbody', {}, body.map((row) => el('tr', {}, row.map((c) => el('td', {}, inline(c)))))))));
    } else if (line.startsWith('- ')) {
      const items = [];
      for (; i < lines.length && lines[i].startsWith('- '); i += 1) items.push(lines[i].slice(2));
      i -= 1;
      view.append(el('ul', {}, items.map((item) => el('li', {}, inline(item)))));
    } else if (line.trim()) {
      view.append(el('p', {}, inline(line)));
    }
  }
}

export function renderFiles() {
  const view = $('view-files');
  view.replaceChildren(el('p', { text: 'The same files the command-line program writes. Tables are tab-separated and open in Excel, R or Python.' }));
  const list = el('ul', { class: 'files' });
  for (const name of OUTPUTS) {
    const content = state.files[name];
    if (content == null) continue;
    const url = URL.createObjectURL(new Blob([content], { type: 'text/plain' }));
    list.append(el('li', {},
      el('a', { href: url, download: `${state.prefix}.${name}`, class: 'file' }, el('b', { text: `${state.prefix}.${name}` }), el('span', { text: FILE_NOTES[name] }), el('small', { text: `${Math.max(1, Math.round(content.length / 1024))} kB` }))));
  }
  view.append(list);
}


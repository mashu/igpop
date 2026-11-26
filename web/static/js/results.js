// The results page: the counts, the filters and the list of alleles.

import { alignmentPanel } from './alignment.js?v=BUILD';
import { $, el, hint, plural, selectTab } from './dom.js?v=BUILD';
import { haplotypePanel } from './haplotype.js?v=BUILD';
import { CONFIRM, rememberVerdicts, renderLimit } from './limit.js?v=BUILD';
import { renderData, renderFiles, renderMethod, renderReport } from './pages.js?v=BUILD';
import { loadPopulations } from './populations.js?v=BUILD';
import { isCandidate, on, show, state } from './state.js?v=BUILD';
import { given, number, parseAlignments, parseTSV } from './tsv.js?v=BUILD';
import { changeWords, detailSection, scoreWords, variantTable } from './variants.js?v=BUILD';
import { HAPLOTYPES, PLACEMENTS, ROUTES, VERDICTS, badge, chips, idWords, scoreText } from './words.js?v=BUILD';

export function display(files, warnings, summary = []) {
  state.files = files;
  state.alleles = parseTSV(files['alleles.tsv']);
  state.placements = null;
  state.calibration = null;
  state.model = undefined;
  // The build most alleles were judged on is the build of the run.
  const builds = new Map();
  for (const a of state.alleles) if (given(a.build)) builds.set(a.build, (builds.get(a.build) || 0) + 1);
  state.build = [...builds].sort((x, y) => y[1] - x[1]).map(([build]) => build)[0] || '';
  state.variants = new Map();
  for (const row of files['variants.tsv'] ? parseTSV(files['variants.tsv']) : []) {
    if (!state.variants.has(row.allele_name)) state.variants.set(row.allele_name, []);
    state.variants.get(row.allele_name).push(row);
  }
  loadPopulations(files);
  state.alignments = parseAlignments(files['alignments.txt']);
  rememberVerdicts(files);
  state.verdict = '';
  state.marker = $('opt-marker').value || '_S';

  const candidates = state.alleles.filter(isCandidate).length;
  $('results-title').textContent = `Results for ${state.prefix}`;
  const traffic = summary.find((line) => line.includes("browser's store"));
  $('results-sub').textContent = `${state.alleles.length} alleles: ${candidates} candidates and ${state.alleles.length - candidates} trusted alleles.${traffic ? ` ${traffic}.` : ''}`;
  renderTiles();
  renderWarnings(warnings);
  renderFilters();
  renderLimit();
  renderAlleles();
  renderData();
  renderMethod();
  renderReport();
  renderFiles();
  selectTab('tab-alleles');
  show('results');
  $('results-title').focus();
}

/** The weakest score of an allele as a number on a continuous colour
 *  scale: from the colour of "needs confirmation" at 0.01 through grey at 1
 *  to the colour of "supported" at 100. */
function scorePill(allele) {
  if (!given(allele.min_score)) return null;
  const score = Number(allele.min_score);
  const t = Math.max(-1, Math.min(1, Math.log10(score) / 2));
  const tone = `color-mix(in srgb, var(${t >= 0 ? '--supported' : '--confirm'}) ${Math.round(100 * Math.abs(t))}%, var(--grey))`;
  return el('span', { class: 'score-pill', style: `--tone:${tone}`, title: `Weakest score among the differences of this allele: ${scoreText(score)} (${scoreWords(score)}). An allele is only as well supported by frequency as its weakest difference.` },
    el('i', { 'aria-hidden': 'true' }), `score ${scoreText(score)}`);
}

export function renderTiles() {
  const tiles = $('tiles');
  tiles.replaceChildren();
  const total = state.alleles.length;
  for (const [verdict, info] of Object.entries(VERDICTS)) {
    const rows = state.alleles.filter((a) => a.verdict === verdict);
    if (!rows.length && state.verdict !== verdict) continue;
    const candidates = rows.filter(isCandidate).length;
    tiles.append(el('button', {
      type: 'button', class: `tile v-${verdict}`, 'aria-pressed': String(state.verdict === verdict), 'data-verdict': verdict, title: info.help,
      onclick: (event) => {
        state.verdict = state.verdict === verdict ? '' : verdict;
        for (const tile of tiles.children) tile.setAttribute('aria-pressed', String(tile.dataset.verdict === state.verdict));
        renderAlleles();
        event.currentTarget.blur();
      },
    },
    el('span', { class: 'tile-count', text: rows.length }),
    el('span', { class: 'tile-label' }, el('i', { 'aria-hidden': 'true', text: info.icon }), info.label),
    el('span', { class: 'tile-sub', text: `${candidates} candidates · ${rows.length - candidates} trusted` }),
    el('span', { class: 'tile-bar', style: `--share:${total ? (100 * rows.length) / total : 0}%` })));
  }
}

function renderWarnings(warnings) {
  const box = $('warnings');
  box.replaceChildren();
  if (!warnings || !warnings.length) return;
  box.append(el('details', { class: 'notice' },
    el('summary', { text: `${warnings.length} thing${warnings.length > 1 ? 's' : ''} to know about this run` }),
    el('ul', {}, warnings.map((w) => {
      // The warning about placement leads to the alleles it is about.
      const unclean = /not placed cleanly/.test(w) ? Object.keys(PLACEMENTS).filter((p) => p !== 'ok' && state.alleles.some((a) => a.locus_status === p)) : [];
      return el('li', {}, w.replace(/`/g, '').replace('; see section 3.', '.'),
        unclean.flatMap((p) => [' ', el('button', { type: 'button', class: 'link', text: `Show: ${PLACEMENTS[p].label.toLowerCase()}`, onclick: () => showPlacement(p) })]));
    }))));
}

function renderFilters() {
  const genes = [...new Set(state.alleles.map((a) => a.gene))].sort((a, b) => a.localeCompare(b, 'en', { numeric: true }));
  $('f-gene').replaceChildren(el('option', { value: '', text: 'all genes' }), ...genes.map((g) => el('option', { value: g, text: g })));
  const statuses = [...new Set(state.alleles.map((a) => a.hap_status))].filter((s) => HAPLOTYPES[s]);
  $('f-hap').replaceChildren(el('option', { value: '', text: 'any' }), ...statuses.map((s) => el('option', { value: s, text: HAPLOTYPES[s].label })));
  const places = Object.keys(PLACEMENTS).filter((p) => state.alleles.some((a) => a.locus_status === p));
  $('f-place').replaceChildren(el('option', { value: '', text: 'any' }),
    ...places.map((p) => el('option', { value: p, text: `${PLACEMENTS[p].label} (${state.alleles.filter((a) => a.locus_status === p).length})` })));
  $('f-search').value = '';
  $('f-kind').value = '';
}

/** Show the alleles with one kind of placement, in the allele list. */
function showPlacement(status) {
  state.verdict = '';
  renderTiles();
  for (const id of ['f-search', 'f-gene', 'f-kind', 'f-hap']) $(id).value = '';
  $('f-place').value = status;
  $('fold-filters').open = true;
  selectTab('tab-alleles');
  renderAlleles();
  $('f-place').scrollIntoView({ block: 'center', behavior: 'smooth' });
}

/** What the folded controls are set to, so that nothing set is hidden. */
function renderFolds() {
  const active = ['f-kind', 'f-hap', 'f-place'].filter((id) => $(id).value).length;
  $('fold-filters-state').textContent = active ? `: ${active} set` : '';
  $('fold-reading-state').textContent = `: supported from a score of ${state.limit}`;
}

/** How many differences of an allele are not supported by themselves,
 *  as the group it falls in: 0 to 4, or 5 for five and more. Alleles
 *  without differences or without a place have none. */
const DOUBT_GROUPS = ['0', '1', '2', '3', '4', '5 or more'];

function doubtOf(allele) {
  const variants = state.variants.get(allele.allele_name) || [];
  if (!variants.length) return -1;
  return Math.min(5, variants.filter((v) => v.verdict !== 'supported').length);
}

/** The alleles by their number of unsupported differences: one bar per
 *  number, split by verdict. A bar is a switch that filters the list. */
function renderDoubt() {
  const shown = ['supported', 'inconclusive', 'needs_confirmation', 'not_assessed'];
  const counts = DOUBT_GROUPS.map(() => Object.fromEntries(shown.map((verdict) => [verdict, 0])));
  for (const allele of state.alleles) {
    const group = doubtOf(allele);
    if (group >= 0 && shown.includes(allele.verdict)) counts[group][allele.verdict] += 1;
  }
  const totals = counts.map((c) => shown.reduce((sum, verdict) => sum + c[verdict], 0));
  const most = Math.max(1, ...totals);
  const all = totals.reduce((sum, n) => sum + n, 0);
  $('fold-doubt').hidden = !all;
  for (const group of [...state.doubt]) if (!totals[group]) state.doubt.delete(group);
  $('fold-doubt-state').textContent = state.doubt.size ? `: showing ${[...state.doubt].sort().map((g) => DOUBT_GROUPS[g]).join(', ')}` : '';
  const bars = DOUBT_GROUPS.map((name, group) => {
    const parts = shown.filter((verdict) => counts[group][verdict]).map((verdict) => `${counts[group][verdict]} ${VERDICTS[verdict].label.toLowerCase()}`);
    const button = el('button', { type: 'button', class: 'doubt-bar', 'aria-pressed': String(state.doubt.has(group)), disabled: totals[group] ? null : '',
      title: totals[group] ? parts.join(', ') : 'no allele',
      'aria-label': `${plural(totals[group], 'allele')} with ${name} unsupported ${group === 1 ? 'difference' : 'differences'}${parts.length ? `: ${parts.join(', ')}` : ''}. Show only these.` },
    el('span', { class: 'doubt-count', text: String(totals[group]) }),
    el('span', { class: 'doubt-stack', style: `height:${(100 * totals[group]) / most}%` },
      shown.filter((verdict) => counts[group][verdict]).reverse().map((verdict) => el('i', { class: `v-${verdict}`, style: `flex:${counts[group][verdict]}` }))),
    el('span', { class: 'doubt-name', text: name }));
    button.addEventListener('click', () => {
      if (state.doubt.has(group)) state.doubt.delete(group); else state.doubt.add(group);
      renderAlleles();
    });
    return button;
  });
  $('doubt').replaceChildren(
    el('p', { class: 'legend' },
      ...shown.filter((verdict) => counts.some((c) => c[verdict])).map((verdict) => el('span', { class: 'key' }, el('i', { class: `swatch v-${verdict}` }), ` ${VERDICTS[verdict].label.toLowerCase()}`)),
      hint('How to read this', el('span', { class: 'help-text' },
        el('span', { text: `Each bar counts the alleles with that many differences that are not supported by themselves at the current limit (${state.limit}). The colours are the verdicts of the alleles.` }),
        el('span', { text: 'An allele with unsupported differences can still be supported, when its exact sequence was found on phased chromosomes.' }),
        el('span', { text: 'Select one or more bars to show only those alleles in the list below; select a bar again to release it.' })))),
    el('div', { class: 'doubt-chart', role: 'group', 'aria-label': 'Alleles by number of unsupported differences' }, bars),
    el('p', { class: 'doubt-axis muted', text: 'unsupported differences in the allele' }),
    state.doubt.size ? el('button', { type: 'button', class: 'link', text: 'Show all again', onclick: () => { state.doubt.clear(); renderAlleles(); } }) : '');
}

function visibleAlleles() {
  const query = $('f-search').value.trim().toLowerCase();
  const gene = $('f-gene').value;
  const kind = $('f-kind').value;
  const hap = $('f-hap').value;
  const place = $('f-place').value;
  return state.alleles.filter((a) => {
    if (place && a.locus_status !== place) return false;
    if (state.verdict && a.verdict !== state.verdict) return false;
    if (state.doubt.size && !state.doubt.has(doubtOf(a))) return false;
    if (gene && a.gene !== gene) return false;
    if (kind === 'candidate' && !isCandidate(a)) return false;
    if (kind === 'trusted' && isCandidate(a)) return false;
    if (hap && a.hap_status !== hap) return false;
    if (query && !`${a.allele_name} ${a.gene} ${a.reason}`.toLowerCase().includes(query)) return false;
    return true;
  });
}

export function renderAlleles() {
  renderDoubt();
  const rows = visibleAlleles();
  renderFolds();
  const list = $('alleles');
  list.replaceChildren();
  $('count').textContent = `Showing ${rows.length} of ${state.alleles.length} alleles.`;
  const order = $('f-order').value;
  if (order) {
    // By weakest score: one list; alleles without a score come last.
    const key = (a) => (given(a.min_score) ? Number(a.min_score) : NaN);
    const sorted = [...rows].sort((x, y) => (Number.isNaN(key(x)) ? 1 : 0) - (Number.isNaN(key(y)) ? 1 : 0) || (order === 'low' ? key(x) - key(y) : key(y) - key(x)));
    list.append(el('section', { class: 'gene' }, sorted.map(alleleCard)));
    if (!rows.length) list.append(el('p', { class: 'muted', text: 'No allele matches these filters.' }));
    return;
  }
  let gene = null;
  let group = null;
  for (const allele of rows) {
    if (allele.gene !== gene) {
      gene = allele.gene;
      group = el('section', { class: 'gene' }, el('h3', { text: gene }));
      list.append(group);
    }
    group.append(alleleCard(allele));
  }
  if (!rows.length) list.append(el('p', { class: 'muted', text: 'No allele matches these filters.' }));
}

function hapPill(allele) {
  const info = HAPLOTYPES[allele.hap_status];
  if (!info || allele.hap_status === 'not_checked') return null;
  const exact = number(allele.hap_exact);
  const text = allele.hap_status === 'found_exact' ? `${info.label} on ${exact} chromosome${exact === 1 ? '' : 's'}` : info.label;
  return el('span', { class: `pill h-${allele.hap_status}`, title: info.help, text });
}

function routeTag(allele) {
  const route = ROUTES[allele.supported_by];
  return route && allele.verdict === 'supported' ? el('span', { class: 'route', title: route.help, text: route.label }) : null;
}

function alleleCard(allele) {
  const body = el('div', { class: 'card-body' });
  const card = el('details', { class: `card v-${allele.verdict}` },
    el('summary', {},
      el('span', { class: 'name', text: allele.allele_name }),
      badge(allele.verdict),
      scorePill(allele),
      routeTag(allele),
      isCandidate(allele) ? el('span', { class: 'tag', text: 'candidate' }) : null,
      given(allele.build) && allele.build !== state.build ? el('span', { class: 'tag', title: `This gene is not on ${state.build}: the allele was judged with the data of ${allele.build}.`, text: allele.build }) : null,
      el('span', { class: 'chips' }, allele.verdict === 'supported' ? [] : chips(allele.filters), hapPill(allele),
        allele.locus_status && allele.locus_status !== 'ok' && PLACEMENTS[allele.locus_status]
          ? el('span', { class: 'pill p-unclean', title: PLACEMENTS[allele.locus_status].help, text: PLACEMENTS[allele.locus_status].label }) : null)),
    body);
  // Details are built when first opened: a full V set has hundreds of cards.
  card.addEventListener('toggle', () => {
    if (card.open && !body.childElementCount) body.append(...alleleDetail(allele));
  });
  return card;
}

/** Why a sequence was not placed cleanly, with the numbers behind it. */
function placementNote(allele) {
  if (!allele.locus_status || allele.locus_status === 'ok' || !PLACEMENTS[allele.locus_status]) return null;
  if (!state.placements) state.placements = new Map((state.files['alignments.tsv'] ? parseTSV(state.files['alignments.tsv']) : []).map((r) => [r.allele, r]));
  const at = state.placements.get(allele.allele_name);
  const lines = [el('b', { text: `${PLACEMENTS[allele.locus_status].label}. ` })];
  if (allele.locus_status === 'ambiguous_locus') {
    if (at && given(at.best_other_gene)) {
      const [own, other] = [Number(at.score), Number(at.best_other_score)];
      lines.push(`Compared with every gene of the list, this sequence fits ${at.best_other_gene} ${other > own ? 'better than' : 'exactly as well as'} ${allele.gene}, the gene of its name (alignment score ${other} against ${own}; ${(100 * Number(at.identity)).toFixed(1)}% identical to ${allele.gene}). `);
    } else lines.push(`This sequence fits another gene of the list at least as well as ${allele.gene}. `);
    lines.push(`It was still compared with ${allele.gene} and judged, but read the verdict with care: the two genes are near-identical copies, so a difference seen here may be the normal sequence of the other gene, and the short reads of the population data are often assigned to the wrong copy as well. If the sequence belongs to the other gene, give it that gene's name and run again.`);
  } else if (allele.locus_status === 'low_identity') {
    lines.push(`The sequence is ${at && given(at.identity) ? `only ${(100 * Number(at.identity)).toFixed(1)}% identical` : 'too different'} to ${allele.gene}, so no differences were called. Check that the name is right, or that the sequence is complete and in one piece.`);
  } else lines.push(PLACEMENTS[allele.locus_status].help);
  return el('p', { class: 'placement-note' }, lines);
}

/** Why the allele has this verdict and not the next one: the rule, with
 *  the numbers of its weakest difference. */
function verdictRule(allele, variants) {
  if (!['inconclusive', 'needs_confirmation'].includes(allele.verdict) || allele.filters === 'sequence_not_in_panel') return null;
  const open = variants.filter((v) => v.verdict === allele.verdict);
  if (!open.length) return null;
  const weakest = open.reduce((low, v) => (given(v.score) && (!given(low.score) || Number(v.score) < Number(low.score)) ? v : low), open[0]);
  const name = given(weakest.rsid) ? weakest.rsid : `the ${changeWords(weakest)} at position ${weakest.aln_pos}`;
  const scored = given(weakest.score) ? `${name} scores ${scoreText(Number(weakest.score))}` : `${name} has no score and ${weakest.class === 'absent' ? 'was never seen' : 'was seen too rarely to tell'}`;
  const rule = allele.verdict === 'needs_confirmation'
    ? `A score of ${CONFIRM} or less means the frequency is what an arbitrary base change looks like, so the allele needs confirmation.`
    : `A score between ${CONFIRM} and ${state.limit} means the frequency cannot tell either way, so the allele is inconclusive.`;
  return el('p', { class: 'muted', text: `An allele is only as well supported as its weakest difference: ${scored}. ${given(weakest.score) ? rule : ''} The limits are ${CONFIRM} for "needs confirmation" and ${state.limit} for "supported"; the haplotype check can still support an allele whose differences fall short.` });
}

function alleleDetail(allele) {
  const placement = placementNote(allele);
  const variants = state.variants.get(allele.allele_name) || [];
  // Genome position of each difference -> opens its evidence in the table.
  const opened = new Map();
  const parts = [detailSection('Verdict', null,
    el('p', { class: 'reason' }, badge(allele.verdict), ' ', idWords(allele.reason).replace(/variant\(s\)/g, variants.length === 1 ? 'difference' : 'differences').replace(/\bvariants\b/g, 'differences')),
    verdictRule(allele, variants),
    placement,
    given(allele.notes) ? el('p', { class: 'muted', text: `Keep in mind: ${allele.notes.replace(/difference\(s\)/g, 'differences').replace(' (column site_quality says why)', '')}.` }) : null)];
  if (variants.length) parts.push(detailSection('Differences from the reference genome', 'Is each difference, taken by itself, carried by people?', variantTable(variants, opened)));
  const hap = haplotypePanel(allele);
  if (hap) parts.push(hap);
  const block = state.alignments.get(allele.gene);
  if (block) parts.push(detailSection('Alignment', 'Where do the differences lie, and what do the other alleles of the gene have there?', alignmentPanel(block, allele.allele_name, opened)));
  return parts;
}

/** The list and the filters of the results page answer to typing. */
export function wireResults() {
  for (const id of ['f-search', 'f-gene', 'f-kind', 'f-hap', 'f-place', 'f-order']) $(id).addEventListener('input', renderAlleles);
}

// Another limit gives other verdicts: everything that shows them follows.
on('verdicts', () => {
  renderTiles();
  renderAlleles();
  renderMethod();
  renderFiles();
});

// igpop web app: page logic. The analysis itself runs in worker.js; this
// file starts it, and turns the result tables into the views of the page.
// The tables are exactly the files the command-line program writes, so the
// page can also display those (drop the .tsv files).

const $ = (id) => document.getElementById(id);
const el = (tag, attrs = {}, ...children) => {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value === false || value == null) continue;
    if (key === 'class') node.className = value;
    else if (key === 'text') node.textContent = value;
    else if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? '' : value);
  }
  for (const child of children.flat()) {
    if (child != null) node.append(child);
  }
  return node;
};
const svg = (tag, attrs = {}, ...children) => {
  const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
  node.append(...children);
  return node;
};

// ---- vocabulary: the fixed words of the tables, in plain language --------

const VERDICTS = {
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
  no_coordinates: 'no coordinates for its gene',
  low_identity: 'too different from its gene',
};
const HAPLOTYPES = {
  found_exact: { label: 'Haplotype found', help: 'Chromosomes with exactly this sequence exist in the phased 1000 Genomes panel.' },
  found_with_extra: { label: 'Found with extra variants', help: 'The combination exists, but always together with further variants inside the allele.' },
  not_found: { label: 'Haplotype not found', help: 'Every variant is in the phased panel, but no chromosome carries them all together.' },
  incomplete: { label: 'Partly testable', help: 'Some variants are not sites of the phased panel, so the full combination cannot be looked up.' },
  untestable: { label: 'Not testable', help: 'None of the variants is a site of the phased panel.' },
  not_checked: { label: 'Not checked', help: 'Not placed, or no phased data for this build.' },
};
const POPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'];
const POP_NAMES = { AFR: 'African', AMR: 'American', EAS: 'East Asian', EUR: 'European', SAS: 'South Asian' };
const OUTPUTS = ['report.md', 'final.tsv', 'alleles.tsv', 'variants.tsv', 'alignments.tsv', 'alignments.txt', 'score_calibration.tsv', 'haplotype_calibration.tsv', 'dbsnp_datasets.tsv'];
const ROUTES = {
  frequency: { label: 'by frequency', help: 'Each difference is an established population variant.' },
  haplotype: { label: 'by haplotype', help: 'The exact sequence is established on phased chromosomes; by frequency alone it would not pass.' },
  'frequency+haplotype': { label: 'by frequency and haplotype', help: 'Each difference is an established variant, and the exact sequence is established on phased chromosomes.' },
};
const FILE_NOTES = {
  'report.md': 'Results and checks in words',
  'final.tsv': 'Alleles that passed',
  'alleles.tsv': 'Every allele: verdict, reasons, haplotype check',
  'variants.tsv': 'Every difference with all population evidence',
  'alignments.tsv': 'Where each allele was placed',
  'alignments.txt': 'Alignments to read by eye',
  'score_calibration.tsv': 'The table behind the score',
  'haplotype_calibration.tsv': 'How often an arbitrary change is found on phased chromosomes',
  'dbsnp_datasets.tsv': 'dbSNP datasets and whether each was used',
};

// ---- parsing ---------------------------------------------------------------

function parseTSV(text) {
  const lines = text.split('\n').filter((line) => line.length);
  const header = lines.shift().split('\t');
  return lines.map((line) => {
    const cells = line.split('\t');
    return Object.fromEntries(header.map((name, i) => [name, cells[i] ?? '.']));
  });
}
const given = (cell) => cell != null && cell !== '.' && cell !== '-' && cell !== '';
const number = (cell) => (given(cell) ? Number(cell) : NaN);

/** Blocks of alignments.txt by gene name. */
function parseAlignments(text) {
  const blocks = new Map();
  for (const block of (text || '').split(/\n(?=>)/)) {
    const match = block.match(/^>(\S+)/);
    if (match) blocks.set(match[1], block.trimEnd());
  }
  return blocks;
}

// ---- state ---------------------------------------------------------------

const state = { prefix: 'igpop', marker: '_S', files: {}, alleles: [], variants: new Map(), alignments: new Map(), verdict: '' };
let worker = null;

function show(section) {
  for (const id of ['start', 'progress', 'error', 'results']) $(id).hidden = id !== section;
}

function fail(message) {
  $('error-text').textContent = message;
  show('error');
  $('error').querySelector('button').focus();
}

// ---- starting a run ------------------------------------------------------

function startWorker() {
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = (event) => {
    const message = event.data;
    if (message.type === 'ready') $('version').textContent = `igpop ${message.version}`;
    if ('stored' in message) showStored(message.stored);
    if (message.type === 'log') logLine(message.message);
    else if (message.type === 'error') fail(message.message);
    else if (message.type === 'done') display(message.output.files, message.output.warnings, message.output.summary);
  };
  worker.onerror = (event) => fail(`The analysis module could not be started (${event.message || 'unknown error'}). A recent Chrome, Edge, Firefox or Safari is needed.`);
}

/** How much downloaded data this browser keeps, with a way to delete it. */
function showStored(bytes) {
  const box = $('stored');
  if (bytes == null) {
    box.textContent = 'This browser does not keep downloaded data between visits (or igpop is open in another tab), so each run downloads again.';
    return;
  }
  const size = bytes < 1e6 ? `${Math.round(bytes / 1e3)} kB` : `${(bytes / 1e6).toFixed(0)} MB`;
  box.replaceChildren(
    `Downloaded data kept in this browser: ${bytes ? size : 'none yet'}. `,
    bytes ? el('button', { type: 'button', class: 'link', text: 'Delete it', onclick: () => worker.postMessage({ type: 'clear' }) }) : '',
  );
}

function logLine(message) {
  const text = message.replace(/https?:\/\/\S+\/([^/\s]+)/g, '$1');
  const log = $('log');
  const last = log.lastElementChild;
  if (text.startsWith('  ') && last) {
    last.append(el('small', { text: ` ${text.trim()}` }));
  } else {
    if (last) last.classList.add('done');
    log.append(el('li', { text }));
  }
  log.scrollTop = log.scrollHeight;
}

async function handleFiles(fileList) {
  const files = [...fileList];
  if (!files.length) return;
  const outputs = files.filter((f) => OUTPUTS.some((suffix) => f.name.endsWith(`.${suffix}`)));
  if (outputs.length) {
    // Output files of an earlier run: display them.
    const loaded = {};
    for (const file of outputs) {
      const suffix = OUTPUTS.find((s) => file.name.endsWith(`.${s}`));
      loaded[suffix] = await file.text();
      state.prefix = file.name.slice(0, -suffix.length - 1);
    }
    if (!loaded['alleles.tsv']) return fail('To browse earlier results, include the .alleles.tsv file (and .variants.tsv for the evidence).');
    return display(loaded, []);
  }
  const file = files[0];
  const fasta = await file.text();
  if (!fasta.trimStart().startsWith('>')) return fail(`${file.name} does not look like a FASTA file (it should start with ">").`);
  state.prefix = file.name.replace(/\.[^.]+$/, '');
  $('progress-name').textContent = file.name;
  $('log').replaceChildren();
  show('progress');
  // Ask the browser not to evict the stored data when disk space is short.
  if (navigator.storage && navigator.storage.persist) navigator.storage.persist().catch(() => {});
  worker.postMessage({
    fasta,
    name: file.name,
    build: $('opt-build').value,
    marker: $('opt-marker').value || '_S',
    haplotypes: $('opt-hap').checked,
    paralog: $('opt-paralog').checked,
  });
}

// ---- results ---------------------------------------------------------------

function display(files, warnings, summary = []) {
  state.files = files;
  state.alleles = parseTSV(files['alleles.tsv']);
  state.variants = new Map();
  for (const row of files['variants.tsv'] ? parseTSV(files['variants.tsv']) : []) {
    if (!state.variants.has(row.allele_name)) state.variants.set(row.allele_name, []);
    state.variants.get(row.allele_name).push(row);
  }
  state.alignments = parseAlignments(files['alignments.txt']);
  state.verdict = '';
  state.marker = $('opt-marker').value || '_S';

  const candidates = state.alleles.filter(isCandidate).length;
  $('results-title').textContent = `Results for ${state.prefix}`;
  const traffic = summary.find((line) => line.includes("browser's store"));
  $('results-sub').textContent = `${state.alleles.length} alleles: ${candidates} candidates and ${state.alleles.length - candidates} trusted alleles.${traffic ? ` ${traffic}.` : ''}`;
  renderTiles();
  renderWarnings(warnings);
  renderFilters();
  renderLegend();
  renderAlleles();
  renderMethod();
  renderReport();
  renderFiles();
  selectTab('tab-alleles');
  show('results');
  $('results-title').focus();
}

// Candidates are named by the marker, as in the command-line program.
const isCandidate = (allele) => allele.allele_name.includes(state.marker);

function badge(verdict) {
  const v = VERDICTS[verdict] || { label: verdict, icon: '' };
  return el('span', { class: `badge v-${verdict}`, title: v.help }, el('i', { 'aria-hidden': 'true', text: v.icon }), v.label);
}

function renderTiles() {
  const tiles = $('tiles');
  tiles.replaceChildren();
  const total = state.alleles.length;
  for (const [verdict, info] of Object.entries(VERDICTS)) {
    const rows = state.alleles.filter((a) => a.verdict === verdict);
    if (!rows.length && verdict === 'not_assessed') continue;
    const candidates = rows.filter(isCandidate).length;
    tiles.append(el('button', {
      type: 'button', class: `tile v-${verdict}`, 'aria-pressed': 'false', 'data-verdict': verdict, title: info.help,
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
    el('ul', {}, warnings.map((w) => el('li', { text: w.replace(/`/g, '') })))));
}

function renderFilters() {
  const genes = [...new Set(state.alleles.map((a) => a.gene))].sort((a, b) => a.localeCompare(b, 'en', { numeric: true }));
  $('f-gene').replaceChildren(el('option', { value: '', text: 'all genes' }), ...genes.map((g) => el('option', { value: g, text: g })));
  const statuses = [...new Set(state.alleles.map((a) => a.hap_status))].filter((s) => HAPLOTYPES[s]);
  $('f-hap').replaceChildren(el('option', { value: '', text: 'any' }), ...statuses.map((s) => el('option', { value: s, text: HAPLOTYPES[s].label })));
  $('f-search').value = '';
  $('f-kind').value = '';
}

function renderLegend() {
  $('legend').replaceChildren(
    el('span', { text: 'Frequency marks: ' }),
    el('span', { class: 'key' }, mark('kgp', 6, 6), ' 1000 Genomes population'),
    el('span', { class: 'key' }, mark('gnomad', 6, 6), ' gnomAD group'),
    el('span', { class: 'key' }, mark('dbsnp', 6, 6), ' dataset in dbSNP'),
    el('span', { class: 'key' }, el('i', { class: 'line' }), ' 1% line'),
  );
}

function visibleAlleles() {
  const query = $('f-search').value.trim().toLowerCase();
  const gene = $('f-gene').value;
  const kind = $('f-kind').value;
  const hap = $('f-hap').value;
  return state.alleles.filter((a) => {
    if (state.verdict && a.verdict !== state.verdict) return false;
    if (gene && a.gene !== gene) return false;
    if (kind === 'candidate' && !isCandidate(a)) return false;
    if (kind === 'trusted' && isCandidate(a)) return false;
    if (hap && a.hap_status !== hap) return false;
    if (query && !`${a.allele_name} ${a.gene} ${a.reason}`.toLowerCase().includes(query)) return false;
    return true;
  });
}

function renderAlleles() {
  const rows = visibleAlleles();
  const list = $('alleles');
  list.replaceChildren();
  $('count').textContent = `Showing ${rows.length} of ${state.alleles.length} alleles.`;
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

function chips(filters) {
  if (!given(filters)) return [];
  return filters.split(',').map((part) => {
    const [code, count] = part.split(':');
    const label = FILTERS[code] || code;
    return el('span', { class: 'chip', text: count && Number(count) > 1 ? `${label} ×${count}` : label });
  });
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
      routeTag(allele),
      isCandidate(allele) ? el('span', { class: 'tag', text: 'candidate' }) : null,
      el('span', { class: 'chips' }, allele.verdict === 'supported' ? [] : chips(allele.filters), hapPill(allele))),
    body);
  // Details are built when first opened: a full V set has hundreds of cards.
  card.addEventListener('toggle', () => {
    if (card.open && !body.childElementCount) body.append(...alleleDetail(allele));
  });
  return card;
}

function alleleDetail(allele) {
  const parts = [el('p', { class: 'reason', text: allele.reason })];
  if (given(allele.notes)) parts.push(el('p', { class: 'muted', text: `Keep in mind: ${allele.notes}.` }));
  const variants = state.variants.get(allele.allele_name) || [];
  if (variants.length) parts.push(variantTable(variants));
  const hap = haplotypePanel(allele);
  if (hap) parts.push(hap);
  const block = state.alignments.get(allele.gene);
  if (block) parts.push(alignmentPanel(block, allele.allele_name));
  return parts;
}

// ---- explaining the score --------------------------------------------------

/** Shows `content` next to `trigger` on hover, keyboard focus or tap.
 *  The explanation is placed on the page, not inside the table, so that a
 *  scrolling table cannot cut it off. */
function explain(trigger, content) {
  explain.count = (explain.count || 0) + 1;
  const tip = el('span', { class: 'hint-text', role: 'tooltip', id: `hint-${explain.count}` }, content);
  trigger.setAttribute('aria-describedby', tip.id);
  const show = () => {
    tip.classList.add('open');
    const at = trigger.getBoundingClientRect();
    const width = tip.offsetWidth;
    const left = Math.min(Math.max(8, at.left + at.width / 2 - width / 2), window.innerWidth - width - 8);
    const below = at.bottom + 8 + tip.offsetHeight < window.innerHeight;
    tip.style.left = `${left}px`;
    tip.style.top = below ? `${at.bottom + 8}px` : `${Math.max(8, at.top - 8 - tip.offsetHeight)}px`;
  };
  const hide = () => tip.classList.remove('open');
  for (const type of ['mouseenter', 'focus']) trigger.addEventListener(type, show);
  for (const type of ['mouseleave', 'blur']) trigger.addEventListener(type, hide);
  trigger.addEventListener('keydown', (event) => { if (event.key === 'Escape') hide(); });
  return el('span', { class: 'hint' }, trigger, tip);
}

/** A small "?" with an explanation. */
function hint(label, content) {
  return explain(el('button', { type: 'button', class: 'hint-button', 'aria-label': label, text: '?' }), content);
}

/** What a score of this size means, in a few words. */
function scoreWords(score) {
  if (score >= 10) return 'typical of real variants';
  if (score >= 3) return 'speaks for a real variant';
  if (score > 1 / 3) return 'says little either way';
  if (score > 0.1) return 'speaks against';
  return 'typical of arbitrary changes';
}

const times = (score) => (score >= 1
  ? `${score >= 10 ? Math.round(score) : score.toFixed(1)} times more often among real variants`
  : `${Math.round(1 / score)} times more often among arbitrary changes`);

/** The score explained, with the levels of this run's calibration. */
function scoreHint() {
  const rows = state.files['score_calibration.tsv'] ? parseTSV(state.files['score_calibration.tsv']).filter((r) => given(r.score)) : [];
  const levels = { none: 'never seen', '>=10%': 'seen in 10% or more of a population' };
  const brief = (score) => (score >= 10 ? Math.round(score) : score >= 1 ? score.toFixed(1) : Number(score.toPrecision(2)));
  // Levels with the same score are one line.
  const lines = [];
  for (const row of rows.reverse()) {
    const name = levels[row.support] || `seen in ${row.support.replace('<', 'under ').replace('-', ' to ')}`;
    const last = lines[lines.length - 1];
    if (last && last.score === row.score) last.names.push(name);
    else lines.push({ score: row.score, names: [name] });
  }
  return el('span', {},
    el('b', { text: 'How to read the score' }),
    el('span', { text: 'Take all variants seen as often as this one. The score says how much more often such a variant belongs to a real allele than is an arbitrary base change (what an error would give).' }),
    lines.length ? el('span', { class: 'hint-table' }, lines.map((line) => el('span', {},
      el('b', { text: brief(Number(line.score)) }),
      el('span', { text: `${line.names.reverse().join(', or ')}: ${scoreWords(Number(line.score))}` })))) : null,
    el('span', { text: '10 or more counts as supported, 0.1 or less as needing confirmation. Around 1 the frequency cannot tell.' }));
}

// ---- variants: one row each, with a frequency strip -------------------------

const LOG_MIN = -6; // strip runs from 1 in a million to 100%
const STRIP = { width: 260, height: 26, pad: 8 };
const xOf = (af) => STRIP.pad + ((Math.max(Math.log10(af), LOG_MIN) - LOG_MIN) / -LOG_MIN) * (STRIP.width - 2 * STRIP.pad);

function mark(kind, x, y) {
  if (kind === 'kgp') return svgIcon(svg('circle', { cx: x, cy: y, r: 4, class: 'm-kgp' }));
  if (kind === 'gnomad') return svgIcon(svg('path', { d: `M${x} ${y - 5}l5 5-5 5-5-5z`, class: 'm-gnomad' }));
  return svgIcon(svg('rect', { x: x - 3.5, y: y - 3.5, width: 7, height: 7, class: 'm-dbsnp' }));
}
function svgIcon(shape) {
  return svg('svg', { viewBox: '0 0 12 12', width: 12, height: 12, 'aria-hidden': 'true' }, shape);
}

/** Every frequency reported for a variant, as {kind, name, af}. */
function frequencies(v) {
  const out = [];
  for (const pop of POPS) {
    const af = number(v[`kgp_AF_${pop}`]);
    if (af > 0) out.push({ kind: 'kgp', name: `1000 Genomes ${POP_NAMES[pop]}`, af });
  }
  const gnomad = number(v.gnomad_max_group_AF);
  if (gnomad > 0) out.push({ kind: 'gnomad', name: `gnomAD ${v.gnomad_max_group}`, af: gnomad });
  if (given(v.dbsnp_freq)) {
    for (const entry of v.dbsnp_freq.split(';')) {
      const [name, af] = entry.split(':');
      if (Number(af) > 0) out.push({ kind: 'dbsnp', name, af: Number(af) });
    }
  }
  return out;
}

const percent = (af) => (af >= 0.01 ? `${(100 * af).toFixed(1)}%` : af >= 0.0001 ? `${(100 * af).toFixed(3)}%` : `${(1e6 * af).toFixed(1)} per million`);

function strip(v) {
  const marks = frequencies(v);
  const best = marks.reduce((top, m) => (m.af > top.af ? m : top), { af: 0 });
  const label = marks.length
    ? `Seen in ${marks.length} population${marks.length > 1 ? 's' : ''}; highest ${percent(best.af)} in ${best.name}.`
    : 'Not seen in any population.';
  const mid = STRIP.height / 2;
  const picture = svg('svg', { viewBox: `0 0 ${STRIP.width} ${STRIP.height}`, class: 'strip', role: 'img', 'aria-label': label });
  picture.append(svg('line', { x1: STRIP.pad, x2: STRIP.width - STRIP.pad, y1: mid, y2: mid, class: 'axis' }));
  for (let power = LOG_MIN; power <= 0; power += 2) {
    picture.append(svg('line', { x1: xOf(10 ** power), x2: xOf(10 ** power), y1: mid - 3, y2: mid + 3, class: 'axis' }));
  }
  picture.append(svg('line', { x1: xOf(0.01), x2: xOf(0.01), y1: 2, y2: STRIP.height - 2, class: 'one-percent' }));
  for (const m of marks) {
    const x = xOf(m.af);
    const shape = m.kind === 'kgp' ? svg('circle', { cx: x, cy: mid, r: 4, class: 'm-kgp' })
      : m.kind === 'gnomad' ? svg('path', { d: `M${x} ${mid - 5}l5 5-5 5-5-5z`, class: 'm-gnomad' })
        : svg('rect', { x: x - 3.5, y: mid - 3.5, width: 7, height: 7, class: 'm-dbsnp' });
    shape.append(svg('title', {}, `${m.name}: ${percent(m.af)}`));
    picture.append(shape);
  }
  return el('div', { class: 'strip-cell' }, picture, el('small', { text: label }));
}

function variantTable(variants) {
  const rsLink = (v) => (given(v.rsid)
    ? el('a', { href: `https://www.ncbi.nlm.nih.gov/snp/${v.rsid}`, target: '_blank', rel: 'noopener', text: v.rsid })
    : el('span', { class: 'muted', text: 'no rsID' }));
  const kind = { SNV: 'base change to', INS: 'insertion,', DEL: 'deletion, leaving' };
  const change = (v) => `${kind[v.type] || v.type} ${v.query_variant}`;
  return el('div', { class: 'scroll' }, el('table', { class: 'variants' },
    el('caption', { text: 'Differences from the reference genome' }),
    el('thead', {}, el('tr', {},
      ...['Variant', 'Position in allele', 'Verdict by frequency', 'Why'].map((h) => el('th', { scope: 'col', text: h })),
      el('th', { scope: 'col' }, 'Score ', hint('What the score means', scoreHint())),
      el('th', { scope: 'col', text: 'Frequency, 1 per million to 100% (log scale)' }))),
    el('tbody', {}, variants.map((v) => el('tr', {},
      el('th', { scope: 'row' }, rsLink(v), el('small', { text: ` ${change(v)}, ${v.contig}:${v.genome_pos}` })),
      el('td', { text: v.aln_pos }),
      el('td', {}, badge(v.verdict)),
      el('td', {}, chips(v.filters)),
      scoreCell(v),
      el('td', {}, strip(v)))))));
}

function scoreCell(v) {
  if (!given(v.score)) return el('td', { class: 'num', text: '–', title: 'Insertions and deletions have no score.' });
  const score = Number(v.score);
  const text = score < 1 ? score.toFixed(3) : score.toFixed(1);
  return el('td', { class: 'num' }, explain(el('span', { class: 'hint-value', tabindex: '0', text }), el('span', {},
    el('b', { text: `${text}: ${scoreWords(score)}` }),
    el('span', { text: `Variants seen as often as this one occur ${times(score)}.` }))));
}

// ---- haplotype check -----------------------------------------------------

function haplotypePanel(allele) {
  const info = HAPLOTYPES[allele.hap_status];
  if (!info || allele.hap_status === 'not_checked') return null;
  const panel = el('div', { class: 'hap' }, el('h4', { text: 'Haplotype check' }), el('p', { text: info.help }));
  const which = given(allele.hap_exact) ? 'exact' : given(allele.hap_carrying) ? 'carrying' : null;
  if (given(allele.hap_untested_str)) {
    panel.append(el('p', { class: 'muted', text: `Tested ${allele.hap_tested} variants. Not in the phased panel: ${allele.hap_untested_str.replace(/;/g, ', ')}.` }));
  }
  if (!which) return panel;
  const total = number(allele[`hap_${which}`]);
  const donors = given(allele.hap_exact_donors) ? `, in ${allele.hap_exact_donors} ${allele.hap_exact_donors === '1' ? 'person' : 'people'}` : '';
  const caption = which === 'exact'
    ? `${total} of ${allele.hap_total} founder chromosomes have exactly this sequence${donors}.`
    : `${total} of ${allele.hap_total} founder chromosomes carry all tested variants.`;
  const shares = POPS.map((pop) => {
    const [count, size] = (allele[`hap_${which}_${pop}`] || '0/0').split('/').map(Number);
    return { pop, count, size, share: size ? count / size : 0 };
  });
  const top = Math.max(...shares.map((s) => s.share), 1e-9);
  panel.append(el('p', { text: caption }),
    el('table', { class: 'bars' }, el('caption', { class: 'sr', text: 'Carrier chromosomes per superpopulation' }),
      el('tbody', {}, shares.map((s) => el('tr', {},
        el('th', { scope: 'row', text: POP_NAMES[s.pop] }),
        el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * s.share) / top}%` })),
        el('td', { class: 'num', text: `${s.count} / ${s.size} (${(100 * s.share).toFixed(s.share < 0.01 ? 2 : 1)}%)` }))))));
  const coincidence = number(allele.hap_coincidence);
  if (which === 'exact' && total > 0 && coincidence >= 0) {
    const share = coincidence < 0.0005 ? 'fewer than 1 in 2000' : `${(100 * coincidence).toFixed(coincidence < 0.01 ? 2 : 1)}% of`;
    panel.append(el('p', { class: 'muted', text: `For comparison: ${share} arbitrary single-base changes of a real chromosome are found on this many chromosomes.` }));
  }
  const p = number(allele.hap_exact_pop_heterogeneity_p);
  if (which === 'exact' && p < 0.001) panel.append(el('p', { class: 'muted', text: 'The frequency differs clearly between populations: this haplotype is population-specific.' }));
  return panel;
}

// ---- alignment -------------------------------------------------------------

function alignmentPanel(block, alleleName) {
  const pre = el('pre', { class: 'alignment', tabindex: '0', 'aria-label': 'Alignment of all alleles of this gene' });
  for (const line of block.split('\n')) {
    const match = line.match(/^(\S+)(\s+)(.*)$/);
    const isRow = match && !line.startsWith('>') && !line.startsWith(' ') && /^[A-Za-z.\- ]+$/.test(match[3]);
    if (!isRow) {
      pre.append(el('span', { class: line.startsWith('>') ? 'aln-title' : 'aln-note', text: `${line}\n` }));
      continue;
    }
    const row = el('span', { class: match[1] === alleleName ? 'aln-row current' : 'aln-row' }, match[1] + match[2]);
    if (match[1] === 'reference') row.append(match[3]);
    else {
      for (const piece of match[3].split(/([A-Za-z-]+)/)) {
        if (piece) row.append(/^[A-Za-z-]+$/.test(piece) ? el('b', { text: piece }) : piece);
      }
    }
    row.append('\n');
    pre.append(row);
  }
  return el('details', { class: 'aln' }, el('summary', { text: 'Alignment of this gene' }),
    el('p', { class: 'muted', text: 'A dot is the same base as the reference; a letter is a difference; a dash is a missing base. This allele is highlighted.' }), pre);
}

// ---- score and checks ------------------------------------------------------

function renderMethod() {
  const view = $('view-method');
  view.replaceChildren();
  const text = state.files['score_calibration.tsv'];
  if (!text) {
    view.append(el('p', { class: 'muted', text: 'No score table in these results.' }));
    return;
  }
  const rows = parseTSV(text);
  const source = { run: 'calibrated on this run', 'built-in': 'the built-in human IGHV table (this run has too few trusted variants of its own)', file: 'a table from an earlier run' }[rows[0].score_source] || rows[0].score_source;
  view.append(
    el('h3', { text: 'How the score was made' }),
    el('p', { text: 'The score of a difference is a likelihood ratio: how many times more often its level of population support is found among the variants of real alleles than among arbitrary single-base changes in the same genes. Above 1 the data speak for a real variant, below 1 against; near 1 they say nothing either way. Real alleles are the named trusted ones and those confirmed as haplotypes (two or more differences found together on phased chromosomes).' }),
    el('p', {}, 'Scores in use: ', el('b', { text: source }), '.'),
    el('div', { class: 'scroll' }, el('table', { class: 'calibration' },
      el('caption', { text: 'Share of arbitrary changes and of trusted variants at each level of support (counts of this run)' }),
      el('thead', {}, el('tr', {}, ...['Support: frequency guaranteed in the best population', 'Arbitrary changes', 'Variants of real alleles', 'Score'].map((h) => el('th', { scope: 'col', text: h })))),
      el('tbody', {}, rows.map((r) => el('tr', {},
        el('th', { scope: 'row', text: r.support }),
        shareCell(r.arbitrary_fraction, r.arbitrary_changes, 'arbitrary'),
        shareCell(r.trusted_fraction, r.trusted_variants, 'trusted'),
        el('td', { class: 'num', text: r.score })))))),
    el('p', { class: 'muted' }, 'By frequency, a difference is ', badge('supported'), ' from a score of 10 or when it is shown to reach 1% in a population, ', badge('needs_confirmation'), ' up to 0.1, and ', badge('inconclusive'), ' in between.'));
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

// ---- report (a small Markdown subset: the one igpop writes) -----------------

function inline(text) {
  const out = document.createDocumentFragment();
  for (const piece of text.split(/(`[^`]+`|\*\*[^*]+\*\*)/)) {
    if (piece.startsWith('`')) out.append(el('code', { text: piece.slice(1, -1) }));
    else if (piece.startsWith('**')) out.append(el('b', { text: piece.slice(2, -2) }));
    else if (piece) out.append(piece);
  }
  return out;
}

function renderReport() {
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

// ---- files -----------------------------------------------------------------

function renderFiles() {
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

// ---- tabs and events -------------------------------------------------------

function selectTab(id) {
  for (const tab of document.querySelectorAll('[role=tab]')) {
    const selected = tab.id === id;
    tab.setAttribute('aria-selected', String(selected));
    tab.tabIndex = selected ? 0 : -1;
    $(tab.getAttribute('aria-controls')).hidden = !selected;
  }
}

function wire() {
  const drop = $('drop');
  const input = $('file');
  drop.addEventListener('click', () => input.click());
  drop.addEventListener('keydown', (event) => {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      input.click();
    }
  });
  input.addEventListener('change', () => handleFiles(input.files).then(() => { input.value = ''; }));
  for (const type of ['dragenter', 'dragover']) {
    document.addEventListener(type, (event) => {
      event.preventDefault();
      drop.classList.add('over');
    });
  }
  for (const type of ['dragleave', 'drop']) document.addEventListener(type, () => drop.classList.remove('over'));
  document.addEventListener('drop', (event) => {
    event.preventDefault();
    if (!$('progress').hidden) return;
    handleFiles(event.dataTransfer.files);
  });
  for (const button of document.querySelectorAll('[data-restart]')) button.addEventListener('click', () => show('start'));
  for (const id of ['f-search', 'f-gene', 'f-kind', 'f-hap']) $(id).addEventListener('input', renderAlleles);
  const tabs = [...document.querySelectorAll('[role=tab]')];
  tabs.forEach((tab, i) => {
    tab.addEventListener('click', () => selectTab(tab.id));
    tab.addEventListener('keydown', (event) => {
      const step = { ArrowRight: 1, ArrowLeft: -1 }[event.key];
      if (!step) return;
      const next = tabs[(i + step + tabs.length) % tabs.length];
      selectTab(next.id);
      next.focus();
    });
  });
}

wire();
// The worker is started with the page, so the analysis module is loaded
// and the stored data counted before the first file is dropped.
startWorker();

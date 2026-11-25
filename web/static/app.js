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
  sequence_not_in_panel: 'each difference is known, but no phased chromosome has this exact sequence',
  no_coordinates: 'no coordinates for its gene',
  low_identity: 'too different from its gene',
};
// How well a sequence could be placed on the gene its name gives.
const PLACEMENTS = {
  ok: { label: 'Placed on its gene', help: 'The sequence fits the gene of its name better than any other gene of the list.' },
  ambiguous_locus: { label: 'Fits another gene too', help: 'The sequence fits another gene of the list at least as well as the gene of its name. The genes are near-identical copies, so its differences, and their frequencies, may belong to the other gene.' },
  low_identity: { label: 'Too different from its gene', help: 'The sequence is too different from the gene of its name to compare; no differences are called.' },
  no_coordinates: { label: 'Gene not in the list', help: 'The gene of its name is not in the gene list, so the sequence could not be placed.' },
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
const OUTPUTS = ['report.md', 'final.tsv', 'alleles.tsv', 'variants.tsv', 'alignments.tsv', 'alignments.txt', 'score_model.tsv', 'score_levels.tsv', 'haplotype_calibration.tsv', 'dbsnp_datasets.tsv'];
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
  'score_model.tsv': 'The curve behind the score',
  'score_levels.tsv': 'Real variants and possible changes counted per level of support',
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

const state = { limit: 10, prefix: 'igpop', marker: '_S', build: '', files: {}, alleles: [], variants: new Map(), alignments: new Map(), verdict: '' };
let worker = null;

function show(section) {
  for (const id of ['start', 'check', 'progress', 'error', 'results']) $(id).hidden = id !== section;
  if (section !== 'progress') clearInterval(run.timer);
}

function fail(message) {
  $('error-text').textContent = message;
  show('error');
  $('error').querySelector('button').focus();
}

// ---- starting a run ------------------------------------------------------

function startWorker() {
  // Every file of the app carries the same build stamp in its address, so
  // that a browser never combines a new page with an older, cached worker
  // or analysis module.
  const stamp = new URL(import.meta.url).search;
  worker = new Worker(new URL(`./worker.js${stamp}`, import.meta.url), { type: 'module' });
  worker.onmessage = (event) => {
    const message = event.data;
    if (message.type === 'ready') {
      if (typeof message.geneList !== 'string') return fail('This page and its analysis module are from different versions. Reload the page (Ctrl+Shift+R, or Cmd+Shift+R on a Mac) and try again.');
      $('version').textContent = `igpop ${message.version}`;
      names.builtin = { coords: message.geneList, aliases: message.aliases };
      renderGeneList();
    }
    if (message.type === 'checked' && message.id === names.asked) renderCheck(message.genes);
    if (message.type === 'check-failed' && message.id === names.asked) checkFailed(message.message);
    if ('stored' in message) showStored(message.stored, message.reading);
    if (message.type === 'log') logLine(message.message);
    else if (message.type === 'error') fail(message.message);
    else if (message.type === 'done') display(message.output.files, message.output.warnings, message.output.summary);
  };
  worker.onerror = (event) => fail(`The analysis module could not be started (${event.message || 'unknown error'}). A recent Chrome, Edge, Firefox or Safari is needed.`);
}

/** How much downloaded data this browser keeps, with a way to delete it. */
function showStored(bytes, reading) {
  const box = $('stored');
  // The same in a few words, for the progress view.
  run.store = bytes == null ? 'Downloaded data cannot be kept in this browser, so this run downloads everything.'
    : reading ? 'igpop is open in another tab: data kept in this browser are used here, and what this tab downloads is kept until it is closed.'
      : '';
  if (bytes == null) {
    box.textContent = 'This browser does not keep downloaded data between visits, so each run downloads again.';
    return;
  }
  const size = bytes < 1e6 ? `${Math.round(bytes / 1e3)} kB` : `${(bytes / 1e6).toFixed(0)} MB`;
  run.kept = bytes ? `${size} of population data are kept in this browser from earlier runs` : 'No population data are kept in this browser yet';
  box.replaceChildren(
    bytes ? `${size} of population data kept in this browser. ` : '',
    reading ? 'igpop is also open in another tab: what this tab downloads is kept only until it is closed. ' : '',
    bytes && !reading ? el('button', { type: 'button', class: 'link', text: 'Delete', onclick: () => worker.postMessage({ type: 'clear' }) }) : '',
  );
}

// The stages of a run, in the order they happen, and the progress lines
// of the analysis module that announce each of them.
const STEPS = [
  { id: 'place', label: 'Placing alleles on their genes', now: 'Reading the reference sequence of each gene and aligning your alleles to it.', starts: /reference:|^aligning/ },
  { id: 'download', label: 'Fetching population data', now: 'Fetching the regions of your genes from 1000 Genomes, gnomAD and dbSNP.', starts: /^downloading|^progress:/ },
  { id: 'read', label: 'Reading frequencies and haplotypes', now: 'Looking up every difference in the population data.', starts: /^(1kGP|gnomAD|dbSNP|phased panel)/ },
  { id: 'judge', label: 'Judging the alleles', now: 'Scoring every difference and writing the results.', starts: /^score:/ },
];
const run = { step: -1, started: 0, timer: null, sources: new Map(), store: '', kept: '' };

function startProgress() {
  run.step = -1;
  run.sources = new Map();
  $('sources').hidden = true;
  $('progress-store').textContent = run.store;
  run.started = Date.now();
  $('steps').replaceChildren(...STEPS.map((step, i) => el('li', { id: `step-${step.id}` },
    el('span', { class: 'step-mark', 'aria-hidden': 'true', text: String(i + 1) }), el('span', { text: step.label }))));
  $('log').replaceChildren(el('li', { class: 'done', text: run.store || run.kept || 'Downloaded data are kept in this browser' }));
  setBar(null);
  setStep(0);
  clearInterval(run.timer);
  const tick = () => {
    const seconds = Math.round((Date.now() - run.started) / 1000);
    $('progress-time').textContent = seconds < 60 ? `${seconds} s` : `${Math.floor(seconds / 60)} min ${seconds % 60} s`;
  };
  tick();
  run.timer = setInterval(tick, 1000);
}

/** Steps only move forward: a later line of an earlier kind changes nothing. */
function setStep(index) {
  if (index <= run.step) return;
  run.step = index;
  STEPS.forEach((step, i) => {
    const item = $(`step-${step.id}`);
    item.className = i < index ? 'done' : i === index ? 'current' : '';
    if (i === index) item.setAttribute('aria-current', 'step');
    else item.removeAttribute('aria-current');
  });
  $('progress-now').textContent = STEPS[index].now;
  setBar(null);
}

/** `share` between 0 and 1, or null while the length of a step is unknown. */
function setBar(share) {
  const bar = $('bar');
  bar.classList.toggle('known', share != null);
  if (share == null) bar.removeAttribute('aria-valuenow');
  else {
    bar.setAttribute('aria-valuenow', String(Math.round(100 * share)));
    bar.firstElementChild.style.width = `${Math.max(2, 100 * share)}%`;
  }
  if (share == null) bar.firstElementChild.style.width = '';
}

/** The source a data file belongs to, in the words of the results. */
function sourceOf(url) {
  const u = url.toLowerCase();
  if (u.includes('gnomad')) return /2\.1/.test(u) ? 'gnomAD v2.1 (GRCh37)' : 'gnomAD v4.1';
  if (/\/snp\/|gcf_/.test(u)) return u.includes('gcf_000001405.25') ? 'dbSNP (GRCh37)' : 'dbSNP';
  if (u.includes('phase3')) return '1000 Genomes phase 3 (GRCh37)';
  if (u.includes('phased')) return 'Phased chromosomes of 1000 Genomes';
  if (/1000g|1000genomes/.test(u)) return '1000 Genomes, 30x';
  return url.split('/').pop();
}

/** The table of sources under the bar: how far each one is, and how much
 *  of it was already in the browser. */
function renderSources() {
  const table = $('sources');
  table.hidden = !run.sources.size || run.step !== 1;
  table.tBodies[0].replaceChildren(...[...run.sources].map(([name, s]) => {
    const done = s.kept + s.downloaded + s.failed;
    const state = !done ? 'waiting'
      : [s.downloaded ? `${s.downloaded} downloaded` : '', s.kept ? `${s.kept} already in this browser` : '', s.failed ? `${s.failed} to be read later` : ''].filter(Boolean).join(', ');
    return el('tr', { class: done >= s.total ? 'complete' : done ? 'active' : '' },
      el('th', { scope: 'row', text: name }),
      el('td', {}, el('span', { class: 'source-track', role: 'img', 'aria-label': `${done} of ${s.total} regions` },
        el('i', { class: 'kept', style: `width:${(100 * s.kept) / s.total}%` }),
        el('i', { class: 'new', style: `width:${(100 * (s.downloaded + s.failed)) / s.total}%` }))),
      el('td', { class: 'num', text: `${done} of ${s.total}` }),
      el('td', { class: 'muted', text: state }));
  }));
}

// The lines of the analysis module in the words of the page. A line that
// matches nothing is shown as it is, without addresses.
const STEP_WORDS = [
  [/^(\S+) reference: .*/, (m) => `Reference genome ${m[1]}`],
  [/^aligning (\d+) alleles to (\d+) loci on (\S+)/, (m) => `${m[1]} alleles aligned to ${m[2]} genes on ${m[3]}`],
  [/^finding the gene of (\d+) alleles by sequence/, (m) => `Gene of ${m[1]} alleles found by sequence`],
  [/^downloading (\d+) regions from (\d+) servers.*/, (m) => `Population data for ${m[1]} regions of the genome`],
  [/^phased panel: .*/, () => 'Phased chromosomes read'],
  [/^(1kGP\w*|gnomAD\S*|dbSNP): .*/, (m) => `${sourceTitle(m[1])} read`],
  [/^(\d+) founders, (\d+) haplotypes/, (m) => `${Number(m[2]).toLocaleString('en')} chromosomes of ${Number(m[1]).toLocaleString('en')} unrelated people`],
  [/^(\d+) regions, (\d+) from cache/, (m) => `${m[1]} regions`],
  [/^score: calibrated from (.*)/, (m) => `Score calibrated from ${m[1].replace(/^this run: /, '')}`],
];

function logLine(message) {
  const planned = message.match(/^plan: (\d+) regions of (\S+)/);
  if (planned) {
    const name = sourceOf(planned[2]);
    const source = run.sources.get(name) || { total: 0, kept: 0, downloaded: 0, failed: 0 };
    source.total += Number(planned[1]);
    run.sources.set(name, source);
    return;
  }
  const text = message.replace(/https?:\/\/\S+\/([^/\s]+)/g, '$1');
  const step = STEPS.findIndex((s) => s.starts.test(text.trim()));
  if (step >= 0) setStep(step);
  const counted = message.match(/^progress: (\d+) of (\d+)(?: regions ready, (\d+) kept[^(]*(?:\((kept|downloaded|failed) (\S+)\))?)?/);
  if (counted) {
    // One line per region: shown as bars, not as text.
    const [done, total, kept] = [Number(counted[1]), Number(counted[2]), Number(counted[3] || 0)];
    const source = counted[5] && run.sources.get(sourceOf(counted[5]));
    if (source) source[counted[4]] += 1;
    if (run.step === 1) {
      setBar(done / total);
      $('progress-now').textContent = `${done} of ${total} regions ready: ${done - kept} downloaded, ${kept} already in this browser.`;
      renderSources();
    }
    return;
  }
  if (run.step !== 1) $('sources').hidden = true;
  const dataset = text.match(/^(1kGP\w*|gnomAD\S*|dbSNP|phased panel):/);
  if (dataset && run.step === 2) $('progress-now').textContent = `Reading ${dataset[1] === 'phased panel' ? 'the phased chromosomes' : sourceTitle(dataset[1])}.`;
  const aligning = text.match(/^aligning (\d+) alleles to (\d+) loci on (\S+)/);
  if (aligning && run.step === 0) $('progress-now').textContent = `Aligning ${aligning[1]} alleles to ${aligning[2]} genes on ${aligning[3]}.`;
  const worded = STEP_WORDS.map(([pattern, words]) => [text.trim().match(pattern), words]).find(([match]) => match);
  const line = worded ? worded[1](worded[0]) : text.trim();
  const log = $('log');
  const last = log.lastElementChild;
  if (text.startsWith('  ') && last) {
    last.append(el('small', { text: `: ${line}` }));
  } else {
    if (last) last.classList.add('done');
    log.append(el('li', { text: line }));
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
  names.pending = { fasta, name: file.name };
  $('check-name').textContent = file.name;
  $('check-summary').textContent = 'Reading the gene names…';
  for (const id of ['check-evidence', 'check-map']) $(id).replaceChildren();
  names.renames = {};
  $('check-run').disabled = true;
  show('check');
  askCheck();
}

// ---- gene names: checked before a run -----------------------------------------

// The gene list and the names of this browser. `rows` is the user's own
// gene list (null = the built-in one); `mapping` says what to do with a
// gene name the list does not have: another gene's name, 'sequence' (find
// the gene by sequence) or 'skip' (leave its alleles out).
const names = { builtin: { coords: '', aliases: '' }, rows: null, mapping: {}, renames: {}, open: new Set(), pending: null, alleles: [], genes: [], asked: 0 };
const KEPT = 'igpop-gene-names';
try {
  const kept = JSON.parse(localStorage.getItem(KEPT) || '{}');
  if (Array.isArray(kept.rows)) names.rows = kept.rows;
  if (kept.mapping && typeof kept.mapping === 'object') names.mapping = kept.mapping;
} catch (error) { /* nothing kept, or storage is off */ }
function keepNames() {
  try { localStorage.setItem(KEPT, JSON.stringify({ rows: names.rows, mapping: names.mapping })); } catch (error) { /* storage is off */ }
}

const COLUMNS = ['gene', 'build', 'contig', 'start', 'end', 'strand'];
/** Rows of a gene list in text form; throws with the line that is wrong. */
function parseGeneList(text) {
  const rows = [];
  text.split(/\r?\n/).forEach((line, i) => {
    if (!line.trim() || line.startsWith('#') || /^gene\b/i.test(line)) return;
    const cells = line.trim().split(/\t|\s{2,}|,/).map((c) => c.trim());
    const [gene, build, contig, start, end, strand = '.'] = cells;
    if (cells.length < 5 || !Number.isInteger(Number(start)) || !Number.isInteger(Number(end)) || Number(start) < 1 || Number(end) < Number(start)) {
      throw new Error(`Line ${i + 1} is not "gene, build, contig, start, end, strand": ${line.trim().slice(0, 80)}`);
    }
    rows.push({ gene, build, contig, start: Number(start), end: Number(end), strand: ['+', '-'].includes(strand) ? strand : '.' });
  });
  if (!rows.length) throw new Error('The list has no genes.');
  return rows;
}
const geneRows = () => names.rows || (names.builtin.coords ? parseGeneList(names.builtin.coords) : []);
const geneListText = (rows) => [COLUMNS.join('\t'), ...rows.map((r) => COLUMNS.map((c) => r[c]).join('\t'))].join('\n');
// What the user can decide for a name: a gene of the list, BY_SEQUENCE
// (find the gene by sequence during the run) or LEFT_OUT.
const BY_SEQUENCE = '?';
const LEFT_OUT = 'skip';
for (const [key, value] of Object.entries(names.mapping)) if (value === 'sequence') names.mapping[key] = BY_SEQUENCE;

/** The decision that applies to an allele: its own, else that of its name. */
const decision = (allele) => names.mapping[allele.name] ?? names.mapping[allele.named];
/** The alias table of a run: the built-in one and the user's decisions. */
const aliasText = () => [names.builtin.aliases, ...Object.entries(names.mapping)
  .filter(([, to]) => to !== LEFT_OUT).map(([from, to]) => `${from}\t${to}`)].join('\n');
/** The FASTA text with the names the user corrected. */
function currentFasta() {
  const { fasta } = names.pending;
  if (!Object.keys(names.renames).length) return fasta;
  return fasta.replace(/^>(\S+)/gm, (line, name) => (names.renames[name] ? `>${names.renames[name]}` : line));
}

function askCheck() {
  if (!names.pending || !worker) return;
  names.asked += 1;
  worker.postMessage({
    type: 'check', id: names.asked, fasta: currentFasta(), name: names.pending.name, build: $('opt-build').value,
    haplotypes: $('opt-hap').checked,
    gnomad: $('opt-gnomad').checked,
    coords: names.rows ? geneListText(names.rows) : undefined, aliases: aliasText(),
  });
}

function checkFailed(message) {
  $('check-summary').textContent = 'The gene list could not be used.';
  $('check-error').textContent = message;
  $('check-error').hidden = false;
  $('check-run').disabled = true;
}

const sortGenes = (a, b) => a.localeCompare(b, 'en', { numeric: true });
const plural = (n, word) => `${n} ${word}${n === 1 ? '' : 's'}`;

/** Where the evidence of the run will come from: one card per group of
 *  genes that share their data, with the sources of each. */
function renderEvidencePlan(groups, genes, sources) {
  const used = (gene) => groups.filter((g) => g.gene === gene.gene && g.state !== 'out');
  const count = (list) => `${plural(list.length, 'gene')}, ${plural(list.reduce((n, gene) => n + used(gene).reduce((k, g) => k + g.alleles.length, 0), 0), 'allele')}`;
  const builds = [...new Set(sources.map((s) => s.build))];
  const cards = builds.map((build, i) => {
    const members = genes.filter((g) => g.build === build && g.evidence && used(g).length);
    const mine = sources.filter((s) => s.build === build);
    if (!members.length) return null;
    return el('div', { class: 'plan-card' },
      el('h3', {}, `Genes on ${build}`, el('small', { text: count(members) })),
      i > 0 ? el('p', { class: 'muted', text: `These genes are missing from ${builds[0]}. They are judged with the older data of ${build}, so read their verdicts with more caution.` }) : null,
      el('ul', { class: 'plan-sources' }, mine.map((s) => el('li', {},
        el('b', { text: s.title }),
        el('span', { class: 'plan-people', text: s.people ? `${s.people.toLocaleString('en')} people` : '' }),
        el('small', { class: 'muted' }, el('span', { class: `plan-kind k-${s.gives}`, text: s.gives === 'haplotypes' ? 'whole sequences' : 'frequencies' }), ` ${s.about}`)))),
      mine.some((s) => s.gives === 'haplotypes') ? null : el('p', { class: 'muted', text: 'The haplotype check is switched off.' }));
  }).filter(Boolean);
  const apart = genes.filter((g) => g.build && !g.evidence && used(g).length);
  if (apart.length) {
    cards.push(el('div', { class: 'plan-card plan-none' },
      el('h3', {}, el('span', {}, el('span', { class: 'plan-warn', 'aria-hidden': 'true', text: '!' }), ' No population data for these genes'), el('small', { text: count(apart) })),
      el('p', { text: `${apart.map((g) => g.gene).join(', ')}: the reference has ${apart.length === 1 ? 'this gene' : 'these genes'} only on a piece of sequence that is not part of a chromosome. The alleles are compared with it, but there are no dependable population data there: an allele that differs from it will be "not assessed".` })));
  }
  $('check-evidence').replaceChildren(cards.length ? el('section', { class: 'plan', 'aria-label': 'Where the evidence will come from' },
    el('h3', { class: 'plan-title', text: 'Where the evidence will come from' }),
    el('div', { class: 'plan-cards' }, cards)) : '');
}

// How a name of the file is connected to a gene of the list.
const LINKS = {
  same: { label: 'same name', arrow: 'solid' },
  rule: { label: 'built-in rule', arrow: 'solid', help: 'The gene list knows these alleles under another gene name.' },
  yours: { label: 'your match', arrow: 'solid' },
  sequence: { label: 'by sequence', arrow: 'dashed', help: 'Each allele will be placed on the gene its sequence fits best. The results say which.' },
  out: { label: 'left out', arrow: 'broken', help: 'These alleles are not analysed.' },
  none: { label: 'not matched', arrow: 'broken', help: 'No gene of the list has this name. Choose a gene, let igpop find it by sequence, correct the name, or add the gene to the list.' },
};

/** How an allele is connected to the gene list: a key of LINKS. */
function linkOf(allele, place) {
  const chosen = decision(allele);
  return chosen === LEFT_OUT ? 'out'
    : allele.gene === BY_SEQUENCE ? 'sequence'
      : !(place && place.build) ? 'none'
        : allele.gene === allele.named ? 'same'
          : chosen === allele.gene ? 'yours' : 'rule';
}

/** The alleles of the file in groups that share name, gene and state. */
function mappingGroups(alleles, genes) {
  const known = new Map(genes.map((g) => [g.gene, g]));
  const groups = new Map();
  for (const allele of alleles) {
    const place = known.get(allele.gene);
    const state = linkOf(allele, place);
    const key = `${allele.named}\n${allele.gene}\n${state}`;
    if (!groups.has(key)) groups.set(key, { key, named: allele.named, gene: allele.gene, state, place, alleles: [] });
    groups.get(key).alleles.push(allele);
  }
  const order = ['none', 'out', 'sequence', 'yours', 'rule', 'same'];
  return [...groups.values()].sort((a, b) => order.indexOf(a.state) - order.indexOf(b.state) || sortGenes(a.named, b.named) || sortGenes(a.gene, b.gene));
}

/** Apply a decision to the alleles of a group: to the whole name if the
 *  group holds every allele of that name, else to each allele. */
function decide(group, value) {
  const whole = names.alleles.filter((a) => a.named === group.named).length === group.alleles.length;
  for (const allele of group.alleles) delete names.mapping[allele.name];
  if (whole) {
    if (value) names.mapping[group.named] = value;
    else delete names.mapping[group.named];
  } else for (const allele of group.alleles) {
    // Without a decision of its own an allele follows its name again.
    if (value) names.mapping[allele.name] = value;
  }
  keepNames();
  askCheck();
}

/** The choices for a group or an allele, as a select. The list of all
 *  genes is long, and there can be hundreds of these selects: the genes
 *  are added when the select is first used. */
function decisionSelect(label, current, similarTo, onChange) {
  const fixed = [['', 'As the name says'], [BY_SEQUENCE, 'Find the gene by sequence'], [LEFT_OUT, 'Leave out of the analysis']];
  const select = el('select', { 'aria-label': label }, fixed.map(([value, text]) => el('option', { value, text })));
  if (current && !fixed.some(([value]) => value === current)) select.append(el('option', { value: current, text: `Alleles of ${current}` }));
  select.value = current || '';
  const fill = () => {
    if (select.dataset.filled) return;
    select.dataset.filled = 'yes';
    const all = [...new Set(geneRows().map((r) => r.gene))].sort(sortGenes);
    const near = all.filter((g) => g.startsWith(similarTo) || similarTo.startsWith(g) || g.replace(/[A-Z]$/, '') === similarTo);
    const chosen = select.value;
    for (const option of [...select.options].slice(fixed.length)) option.remove();
    if (near.length) select.append(el('optgroup', { label: 'Genes with a similar name' }, near.map((n) => el('option', { value: n, text: `Alleles of ${n}` }))));
    select.append(el('optgroup', { label: 'All genes of the list' }, all.filter((n) => !near.includes(n)).map((n) => el('option', { value: n, text: `Alleles of ${n}` }))));
    select.value = chosen;
  };
  for (const event of ['pointerdown', 'focus', 'keydown']) select.addEventListener(event, fill, { once: true });
  select.addEventListener('change', () => onChange(select.value));
  return select;
}

const short = (allele) => allele.name.slice(allele.named.length) || allele.name;
const original = (name) => Object.keys(names.renames).find((from) => names.renames[from] === name);

function mappingRow(group, sources) {
  const link = LINKS[group.state];
  const whole = names.alleles.filter((a) => a.named === group.named).length === group.alleles.length;
  const chosen = whole ? names.mapping[group.named] : (group.alleles.every((a) => names.mapping[a.name] === names.mapping[group.alleles[0].name]) ? names.mapping[group.alleles[0].name] : undefined);
  const target = group.state === 'sequence' ? el('span', { class: 'muted', text: 'found during the run' })
    : group.state === 'out' ? el('span', { class: 'muted', text: 'not analysed' })
      : group.state === 'none' ? el('span', { class: 'map-missing', text: 'no gene of the list' })
        : el('span', {}, el('b', { text: group.gene }), el('small', { class: 'muted', text: !group.place.evidence ? ' not on a chromosome: no population data'
          : ` ${group.place.build}${sources.some((x) => x.build === group.place.build) && group.place.build !== (sources[0] || {}).build ? ', older data' : ''}` }));
  const toggle = el('button', { type: 'button', class: 'link', text: 'correct names', onclick: () => showNameEditor(group.named) });
  return [el('tr', { class: `map-row m-${group.state}` },
    el('th', { scope: 'row' }, group.named,
      el('small', { class: 'muted', text: ` ${plural(group.alleles.length, 'allele')}: ${group.alleles.slice(0, 4).map(short).join(', ')}${group.alleles.length > 4 ? ', …' : ''}` }),
      group.state === 'rule' ? el('small', { class: 'map-rule', text: `The gene list has no gene called ${group.named}. It knows these alleles as ${group.gene}, so they are sent there automatically.` }) : null),
    el('td', { class: 'map-link' }, el('span', { class: `arrow a-${link.arrow}`, 'aria-hidden': 'true' }), el('span', { class: 'map-how', title: link.help || '', text: link.label })),
    el('td', {}, target),
    el('td', {}, decisionSelect(`What to do with ${group.named}`, chosen, group.named, (value) => decide(group, value)), ' ', toggle))];
}

// ---- the names of the sequences, corrected one by one or all at once ----------

// One line per sequence, in the order of the file. The lines are built
// once per file and only their content changes afterwards, so that a line
// never moves or disappears while it is being edited.
const editor = { lines: new Map(), built: null };

/** Give the sequence first called `first` the name `next`; false if refused. */
function rename(first, next) {
  const now = names.renames[first] || first;
  if (!next || next === now) return true;
  if (names.alleles.some((a) => a.name === next)) {
    $('check-error').textContent = `Another sequence of the file is already called ${next}.`;
    $('check-error').hidden = false;
    return false;
  }
  // A decision made for the sequence follows its name.
  if (now in names.mapping) {
    names.mapping[next] = names.mapping[now];
    delete names.mapping[now];
    keepNames();
  }
  if (next === first) delete names.renames[first];
  else names.renames[first] = next;
  return true;
}

function buildNameEditor() {
  editor.built = names.pending;
  editor.lines = new Map();
  const body = el('tbody');
  for (const allele of names.alleles) {
    const first = allele.name;
    const input = el('input', { type: 'text', value: first, size: Math.max(20, first.length + 4), 'aria-label': `Name of ${first}`, spellcheck: 'false' });
    input.addEventListener('change', () => {
      if (rename(first, input.value.trim().split(/\s/)[0])) askCheck();
      else input.value = names.renames[first] || first;
    });
    const status = el('td', { class: 'name-status' });
    const was = el('td', { class: 'muted' });
    const choice = el('td');
    const row = el('tr', {}, el('td', {}, input), status, choice, was);
    editor.lines.set(first, { row, input, status, choice, was });
    body.append(row);
  }
  $('name-table').replaceChildren(
    el('thead', {}, el('tr', {}, ...['Name used for the run', 'Leads to', 'Or decide for this sequence', ''].map((h) => el('th', { scope: 'col', text: h })))), body);
}

/** Bring every line of the editor up to date with the last check. */
function updateNameEditor() {
  if (editor.built !== names.pending) buildNameEditor();
  const known = new Map(names.genes.map((g) => [g.gene, g]));
  const show = $('name-show').value;
  const wanted = $('name-filter').value.trim().toLowerCase();
  let changed = 0;
  let shown = 0;
  for (const allele of names.alleles) {
    const first = original(allele.name) || allele.name;
    const line = editor.lines.get(first);
    if (!line) continue;
    const state = linkOf(allele, known.get(allele.gene));
    if (document.activeElement !== line.input) line.input.value = allele.name;
    line.row.className = `n-${state}${first !== allele.name ? ' n-changed' : ''}`;
    const shownAs = `${state}\n${allele.named}\n${allele.gene}`;
    if (line.status.dataset.shown !== shownAs) line.status.replaceChildren(
      el('span', { class: `name-mark k-${state}`, 'aria-hidden': 'true', text: state === 'none' ? '×' : state === 'out' ? '–' : state === 'sequence' ? '?' : '✓' }),
      state === 'none' ? `no gene called ${allele.named}` : state === 'out' ? 'left out' : state === 'sequence' ? 'found by sequence during the run'
        : `${allele.gene}${state === 'rule' ? ' (built-in rule)' : state === 'yours' ? ' (your match)' : ''}`);
    line.status.dataset.shown = shownAs;
    // The select of a line is built once per name and decision, not on every check.
    const made = `${allele.name}\n${names.mapping[allele.name] || ''}`;
    if (line.choice.dataset.made !== made) {
      line.choice.dataset.made = made;
      line.choice.replaceChildren(decisionSelect(`What to do with ${allele.name}`, names.mapping[allele.name], allele.named, (value) => {
        if (value) names.mapping[allele.name] = value;
        else delete names.mapping[allele.name];
        keepNames();
        askCheck();
      }));
    }
    line.was.replaceChildren(...(first !== allele.name ? [`was ${first} `, el('button', { type: 'button', class: 'link', text: 'undo', onclick: () => { rename(first, first); askCheck(); } })] : []));
    if (first !== allele.name) changed += 1;
    const visible = (!show || (show === 'open' && state === 'none') || (show === 'changed' && first !== allele.name))
      && (!wanted || allele.name.toLowerCase().includes(wanted) || first.toLowerCase().includes(wanted));
    // A line being edited stays, even if its new name no longer fits the filter.
    line.row.hidden = !visible && document.activeElement !== line.input;
    if (!line.row.hidden) shown += 1;
  }
  $('name-editor-state').textContent = changed ? `: ${plural(changed, 'name')} changed` : '';
  previewReplace(shown);
}

/** Open the editor on the sequences of one gene name. */
function showNameEditor(named) {
  $('name-show').value = '';
  $('name-filter').value = named;
  $('name-editor').open = true;
  updateNameEditor();
  $('name-editor').scrollIntoView({ block: 'start', behavior: 'smooth' });
}

/** What "find and replace" would do to the names shown: [first, next] pairs,
 *  or an error message for a pattern that cannot be read. */
function replacements() {
  const form = new FormData($('name-replace'));
  const find = String(form.get('find') || '');
  if (!find) return [];
  let pattern;
  try {
    pattern = form.get('regex') ? new RegExp(find, 'g') : null;
  } catch (error) {
    return `That is not a valid regular expression: ${error.message}`;
  }
  const out = [];
  for (const [first, line] of editor.lines) {
    if (line.row.hidden) continue;
    const now = names.renames[first] || first;
    const next = (pattern ? now.replace(pattern, String(form.get('replace') || '')) : now.split(find).join(String(form.get('replace') || ''))).trim().split(/\s/)[0];
    if (next && next !== now) out.push([first, next]);
  }
  return out;
}

function previewReplace(shown) {
  const planned = replacements();
  $('name-replace-preview').textContent = typeof planned === 'string' ? planned
    : new FormData($('name-replace')).get('find') ? `${plural(planned.length, 'name')} of the ${shown ?? [...editor.lines.values()].filter((l) => !l.row.hidden).length} shown would change${planned.length ? `, for example ${names.renames[planned[0][0]] || planned[0][0]} to ${planned[0][1]}` : ''}.` : '';
}

function wireNameEditor() {
  for (const id of ['name-show', 'name-filter']) $(id).addEventListener('input', updateNameEditor);
  $('name-replace').addEventListener('input', () => previewReplace());
  $('name-replace').addEventListener('submit', (event) => {
    event.preventDefault();
    const planned = replacements();
    if (typeof planned === 'string' || !planned.length) return;
    // All or nothing: two sequences must not end up with one name.
    const after = new Map(names.alleles.map((a) => [original(a.name) || a.name, a.name]));
    for (const [first, next] of planned) after.set(first, next);
    const taken = new Set();
    for (const name of after.values()) {
      if (taken.has(name)) {
        $('check-error').textContent = `Nothing was replaced: two sequences would both be called ${name}.`;
        $('check-error').hidden = false;
        return;
      }
      taken.add(name);
    }
    $('check-error').hidden = true;
    for (const [first, next] of planned) {
      const now = names.renames[first] || first;
      if (now in names.mapping) {
        names.mapping[next] = names.mapping[now];
        delete names.mapping[now];
      }
      if (next === first) delete names.renames[first];
      else names.renames[first] = next;
    }
    keepNames();
    askCheck();
  });
}

/** The check of the file: how its names connect to the gene list. Called
 *  with a new answer of the analysis module, or without to draw again. */
function renderCheck(result) {
  if (result) Object.assign(names, { alleles: result.alleles, genes: result.genes, sources: result.sources });
  const { alleles, genes, sources = [] } = names;
  $('check-error').hidden = true;
  const groups = mappingGroups(alleles, genes);
  renderEvidencePlan(groups, genes, sources);
  const of = (state) => groups.filter((g) => g.state === state);
  const allelesOf = (list) => list.reduce((n, g) => n + g.alleles.length, 0);
  const open = of('none');
  $('check-summary').textContent = `${plural(alleles.length, 'allele')} with ${plural(new Set(alleles.map((a) => a.named)).size, 'gene name')}. `
    + (open.length ? `${plural(new Set(open.map((g) => g.named)).size, 'name')} (${plural(allelesOf(open), 'allele')}) not matched to a gene yet. ` : 'Every name is matched to a gene. ')
    + 'Nothing has been downloaded yet.';

  const tally = ['none', 'sequence', 'out', 'yours', 'rule', 'same'].filter((state) => of(state).length)
    .map((state) => el('span', { class: `map-tally t-${state}` }, el('b', { text: String(allelesOf(of(state))) }), ` ${LINKS[state].label}`));
  const head = el('thead', {}, el('tr', {}, ...['Name in your file', '', 'Gene of the list', 'Change'].map((h) => el('th', { scope: 'col', text: h }))));
  const shown = groups.filter((g) => g.state !== 'same');
  const same = of('same');
  $('check-map').replaceChildren(el('section', { class: 'map', 'aria-label': 'How the names of your file match the gene list' },
    el('h3', { class: 'plan-title', text: 'How the names of your file match the gene list' }),
    el('p', { class: 'map-tallies' }, 'Alleles: ', tally),
    of('rule').length ? el('p', { class: 'map-note' }, el('b', { text: 'Some names were redirected for you. ' }),
      `${[...new Set(of('rule').map((g) => g.named))].join(', ')}: the gene list has no gene of ${new Set(of('rule').map((g) => g.named)).size === 1 ? 'this name' : 'these names'}, but it knows the alleles under other gene names, shown below. Check that this is what you mean; you can redirect a line or correct the names.`) : null,
    shown.length ? el('div', { class: 'scroll' }, el('table', { class: 'map-table' }, head, el('tbody', {}, shown.flatMap((g) => mappingRow(g, sources))))) : null,
    open.length ? el('p', { class: 'map-open' }, `${plural(allelesOf(open), 'allele')} not matched: without a choice they are reported as "not placed". `,
      el('button', { type: 'button', class: 'ghost', text: 'Find all of them by sequence', onclick: () => { for (const g of open) decide(g, BY_SEQUENCE); } })) : null,
    same.length ? el('details', { class: 'fold', open: shown.length ? null : true }, el('summary', { text: `${plural(same.length, 'name')} found in the gene list as written (${plural(allelesOf(same), 'allele')})` }),
      el('div', { class: 'scroll map-same' }, el('table', { class: 'map-table' }, head.cloneNode(true), el('tbody', {}, same.flatMap((g) => mappingRow(g, sources)))))) : null));
  updateNameEditor();
  $('check-run').disabled = false;
  $('check-run').textContent = open.length ? `Run the analysis (${plural(allelesOf(open), 'allele')} not matched)` : 'Run the analysis';
}

/** Start the run on the file that was checked. */
function runChecked() {
  const { name } = names.pending;
  // Alleles the user chose to leave out are taken from the FASTA text.
  const skipped = new Set(names.alleles.filter((a) => decision(a) === LEFT_OUT).map((a) => a.name));
  const fasta = currentFasta();
  const kept = skipped.size
    ? fasta.split(/^(?=>)/m).filter((record) => !skipped.has(record.slice(1).split(/\s/)[0])).join('')
    : fasta;
  if (!kept.trim()) return checkFailed('Every allele of the file is left out: there is nothing to analyse.');
  $('progress-name').textContent = name;
  startProgress();
  show('progress');
  // Ask the browser not to evict the stored data when disk space is short.
  if (navigator.storage && navigator.storage.persist) navigator.storage.persist().catch(() => {});
  worker.postMessage({
    fasta: kept,
    name,
    build: $('opt-build').value,
    marker: $('opt-marker').value || '_S',
    haplotypes: $('opt-hap').checked,
    gnomad: $('opt-gnomad').checked,
    paralog: $('opt-paralog').checked,
    coords: names.rows ? geneListText(names.rows) : undefined,
    aliases: aliasText(),
  });
}

// ---- the gene list, shown and edited ---------------------------------------------

function renderGeneList() {
  const rows = geneRows();
  $('gene-editor-state').textContent = names.rows ? `your own, ${plural(rows.length, 'gene')}` : `built-in, ${plural(rows.length, 'gene')} of the human IGH locus`;
  $('gene-reset').hidden = !names.rows;
  const wanted = $('gene-filter').value.trim().toLowerCase();
  const shown = rows.map((row, i) => ({ row, i })).filter(({ row }) => !wanted || row.gene.toLowerCase().includes(wanted))
    .sort((a, b) => sortGenes(a.row.gene, b.row.gene) || a.row.build.localeCompare(b.row.build));
  $('gene-table').replaceChildren(
    el('thead', {}, el('tr', {}, ...['Gene', 'Build', 'Contig', 'First base', 'Last base', 'Strand', ''].map((h) => el('th', { scope: 'col', text: h })))),
    el('tbody', {}, shown.map(({ row, i }) => el('tr', {},
      el('th', { scope: 'row', text: row.gene }),
      el('td', { text: row.build }), el('td', { text: row.contig }),
      el('td', { class: 'num', text: row.start.toLocaleString('en') }), el('td', { class: 'num', text: row.end.toLocaleString('en') }),
      el('td', { text: row.strand === '.' ? 'either' : row.strand }),
      el('td', {}, el('button', { type: 'button', class: 'link', 'aria-label': `Remove ${row.gene} on ${row.build}`, text: 'remove', onclick: () => setGeneList(rows.filter((_, k) => k !== i)) }))))));
}

/** Use `rows` as the gene list (null = the built-in one) and check again. */
function setGeneList(rows) {
  names.rows = rows;
  keepNames();
  renderGeneList();
  askCheck();
}

function wireGeneList() {
  $('check-run').addEventListener('click', runChecked);
  $('example').addEventListener('click', async () => {
    try {
      const response = await fetch(new URL(`./example_KIARVA_IGHV.fasta${new URL(import.meta.url).search}`, import.meta.url));
      if (!response.ok) throw new Error(`status ${response.status}`);
      await handleFiles([new File([await response.text()], 'KIARVA_IGHV.fasta', { type: 'text/plain' })]);
    } catch (error) {
      fail(`The example could not be loaded (${error.message}).`);
    }
  });
  for (const id of ['opt-build', 'opt-hap', 'opt-gnomad']) $(id).addEventListener('change', askCheck);
  $('gene-filter').addEventListener('input', renderGeneList);
  $('gene-reset').addEventListener('click', () => setGeneList(null));
  $('gene-upload').addEventListener('click', () => $('gene-file').click());
  $('gene-file').addEventListener('change', async () => {
    const file = $('gene-file').files[0];
    $('gene-file').value = '';
    if (!file) return;
    try {
      setGeneList(parseGeneList(await file.text()));
      $('check-error').hidden = true;
    } catch (error) {
      $('check-error').textContent = `${file.name}: ${error.message}`;
      $('check-error').hidden = false;
    }
  });
  $('gene-download').addEventListener('click', () => {
    const link = el('a', { href: URL.createObjectURL(new Blob([`${geneListText(geneRows())}\n`], { type: 'text/tab-separated-values' })), download: 'igpop_gene_list.tsv' });
    link.click();
    URL.revokeObjectURL(link.href);
  });
  $('gene-add').addEventListener('submit', (event) => {
    event.preventDefault();
    const form = new FormData(event.target);
    const row = { gene: form.get('gene').trim(), build: form.get('build'), contig: form.get('contig').trim(), start: Number(form.get('start')), end: Number(form.get('end')), strand: form.get('strand') };
    const rows = geneRows();
    if (row.end < row.start) {
      $('check-error').textContent = 'The last base of a gene cannot lie before its first base.';
      $('check-error').hidden = false;
      return;
    }
    if (rows.some((r) => r.gene === row.gene && r.build === row.build)) {
      $('check-error').textContent = `${row.gene} is already in the list for ${row.build}. Remove it first to change it.`;
      $('check-error').hidden = false;
      return;
    }
    $('check-error').hidden = true;
    setGeneList([...rows, row]);
    event.target.reset();
  });
}

// ---- results ---------------------------------------------------------------

function display(files, warnings, summary = []) {
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
  renderPriorNote();
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

// Candidates are named by the marker, as in the command-line program.
const isCandidate = (allele) => allele.allele_name.includes(state.marker);

function badge(verdict) {
  const v = VERDICTS[verdict] || { label: verdict, icon: '' };
  return el('span', { class: `badge v-${verdict}`, title: v.help }, el('i', { 'aria-hidden': 'true', text: v.icon }), v.label);
}

// ---- the limit for "supported": a setting, not part of the evidence ---------

const LIMITS = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100];
const RUN_LIMIT = 10; // the limit the analysis itself was run with
const CONFIRM = 0.1;
const FILTER_ORDER = ['no_population_data', 'never_observed', 'single_observation', 'rare_in_every_population', 'common_not_demonstrated', 'low_score', 'intermediate_score', 'not_scored'];
const SCORE_FILTERS = ['low_score', 'intermediate_score', 'not_scored'];

/** Keep what the analysis decided, so that the limit can be moved and
 *  moved back. */
function rememberVerdicts(files) {
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
        const shown = given(v.run.filters) ? v.run.filters.split(',').filter((f) => !SCORE_FILTERS.includes(f)) : ['common_not_demonstrated'];
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
function renderLimit() {
  const model = scoreModel();
  $('limit-setting').hidden = !model;
  if (!model) return;
  $('limit').value = String(LIMITS.indexOf(state.limit));
  $('limit-value').textContent = String(state.limit);
  $('limit-reset').hidden = state.limit === RUN_LIMIT;
  const at = model.supportAt(state.limit);
  const supported = state.alleles.filter((a) => a.verdict === 'supported').length;
  const before = state.alleles.filter((a) => a.run.verdict === 'supported').length;
  const change = supported === before ? '' : ` (${supported > before ? `${supported - before} more` : `${before - supported} fewer`} than at ${RUN_LIMIT})`;
  $('limit-note').textContent = Number.isFinite(at)
    ? `A score of ${state.limit} is reached at a support of ${percent(at)}. ${pct(model.chance(at))} of all possible single-base changes in these genes score this high, which is how often an arbitrary change would pass. ${supported} alleles are supported${change}. ${state.limit === RUN_LIMIT ? '10 is the usual mark for strong evidence; the scores themselves do not depend on it.' : 'Verdicts, counts and the tables under Files follow this limit; the report describes the run at 10.'}`
    : `No support reaches a score of ${state.limit} in this run. ${supported} alleles are supported${change}.`;
}

function wireLimit() {
  const move = (limit) => {
    applyLimit(limit);
    renderTiles();
    renderLimit();
    renderAlleles();
    renderMethod();
    renderFiles();
  };
  $('limit').addEventListener('input', () => {
    // While dragging, only the number moves; the list follows on release.
    $('limit-value').textContent = String(LIMITS[Number($('limit').value)]);
  });
  $('limit').addEventListener('change', () => move(LIMITS[Number($('limit').value)]));
  $('limit-reset').addEventListener('click', () => move(RUN_LIMIT));
  $('f-order').addEventListener('input', renderAlleles);
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

function renderTiles() {
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

/** The key to the frequency lines, shown with each table of differences. */
function frequencyLegend() {
  return el('div', { class: 'legend' },
    el('span', { class: 'key' }, mark('kgp', 6, 6), ' 1000 Genomes population'),
    el('span', { class: 'key' }, mark('gnomad', 6, 6), ' gnomAD group'),
    el('span', { class: 'key' }, mark('dbsnp', 6, 6), ' dataset in dbSNP'),
    el('span', { class: 'key' }, el('i', { class: 'line' }), ' 1% line'));
}

/** What the folded controls are set to, so that nothing set is hidden. */
function renderFolds() {
  const active = ['f-kind', 'f-hap', 'f-place'].filter((id) => $(id).value).length;
  $('fold-filters-state').textContent = active ? `: ${active} set` : '';
  $('fold-reading-state').textContent = `: supported from a score of ${state.limit}${state.prior ? `, probability shown` : ''}`;
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
    if (gene && a.gene !== gene) return false;
    if (kind === 'candidate' && !isCandidate(a)) return false;
    if (kind === 'trusted' && isCandidate(a)) return false;
    if (hap && a.hap_status !== hap) return false;
    if (query && !`${a.allele_name} ${a.gene} ${a.reason}`.toLowerCase().includes(query)) return false;
    return true;
  });
}

/** The one setting that turns scores into probabilities, for the whole
 *  result: off unless the user turns it on. */
function wirePrior() {
  const [on, box, value, note] = [$('prior-on'), $('prior-box'), $('prior-value'), $('prior-note')];
  try {
    const kept = JSON.parse(localStorage.getItem('igpop-prior') || 'null');
    if (kept) { on.checked = Boolean(kept.on); value.value = String(kept.percent); }
  } catch (error) { /* storage is off */ }
  const apply = () => {
    const percentage = Math.min(99, Math.max(1, Number(value.value) || 50));
    box.hidden = !on.checked;
    state.prior = on.checked ? percentage / 100 : null;
    try { localStorage.setItem('igpop-prior', JSON.stringify({ on: on.checked, percent: percentage })); } catch (error) { /* storage is off */ }
    if (state.alleles.length) renderAlleles();
  };
  on.addEventListener('change', apply);
  value.addEventListener('change', apply);
  note.addEventListener('click', (event) => {
    const use = event.target.closest('button[data-percent]');
    if (!use) return;
    value.value = use.dataset.percent;
    apply();
  });
  box.hidden = !on.checked;
  state.prior = on.checked ? Math.min(99, Math.max(1, Number(value.value) || 50)) / 100 : null;
}

/** Something to hold on to when choosing the setting: how many candidates
 *  of this run turned out to exist on phased chromosomes. */
function renderPriorNote() {
  const note = $('prior-note');
  const looked = state.alleles.filter((a) => a.trusted === 'no' && ['found_exact', 'found_with_extra', 'not_found'].includes(a.hap_status));
  const found = looked.filter((a) => a.hap_status === 'found_exact').length;
  note.replaceChildren('This is your judgement of the method that produced the alleles; the score then moves it up or down.');
  if (looked.length >= 20) {
    const share = Math.min(99, Math.max(1, Math.round((100 * found) / looked.length)));
    note.append(` For orientation: of the ${looked.length} candidates in this run whose whole sequence could be looked up on phased chromosomes, ${found} were found there. `,
      el('button', { type: 'button', class: 'link', 'data-percent': String(share), text: `Use ${share}%` }));
  }
}

function renderAlleles() {
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

function alleleDetail(allele) {
  const parts = [el('p', { class: 'reason', text: allele.reason })];
  const placement = placementNote(allele);
  if (placement) parts.push(placement);
  if (given(allele.notes)) parts.push(el('p', { class: 'muted', text: `Keep in mind: ${allele.notes}.` }));
  const variants = state.variants.get(allele.allele_name) || [];
  // Genome position of each difference -> opens its evidence in the table.
  const opened = new Map();
  if (variants.length) parts.push(variantTable(variants, opened));
  const hap = haplotypePanel(allele);
  if (hap) parts.push(hap);
  const block = state.alignments.get(allele.gene);
  if (block) parts.push(alignmentPanel(block, allele.allele_name, opened));
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

/** The score explained, with a few points of this run's curve. */
function scoreHint() {
  const model = scoreModel();
  const points = model ? [[0, 'never seen in any dataset'], [1e-4, 'shown to reach 0.01% in a population'], [1e-3, 'shown to reach 0.1%'], [1e-2, 'shown to reach 1%'], [1e-1, 'shown to reach 10%']] : [];
  return el('span', {},
    el('b', { text: 'How to read the score' }),
    el('span', { text: 'The score says how much more typical the frequency of this base is of the variants of real alleles than of an arbitrary base change (what an error would give). It rises smoothly with the frequency the data guarantee:' }),
    points.length ? el('span', { class: 'hint-table' }, points.map(([support, words]) => el('span', {},
      el('b', { text: scoreText(model.score(support)) }),
      el('span', { text: `${words}: ${scoreWords(model.score(support))}` })))) : null,
    el('span', { text: `${state.limit} or more counts as supported (the limit you can set above the list), 0.1 or less as needing confirmation. These limits are conventions for strong evidence; around 1 the frequency cannot tell.` }));
}

// ---- variants: one row each, with a frequency strip -------------------------

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

// Frequencies are always percentages: one decimal from 1% up, two
// significant digits below (0.70%, 0.0074%).
const percent = (af) => (!(af > 0) ? '0%' : af >= 0.01 ? `${(100 * af).toFixed(1)}%` : `${Number((100 * af).toPrecision(2))}%`);

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
    const shape = m.kind === 'kgp' ? svg('circle', { cx: x, cy: mid, r: 4, class: 'm-kgp' })
      : m.kind === 'gnomad' ? svg('path', { d: `M${x} ${mid - 5}l5 5-5 5-5-5z`, class: 'm-gnomad' })
        : svg('rect', { x: x - 3.5, y: mid - 3.5, width: 7, height: 7, class: 'm-dbsnp' });
    shape.append(svg('title', {}, `${m.name}: ${percent(m.af)}`));
    picture.append(shape);
  }
  return el('div', { class: 'strip-cell' }, picture, el('small', { text: label }));
}

const GNOMAD_NAMES = { afr: 'African / African American', ami: 'Amish', amr: 'Admixed American', asj: 'Ashkenazi Jewish', eas: 'East Asian', fin: 'Finnish', mid: 'Middle Eastern', nfe: 'European (non-Finnish)', sas: 'South Asian' };
const sourceTitle = (name) => ({ '1kGP_30x': '1000 Genomes, 30x', '1kGP_phase3': '1000 Genomes, phase 3', dbSNP: 'dbSNP' }[name] || name.replace('gnomAD_v', 'gnomAD v'));
/** A stratum such as `1kGP_30x:AFR` in words. */
function stratumWords(stratum) {
  const [source, population] = stratum.split(':');
  if (source.startsWith('1kGP')) return `${POP_NAMES[population] || population} people of ${sourceTitle(source)}`;
  if (source.startsWith('gnomAD')) return `the ${GNOMAD_NAMES[population] || population} group of ${sourceTitle(source)}`;
  return `the study ${population}`;
}

/** Every population that looked at a variant, with its counts:
 *  [{ source, total: {ac, an}, note, rows: [{ name, ac, an, estimated }] }]. */
function populations(v) {
  const out = [];
  const kgp = POPS.map((pop) => {
    const [ac, an] = (v[`kgp_AC/AN_${pop}`] || '').split('/').map(Number);
    return { name: POP_NAMES[pop], ac, an };
  }).filter((r) => r.an > 0);
  if (kgp.length) {
    const [het, hom] = [number(v.kgp_het), number(v.kgp_hom)];
    out.push({
      source: v.build === 'GRCh37' ? '1000 Genomes, phase 3' : '1000 Genomes, 30x',
      total: { ac: number(v.kgp_AC), an: number(v.kgp_AN) },
      note: het >= 0 && hom >= 0 && (het || hom) ? `${het} ${het === 1 ? 'person' : 'people'} with one copy, ${hom} with two` : '',
      rows: kgp,
    });
  }
  if (given(v.gnomad_groups)) {
    out.push({
      source: v.build === 'GRCh37' ? 'gnomAD v2.1 genomes' : 'gnomAD v4.1 genomes',
      total: { ac: number(v.gnomad_AC), an: number(v.gnomad_AN) },
      note: number(v.gnomad_nhomalt) > 0 ? `${v.gnomad_nhomalt} people with two copies` : '',
      rows: v.gnomad_groups.split(';').map((entry) => {
        const [group, counts] = entry.split(':');
        const [ac, an] = counts.split('/').map(Number);
        return { name: GNOMAD_NAMES[group] || group, ac, an };
      }),
    });
  }
  if (given(v.dbsnp_freq)) {
    const sizes = new Map((state.files['dbsnp_datasets.tsv'] ? parseTSV(state.files['dbsnp_datasets.tsv']) : []).map((d) => [d.dataset, Number(d.implied_chromosomes)]));
    out.push({
      source: 'Studies collected by dbSNP',
      note: 'dbSNP gives frequencies only; the counts are worked out from the estimated size of each study',
      rows: v.dbsnp_freq.split(';').map((entry) => {
        const [name, af] = entry.split(':');
        const an = sizes.get(name) || 0;
        return { name, ac: an ? Math.round(Number(af) * an) : NaN, an, af: Number(af), estimated: true };
      }),
    });
  }
  return out;
}

/** The level of support a guaranteed frequency falls in, as the
 *  calibration table names it. */
function supportLevel(support) {
  if (!(support > 0)) return 'none';
  return ['<0.001%', '0.001-0.01%', '0.01-0.1%', '0.1-1%', '1-10%', '>=10%'][[1e-5, 1e-4, 1e-3, 1e-2, 1e-1].filter((edge) => support >= edge).length];
}

/** Why a difference has its verdict and score, step by step. */
function whyScore(v) {
  const steps = [];
  const support = number(v.support);
  if (!(support > 0)) steps.push('No population of any source has seen this base, so no frequency is guaranteed.');
  else if (given(v.support_from)) steps.push(`The strongest evidence comes from ${stratumWords(v.support_from)}. Its counts guarantee, with 95% confidence, a frequency of at least ${percent(support)} there. This guaranteed frequency is what igpop calls the support; it is lower than the observed frequency, the more so the fewer people were examined.`);
  else steps.push(`The counts guarantee, with 95% confidence, a frequency of at least ${percent(support)} in some population.`);
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
    ? `Apart from the score, a test shows that the frequency reaches 1% in ${stratumWords(v.class_basis)}${given(v.class_p) ? ` (p = ${Number(v.class_p).toPrecision(2)} after correction for the number of populations)` : ''}. A common variant is supported whatever its score.`
    : 'dbSNP marks this base as common. A common variant is supported whatever its score.');
  else if (support > 0) steps.push(`No population shows a frequency of 1% or more${given(v.class_p) ? ` (best test: p = ${Number(v.class_p).toPrecision(2)})` : ''}, so the verdict follows the score: ${state.limit} or more is supported, 0.1 or less needs confirmation, in between inconclusive.`);
  return steps;
}

const scoreText = (x) => (x >= 10 ? x.toFixed(0) : x >= 1 ? x.toFixed(1) : x >= 0.1 ? x.toFixed(2) : x.toPrecision(1));

/** The curve the scores of this run come from (score_model.tsv), with the
 *  real variants it was drawn through. Worked out once per result. */
function scoreModel() {
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
function comparisonChart(v) {
  const model = scoreModel();
  if (!model) return null;
  const support = v ? number(v.support) : NaN;
  const [W, H] = [700, 318];
  const [slot0, slot1, x0, x1] = [52, 104, 122, 688];
  const xOf = (s) => (s > 0 ? x0 + ((Math.log10(Math.min(1, Math.max(1e-6, s))) + 6) / 6) * (x1 - x0) : (slot0 + slot1) / 2);
  const [top, bottom, yLow, yHigh] = [26, 168, Math.log10(0.002), Math.log10(300)];
  const yOf = (score) => bottom - ((Math.log10(Math.min(300, Math.max(0.002, score))) - yLow) / (yHigh - yLow)) * (bottom - top);
  const picture = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'compare', role: 'img',
    'aria-label': v
      ? `The score rises smoothly with support. This difference has a support of ${support > 0 ? percent(support) : 'zero (never seen)'} and a score of ${scoreText(Number(v.score))}.`
      : 'The score as a smooth curve of support, with the real variants and all possible changes it was made from.' });
  const text = (x, y, cls, words, anchor = 'middle') => svg('text', { x, y, class: cls, 'text-anchor': anchor }, words);

  // Frame: the slot for "never seen", then the axis of support.
  picture.append(svg('rect', { x: slot0, y: top, width: slot1 - slot0, height: bottom - top, class: 'panel' }));
  picture.append(svg('rect', { x: x0, y: top, width: x1 - x0, height: bottom - top, class: 'panel' }));
  for (const tick of [0.01, 0.1, 1, 10, 100]) {
    picture.append(svg('line', { x1: slot0, x2: x1, y1: yOf(tick), y2: yOf(tick), class: tick === 0.1 ? 'limit low' : 'grid' }));
    picture.append(text(slot0 - 5, yOf(tick) + 3, 'tick', String(tick), 'end'));
  }
  picture.append(svg('line', { x1: slot0, x2: x1, y1: yOf(state.limit), y2: yOf(state.limit), class: 'limit high' }));
  picture.append(text(x0 + 6, yOf(state.limit) - 4, 'limit-name high', `supported from ${state.limit}`, 'start'));
  picture.append(text(x1 - 4, yOf(0.1) + 11, 'limit-name low', 'needs confirmation up to 0.1', 'end'));
  picture.append(text(slot0, top - 9, 'axis-name', 'score', 'start'));
  for (let decade = -6; decade <= 0; decade += 1) {
    const x = xOf(10 ** decade);
    picture.append(svg('line', { x1: x, x2: x, y1: top, y2: H - 30, class: 'grid' }));
    if (decade > -6) picture.append(text(x, H - 18, 'tick', `${Number((100 * 10 ** decade).toPrecision(1))}%`, decade === 0 ? 'end' : 'middle'));
  }
  picture.append(text((slot0 + slot1) / 2, H - 18, 'tick', 'never seen'));
  picture.append(text((x0 + x1) / 2, H - 4, 'axis-name', 'support: the frequency the data guarantee in the best population'));

  // The curve, with the range the real variants leave open.
  const samples = Array.from({ length: 121 }, (_, i) => 10 ** (-6 + i / 20));
  const band = [...samples.map((s) => `${xOf(s).toFixed(1)},${yOf(model.interval(s)[1]).toFixed(1)}`),
    ...samples.reverse().map((s) => `${xOf(s).toFixed(1)},${yOf(model.interval(s)[0]).toFixed(1)}`)];
  samples.reverse();
  picture.append(svg('polygon', { points: band.join(' '), class: 'band' }));
  picture.append(svg('polyline', { points: samples.map((s) => `${xOf(s).toFixed(1)},${yOf(model.score(s)).toFixed(1)}`).join(' '), class: 'curve' }));
  picture.append(svg('line', { x1: slot0 + 8, x2: slot1 - 8, y1: yOf(model.never), y2: yOf(model.never), class: 'curve' }));

  // The real variants: one dot each, at its own support.
  const dotBase = 222;
  picture.append(text(slot0 - 5, dotBase - 14, 'row-name', 'real', 'end'));
  picture.append(text(slot0 - 5, dotBase - 4, 'row-name', 'variants', 'end'));
  const stack = new Map();
  for (const variant of model.real) {
    const column = Math.round(xOf(variant.support) / 4.4);
    const n = stack.get(column) || 0;
    stack.set(column, n + 1);
    const dot = svg('circle', { cx: column * 4.4, cy: dotBase - 3 - (n % 10) * 4.2, r: 1.9, class: 'real' });
    dot.append(svg('title', {}, `${variant.label}: at least ${percent(variant.support)}, score ${scoreText(model.score(variant.support))}`));
    picture.append(dot);
  }
  picture.append(svg('line', { x1: slot0, x2: x1, y1: dotBase, y2: dotBase, class: 'base' }));

  // All possible changes: their share in each eighth of a step of ten.
  const barBase = H - 32;
  const shares = model.tail.map((t, k) => t.share - (model.tail[k + 1] ? model.tail[k + 1].share : 0));
  const tallest = Math.max(...shares, 1e-9);
  picture.append(text(slot0 - 5, barBase - 14, 'row-name', 'possible', 'end'));
  picture.append(text(slot0 - 5, barBase - 4, 'row-name', 'changes', 'end'));
  const step = (x1 - x0) / 48;
  shares.forEach((share, k) => {
    if (!(share > 0) || k === shares.length - 1) return;
    const tall = Math.max(1, (share / tallest) * 40);
    const bar = svg('rect', { x: xOf(model.tail[k].support) + 0.5, y: barBase - tall, width: step - 1, height: tall, class: 'arbitrary' });
    bar.append(svg('title', {}, `${pct(share)} of all possible changes have a support from ${percent(model.tail[k].support)}`));
    picture.append(bar);
  });
  picture.append(svg('rect', { x: slot0 + 10, y: barBase - 40, width: slot1 - slot0 - 20, height: 40, class: 'arbitrary' }));
  picture.append(text((slot0 + slot1) / 2, barBase - 44, 'count', pct(1 - model.seenShare)));

  if (v) {
    const at = xOf(support);
    picture.append(svg('line', { x1: at, x2: at, y1: top, y2: barBase, class: 'this' }));
    picture.append(svg('circle', { cx: at, cy: yOf(Number(v.score)), r: 4.5, class: 'this-dot' }));
    picture.append(text(Math.min(Math.max(at, 70), W - 70), top - 9, 'this-label', `this difference: score ${scoreText(Number(v.score))}`));
  }
  return picture;
}

const pct = (x) => (x >= 0.995 ? '100%' : `${Number((100 * x).toPrecision(2))}%`);
/** A probability in words: never a flat 0% or 100%. */
const chanceText = (p) => (p > 0.99 ? 'more than 99%' : p < 0.01 ? 'less than 1%' : `${Math.round(100 * p)}%`);

/** The probability that a difference is real, if a share `prior` of such
 *  differences is real before looking at population data. */
const afterEvidence = (prior, score) => (prior * score) / (prior * score + (1 - prior));

/** Everything behind one difference: where it is, who carries it, and why
 *  it has its score. */
function evidencePanel(v) {
  const groups = populations(v);
  const top = Math.max(1e-12, ...groups.flatMap((g) => g.rows.map((r) => (r.af ?? (r.an ? r.ac / r.an : 0)))));
  const bar = (r) => {
    const af = r.af ?? (r.an ? r.ac / r.an : 0);
    return el('tr', { class: af > 0 ? '' : 'zero' },
      el('th', { scope: 'row', text: r.name }),
      el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * af) / top}%` })),
      el('td', { class: 'num', text: Number.isFinite(r.ac) && r.an ? `${r.estimated ? 'about ' : ''}${r.ac.toLocaleString('en')} of ${r.an.toLocaleString('en')}` : '' }),
      el('td', { class: 'num', text: percent(af) }));
  };
  const kind = { SNV: 'single-base change', INS: 'insertion', DEL: 'deletion' }[v.type] || v.type;
  return el('div', { class: 'evidence' },
    el('div', { class: 'evidence-part' },
      el('h4', { text: 'Where' }),
      el('p', {}, `Base ${v.aln_pos} of your sequence (counted from its first base), ${v.build || ''} ${v.contig}:${Number(v.genome_pos).toLocaleString('en')} on the genome. `,
        `A ${kind}: the reference genome has ${v.vcf_ref}, this allele has ${v.vcf_alt.split(',').length > 1 ? v.query_variant : v.vcf_alt}`,
        el('span', { class: 'muted', text: ' (written on the forward strand of the genome, as in the population files; the alignment below shows the gene in its own direction).' }))),
    el('div', { class: 'evidence-part' },
      el('h4', { text: 'Who carries it' }),
      groups.length ? groups.map((g) => el('table', { class: 'people' },
        el('caption', {}, el('b', { text: g.source }), g.total && g.total.an > 0 ? `: ${g.total.ac.toLocaleString('en')} of ${g.total.an.toLocaleString('en')} chromosomes carry it` : '',
          g.note ? el('small', { class: 'muted', text: ` ${g.note}` }) : null),
        el('thead', { class: 'sr' }, el('tr', {}, ...['Population', 'Share', 'Chromosomes with it', 'Frequency'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, g.rows.sort((a, b) => (b.af ?? b.ac / b.an) - (a.af ?? a.ac / a.an)).map(bar))))
        : el('p', { class: 'muted', text: 'No source has a record of this base at this position.' }),
      groups.length ? el('p', { class: 'muted', text: 'Each person has two chromosomes, so the counts are of chromosomes. Bars are drawn relative to the population where the base is most frequent.' }) : null),
    el('div', { class: 'evidence-part' },
      el('h4', { text: 'Why this verdict and score' }),
      el('ol', {}, whyScore(v).map((step) => el('li', { text: step })))),
    comparison(v));
}

/** The part of the evidence panel that shows what a single-base difference
 *  was compared with, and what its score is worth. */
function comparison(v) {
  if (!given(v.score)) return null;
  const model = scoreModel();
  if (!model) return null;
  const score = Number(v.score);
  const [low, high] = given(v.score_low) ? [Number(v.score_low), Number(v.score_high)] : model.interval(number(v.support));
  const reading = score >= state.limit ? (state.limit >= 10 ? 'That is strong evidence for a real variant.' : 'That is above the limit you set for a supported difference.')
    : score > 1 ? 'That speaks for a real variant, the more so the higher the score, but is not strong evidence by itself. "Inconclusive" means exactly that, not "ruled out".'
      : score > 0.1 ? 'That speaks mildly against a real variant without ruling it out.'
        : 'That is what an arbitrary change looks like: the difference needs confirmation from other data.';
  return el('div', { class: 'evidence-part' },
    el('h4', { text: 'What it was compared with' }),
    el('p', {}, model.own
      ? `Blue dots: the ${model.real.length} distinct single-base variants of the ${model.alleles} alleles (${plural(model.genes, 'gene')}) that this run accepts as real, namely the database alleles and the candidates whose whole sequence was found on several phased chromosomes. Point at a dot to see which variant it is. Grey bars: all ${model.counts.arbitrary.toLocaleString('en')} single-base changes that are possible in the genes of the run, which is what an arbitrary error looks like. The line is the score: how much more typical a support is of the dots than of the bars. It is fitted to all of them at once, so it has no steps.`
      : `The scores of this run come from igpop's built-in curve, because the run has too few established alleles to make its own. The curve was made from ${model.counts.trusted} real variants (blue dots) and ${model.counts.arbitrary.toLocaleString('en')} possible changes (grey bars) in human IGHV genes.`),
    el('div', { class: 'scroll' }, comparisonChart(v)),
    v.score_basis === 'borrowed' ? el('p', { text: 'This difference is an insertion or deletion, and the picture is made from single-base changes. Its score is read from the same line as a lower bound.' }) : null,
    el('p', {}, `This difference has a score of ${scoreText(score)}. The shaded band shows how well the real variants pin the line down here: between ${scoreText(low)} and ${scoreText(high)}. ${reading}`),
    state.prior
      ? el('p', {}, `With your setting that ${Math.round(100 * state.prior)}% of such differences are real before looking at population data, this one is real with a probability of `, el('b', { text: chanceText(afterEvidence(state.prior, score)) }), ` (between ${chanceText(afterEvidence(state.prior, low))} and ${chanceText(afterEvidence(state.prior, high))}, given how well the score is known).`)
      : null,
    el('p', { class: 'muted', text: 'The limits for the verdict (10 for "supported" unless you move it above the list of alleles, 0.1 for "needs confirmation") are conventions for strong evidence, as in other fields that use likelihood ratios; they are not measured. Because the score is smooth, a difference just below a limit is nearly as well supported as one just above it: read the number, not only the verdict.' }));
}

/** The population data of one difference, briefly, for a popup: per
 *  source the chromosomes that carry it, and the populations where it is
 *  most frequent as bars on the scale of the frequency line. */
function peoplePopup(v) {
  const groups = populations(v);
  const support = number(v.support);
  const share = (r) => r.af ?? (r.an ? r.ac / r.an : 0);
  // The same scale as the frequency line: one in a million to 100%.
  const width = (af) => `${Math.max(2, Math.min(100, ((Math.log10(af) - LOG_MIN) / -LOG_MIN) * 100))}%`;
  const kind = (source) => (/1000 Genomes/.test(source) ? 'kgp' : /gnomAD/.test(source) ? 'gnomad' : 'dbsnp');
  const blocks = groups.map((g) => {
    const seen = [...g.rows].filter((r) => share(r) > 0).sort((x, y) => share(y) - share(x));
    const top = seen.slice(0, 3);
    return el('span', { class: `tip-source s-${kind(g.source)}` },
      el('span', { class: 'tip-head' }, el('b', { text: g.source }),
        g.total && g.total.an > 0 ? el('span', { text: `${g.total.ac.toLocaleString('en')} of ${g.total.an.toLocaleString('en')}` }) : null),
      top.length ? el('span', { class: 'tip-rows' }, top.flatMap((r) => [
        el('span', { class: 'tip-name', text: r.name }),
        el('span', { class: 'tip-track' }, el('i', { style: `width:${width(share(r))}` })),
        el('span', { class: 'tip-count', text: Number.isFinite(r.ac) && r.an ? `${r.estimated ? '~' : ''}${r.ac.toLocaleString('en')} of ${r.an.toLocaleString('en')}` : '' }),
        el('span', { class: 'tip-share', text: percent(share(r)) }),
      ])) : el('span', { class: 'tip-none', text: 'not seen' }),
      seen.length > top.length ? el('span', { class: 'tip-none', text: `and ${plural(seen.length - top.length, 'more population')}` }) : null);
  });
  return el('span', { class: 'people-tip' },
    blocks.length ? blocks : el('span', { text: 'No source has a record of this base at this position.' }),
    blocks.length ? el('span', { class: 'tip-none', text: 'Counts are chromosomes. Bars run from 0.0001% to 100% in steps of ten; the line marks 1%.' }) : null,
    support > 0 || given(v.score) ? el('span', { class: 'tip-foot' },
      support > 0 ? el('span', {}, 'Guaranteed frequency ', el('b', { text: `at least ${percent(support)}` })) : el('span', { text: 'Never seen' }),
      given(v.score) ? el('span', {}, 'Score ', el('b', { text: scoreText(Number(v.score)) })) : null) : null);
}

function variantTable(variants, opened) {
  const rsLink = (v) => (given(v.rsid)
    ? el('a', { href: `https://www.ncbi.nlm.nih.gov/snp/${v.rsid}`, target: '_blank', rel: 'noopener', text: v.rsid })
    : el('span', { class: 'muted', text: 'no rsID' }));
  const kind = { SNV: 'base change to', INS: 'insertion,', DEL: 'deletion, leaving' };
  const change = (v) => `${kind[v.type] || v.type} ${v.query_variant}`;
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
      el('td', { text: v.aln_pos }),
      el('td', {}, badge(v.verdict)),
      el('td', {}, chips(v.filters), siteQuality(v)),
      scoreCell(v),
      el('td', {}, explain(el('span', { class: 'hint-value plain', tabindex: '0' }, strip(v)), peoplePopup(v))),
      el('td', {}, button)), detail];
  });
  return el('div', { class: 'scroll' }, el('table', { class: 'variants' },
    el('caption', { text: 'Differences from the reference genome' }),
    el('thead', {}, el('tr', {},
      ...['Variant', 'Position in allele', 'Verdict by frequency', 'Why'].map((h) => el('th', { scope: 'col', text: h })),
      el('th', { scope: 'col' }, 'Score ', hint('What the score means', scoreHint())),
      el('th', { scope: 'col' }, 'How common in people ', hint('How to read the frequency line', stripHint())),
      el('th', { scope: 'col', text: 'Details' }))),
    el('tbody', {}, rows)), frequencyLegend());
}

// Why the 1000 Genomes genotypes at a position are less reliable: the
// project's own criteria for a site, in plain words.
const SITE_QC = [
  [/^FILTER=/, 'the project\u2019s variant filter did not pass this call: the reads that support it look less like a true variant than usual'],
  [/^missingness/, 'more than 5% of the people have no genotype here, usually because too few reads could be placed'],
  [/^HWE/, 'the numbers of carriers of one and of two copies do not fit together in any population, as happens when reads of a second copy of the gene are mixed in'],
  [/^mendelian/, 'in more than 5% of the families the child\u2019s genotype does not fit the parents\u2019'],
];

/** A note under the reasons when the position fails 1000 Genomes site checks. */
function siteQuality(v) {
  if (!given(v.kgp_site_qc) || v.kgp_site_qc === 'PASS') return null;
  const failed = v.kgp_site_qc.split(/[;,]/).filter(Boolean);
  const words = failed.map((item) => (SITE_QC.find(([pattern]) => pattern.test(item)) || [null, item])[1]);
  return el('span', { class: 'qc' }, explain(
    el('span', { class: 'chip chip-qc hint-value', tabindex: '0', text: 'less reliable position' }),
    el('span', {},
      el('b', { text: 'The 1000 Genomes genotypes are less reliable here' }),
      el('span', { text: 'The project checks every position of the genome. This one fails:' }),
      el('span', { class: 'hint-list' }, words.map((w) => el('span', { text: `\u2022 ${w}` }))),
      el('span', { text: 'This is common in antibody genes, which exist in several similar copies. It does not mean the variant is wrong. It means its frequency in 1000 Genomes is less certain, and that the haplotype check cannot use this position.' }))));
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
    given(v.chance) && number(v.support) > 0 ? el('span', { text: `${pct(Number(v.chance))} of all possible changes are supported at least as well.` }) : null)),
  state.prior ? el('small', { class: 'probability', text: `${chanceText(afterEvidence(state.prior, score)).replace('more than ', '>').replace('less than ', '<')} real` }) : null);
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

/** One gene block of alignments.txt as columns: the positions on the
 *  genome, and one string per row (reference first). */
function parseBlock(block) {
  const lines = block.split('\n');
  const where = lines[0].match(/^>(\S+)\s+(\S+)\s+(\S+):(\d+)-(\d+)/) || [];
  const out = { title: lines[0].slice(1), build: where[2] || '', contig: where[3] || '', positions: [], rows: new Map(), notes: new Map() };
  let [width, first] = [0, 0];
  for (const line of lines.slice(1)) {
    const chunk = line.match(/^\s+(\d+) \.\. (\d+)\s*$/);
    if (chunk) {
      const [from, to] = [Number(chunk[1]), Number(chunk[2])];
      width = Math.abs(to - from) + 1;
      for (let k = 0; k < width; k += 1) out.positions.push(from + (to >= from ? k : -k));
      continue;
    }
    const note = line.match(/^ {2}(\S+)\s+(identity .*)$/);
    if (note) {
      out.notes.set(note[1], note[2]);
      continue;
    }
    const row = line.match(/^(\S+)(\s+)/);
    if (!row || !width || line.startsWith(' ')) continue;
    // The reference row comes first and has no blanks: it shows in which
    // column of the text the sequences begin.
    if (row[1] === 'reference') first = row[0].length;
    const before = out.positions.length - width;
    const seq = (out.rows.get(row[1]) || '').padEnd(before, ' ') + line.slice(first).padEnd(width, ' ').slice(0, width);
    out.rows.set(row[1], seq);
  }
  return out;
}

const BASE_WORDS = { '.': 'same as the reference', '-': 'base missing (deletion)', ' ': 'not covered by this allele' };

/** The alignment of one gene: every allele under the reference, one column
 *  per base. Pointing at a column highlights it and says what is there. */
function alignmentPanel(block, alleleName, opened = new Map()) {
  const aln = parseBlock(block);
  const columns = aln.positions.length;
  const names = ['reference', ...(aln.rows.has(alleleName) ? [alleleName] : []), ...[...aln.rows.keys()].filter((n) => n !== 'reference' && n !== alleleName)];
  const ruler = Array.from({ length: columns }, (_, c) => ((c + 1) % 10 === 0 ? '|' : (c + 1) % 5 === 0 ? ':' : '·')).join('');
  const numbers = Array.from({ length: Math.floor(columns / 10) }, (_, k) => String(10 * (k + 1)).padStart(10, ' ')).join('');
  const line = (name, cls, content) => el('div', { class: `aln-line ${cls}`, 'data-name': name }, content);
  const coloured = (seq) => {
    const parts = [];
    for (const piece of seq.split(/([A-Za-z-])/)) {
      if (!piece) continue;
      parts.push(/^[A-Za-z-]$/.test(piece) ? el('b', { class: `nt nt-${piece === '-' ? 'gap' : piece.toUpperCase()}`, text: piece }) : piece.replace(/\./g, '·'));
    }
    return parts;
  };
  const marker = el('div', { class: 'aln-column', hidden: true });
  const seqs = el('div', { class: 'aln-seqs', style: `width:${columns}ch` }, marker,
    line('', 'aln-ruler', numbers), line('', 'aln-ruler', ruler),
    names.map((name) => line(name, name === 'reference' ? 'aln-ref' : name === alleleName ? 'current' : '',
      name === 'reference' ? aln.rows.get(name) : coloured(aln.rows.get(name) || ''))));
  const labels = el('div', { class: 'aln-names' },
    line('', 'aln-ruler', 'position in gene'), line('', 'aln-ruler', ' '),
    names.map((name) => line(name, name === 'reference' ? 'aln-ref' : name === alleleName ? 'current' : '', name)));
  const tip = el('div', { class: 'aln-tip', role: 'status', hidden: true });
  const grid = el('div', { class: 'aln-grid', tabindex: '0', role: 'group', 'aria-label': `Alignment of ${names.length - 1} alleles of this gene. Use the left and right arrow keys to move along it.` }, labels, seqs);
  const at = { col: -1, name: alleleName };

  const describe = (col, name) => {
    const pos = aln.positions[col];
    const reference = aln.rows.get('reference')[col];
    const row = name && name !== 'reference' && aln.rows.has(name) ? name : alleleName;
    const base = (aln.rows.get(row) || '')[col] || ' ';
    const differing = [...aln.rows].filter(([n, seq]) => n !== 'reference' && /[A-Za-z-]/.test(seq[col] || '')).length;
    const parts = [
      el('b', { text: `Base ${col + 1} of the gene` }),
      el('span', { class: 'muted', text: `${aln.build} ${aln.contig}:${pos.toLocaleString('en')}` }),
      el('span', {}, 'Reference: ', el('b', { class: `nt nt-${reference}`, text: reference })),
      el('span', {}, `${row}: `, BASE_WORDS[base] ? BASE_WORDS[base] : el('b', { class: `nt nt-${base.toUpperCase()}`, text: base })),
    ];
    // The difference of this allele at this position, if it has one.
    const variant = (state.variants.get(row) || []).find((v) => (v.type === 'SNV'
      ? Number(v.genome_pos) === pos
      : Math.abs(Number(v.genome_pos) - pos) <= Math.max(v.vcf_ref.length, v.vcf_alt.length)));
    if (variant && !BASE_WORDS[base] || (variant && base === '-')) {
      const seen = frequencies(variant);
      const best = seen.reduce((top, m) => (m.af > top.af ? m : top), { af: 0 });
      const counted = populations(variant).filter((g) => g.total && g.total.an > 0)
        .map((g) => `${g.total.ac.toLocaleString('en')} of ${g.total.an.toLocaleString('en')} chromosomes in ${g.source}`);
      parts.push(
        el('span', {}, given(variant.rsid) ? variant.rsid : 'no rsID', ' ', badge(variant.verdict)),
        el('span', { text: seen.length ? `Seen in ${seen.length} population${seen.length > 1 ? 's' : ''}, most often ${percent(best.af)} in ${best.name}.` : 'Not seen in any population.' }),
        counted.length ? el('span', { text: `Carried by ${counted.join('; ')}.` }) : null,
        given(variant.score) ? el('span', { text: `Score ${Number(variant.score) < 1 ? Number(variant.score).toFixed(3) : Number(variant.score).toFixed(1)}: ${scoreWords(Number(variant.score))}.` }) : null,
        row === alleleName && opened.has(Number(variant.genome_pos)) ? el('span', { class: 'aln-more', text: 'Click (or press Enter) for the populations and the reasons.' }) : null);
    }
    parts.push(el('span', { class: 'muted', text: differing ? `${differing} of ${aln.rows.size - 1} alleles differ from the reference here.` : 'No allele differs from the reference here.' }));
    return parts;
  };
  const point = (col, name, scroll) => {
    if (col < 0 || col >= columns) return;
    if (col === at.col && name === at.name) return;
    Object.assign(at, { col, name });
    marker.hidden = false;
    marker.style.left = `${col}ch`;
    for (const row of grid.querySelectorAll('.aln-line.over')) row.classList.remove('over');
    if (name) for (const row of grid.querySelectorAll(`.aln-line[data-name="${CSS.escape(name)}"]`)) row.classList.add('over');
    tip.replaceChildren(...describe(col, name).filter(Boolean));
    tip.hidden = false;
    if (scroll) marker.scrollIntoView({ block: 'nearest', inline: 'center' });
    const box = marker.getBoundingClientRect();
    const frame = grid.getBoundingClientRect();
    const left = Math.min(Math.max(8, box.left + box.width / 2 - tip.offsetWidth / 2), window.innerWidth - tip.offsetWidth - 8);
    tip.style.left = `${left}px`;
    tip.style.top = frame.bottom + 8 + tip.offsetHeight < window.innerHeight ? `${frame.bottom + 8}px` : `${Math.max(8, frame.top - 8 - tip.offsetHeight)}px`;
  };
  const leave = () => {
    Object.assign(at, { col: -1 });
    marker.hidden = true;
    tip.hidden = true;
    for (const row of grid.querySelectorAll('.aln-line.over')) row.classList.remove('over');
  };
  seqs.addEventListener('mousemove', (event) => {
    const box = seqs.getBoundingClientRect();
    const row = event.target.closest('.aln-line');
    point(Math.floor(((event.clientX - box.left) / box.width) * columns), row ? row.dataset.name : alleleName, false);
  });
  // A click on a difference of this allele opens its evidence in the table above.
  const openAt = (col) => {
    const pos = aln.positions[col];
    const variant = (state.variants.get(alleleName) || []).find((v) => (v.type === 'SNV'
      ? Number(v.genome_pos) === pos
      : Math.abs(Number(v.genome_pos) - pos) <= Math.max(v.vcf_ref.length, v.vcf_alt.length)));
    const open = variant && /[A-Za-z-]/.test((aln.rows.get(alleleName) || '')[col] || '') && opened.get(Number(variant.genome_pos));
    if (open) { leave(); open(); }
  };
  seqs.addEventListener('click', (event) => {
    const box = seqs.getBoundingClientRect();
    openAt(Math.floor(((event.clientX - box.left) / box.width) * columns));
  });
  seqs.addEventListener('mouseleave', () => { if (document.activeElement !== grid) leave(); });
  grid.addEventListener('blur', leave);
  grid.addEventListener('scroll', () => { if (document.activeElement !== grid) leave(); });
  grid.addEventListener('keydown', (event) => {
    const step = { ArrowRight: 1, ArrowLeft: -1 }[event.key];
    if (event.key === 'Escape') leave();
    if (event.key === 'Enter' && at.col >= 0) openAt(at.col);
    if (!step) return;
    event.preventDefault();
    // The arrow keys walk through the differences of this allele first.
    const own = aln.rows.get(alleleName) || '';
    let col = at.col < 0 ? (step > 0 ? -1 : columns) : at.col;
    if (event.shiftKey) col += step;
    else do { col += step; } while (col >= 0 && col < columns && !/[A-Za-z-]/.test(own[col] || ''));
    point(Math.min(Math.max(col, 0), columns - 1), alleleName, true);
  });
  const own = aln.notes.get(alleleName);
  return el('details', { class: 'aln' }, el('summary', { text: 'Alignment of this gene' }),
    el('p', { class: 'muted', text: 'Every allele of the gene under the reference, read 5′ to 3′. A dot is the same base as the reference, a coloured letter is a difference, a dash is a missing base. This allele is the highlighted row. Point at a column to see what is there; with the keyboard, the arrow keys jump between the differences of this allele (hold Shift to move one base).' }),
    own ? el('p', { class: 'muted', text: `${alleleName}: ${own.replace(/\s{2,}/g, ', ')}.` }) : null,
    el('p', { class: 'aln-where', text: aln.title.replace(/\s{2,}/g, ' · ') }),
    grid, tip);
}

// ---- score and checks ------------------------------------------------------

// ---- the data behind the results ---------------------------------------------

/** What each data source contributed to this run, from the result tables:
 *  how many of the differences it has seen and for how many it settled
 *  that they are common. */
function renderData() {
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
  const title = (name) => ({ '1kGP_30x': '1000 Genomes, 30x', '1kGP_phase3': '1000 Genomes, phase 3', dbSNP: 'Studies collected by dbSNP' }[name] || name.replace('gnomAD_v', 'gnomAD v').replace(/$/, name.startsWith('gnomAD') ? ' genomes' : ''));
  const sections = builds.map((build) => {
    const mine = all.filter((v) => v.build === build);
    const assessed = mine.filter((v) => v.class !== 'not_assessed');
    const sources = [...new Set(mine.flatMap((v) => (given(v.evidence) ? v.evidence.split(';') : [])))];
    const chromosomes = (name) => {
      const column = name.startsWith('1kGP') ? 'kgp_AN' : name.startsWith('gnomAD') ? 'gnomad_AN' : null;
      const most = column ? Math.max(0, ...mine.map((v) => number(v[column])).filter((x) => x > 0)) : 0;
      return most ? `${most.toLocaleString('en')} chromosomes` : '';
    };
    const rows = sources.map((name) => ({
      name,
      seen: mine.filter((v) => given(v.evidence) && v.evidence.split(';').includes(name)).length,
      common: mine.filter((v) => v.class === 'common' && given(v.class_basis) && v.class_basis.split(':')[0] === name).length,
    })).sort((a, b) => b.seen - a.seen);
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
        el('thead', {}, el('tr', {}, ...['Data source', 'Size', `Differences it has seen (of ${mine.length})`, 'Shown by it to be common'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, rows.map((r) => el('tr', {},
          el('th', { scope: 'row', text: title(r.name) }),
          el('td', { class: 'muted', text: r.name === 'dbSNP' ? 'many studies, see below' : chromosomes(r.name) }),
          el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * r.seen) / mine.length}%` }), ` ${r.seen}`),
          el('td', { class: 'num', text: String(r.common) })))),
        el('tfoot', {}, el('tr', {}, el('td', { colspan: '4', class: 'muted', text: `${never} of the differences ${never === 1 ? 'was' : 'were'} seen in no source${mine.length - assessed.length ? `; ${mine.length - assessed.length} could not be looked up at all` : ''}. A difference counts once per source; "shown to be common" names the source whose population proved a frequency of 1% or more.` }))))),
      panel ? el('p', {}, el('b', { text: 'Whole sequences: ' }), `the phased panel has ${panel.toLocaleString('en')} chromosomes. ${looked.length} of ${alleles.length} alleles could be looked up completely, and ${found.length} of those were found with their exact sequence.`) : null);
  });
  view.append(el('p', { text: 'Which data the verdicts rest on, and how much each source had to say about the alleles of this run.' }), ...sections);

  // The studies behind dbSNP, by how many of the differences they report.
  const studies = new Map();
  for (const v of all) {
    if (!given(v.dbsnp_freq)) continue;
    for (const entry of v.dbsnp_freq.split(';')) {
      const [name, af] = entry.split(':');
      if (Number(af) > 0) studies.set(name, (studies.get(name) || 0) + 1);
    }
  }
  const sizes = new Map((state.files['dbsnp_datasets.tsv'] ? parseTSV(state.files['dbsnp_datasets.tsv']) : []).map((d) => [d.dataset, d]));
  if (studies.size) {
    const top = Math.max(...studies.values());
    view.append(el('details', { class: 'fold' }, el('summary', { text: `The ${studies.size} studies behind dbSNP` }),
      el('p', { class: 'muted', text: 'dbSNP passes on frequencies from these studies without the number of people. igpop estimates each size from the smallest frequency the study reports; a larger study can demonstrate rarer variants.' }),
      el('div', { class: 'scroll' }, el('table', { class: 'data-table' },
        el('thead', {}, el('tr', {}, ...['Study', 'Estimated size', 'Differences it reports'].map((h) => el('th', { scope: 'col', text: h })))),
        el('tbody', {}, [...studies].sort((a, b) => b[1] - a[1]).map(([name, seen]) => el('tr', {},
          el('th', { scope: 'row', text: name }),
          el('td', { class: 'muted', text: sizes.has(name) && Number(sizes.get(name).implied_chromosomes) ? `${Number(sizes.get(name).implied_chromosomes).toLocaleString('en')} chromosomes` : '' }),
          el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * seen) / top}%` }), ` ${seen}`))))))));
  }
}

function renderMethod() {
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
      el('thead', {}, el('tr', {}, ...['Support: frequency guaranteed in the best population', 'Score', 'How well the score is known', 'Possible changes supported at least as well'].map((h) => el('th', { scope: 'col', text: h })))),
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
  for (const id of ['f-search', 'f-gene', 'f-kind', 'f-hap', 'f-place']) $(id).addEventListener('input', renderAlleles);
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
wireGeneList();
wirePrior();
wireLimit();
wireNameEditor();
// The worker is started with the page, so the analysis module is loaded
// and the stored data counted before the first file is dropped.
startWorker();

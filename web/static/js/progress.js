// A run in progress: its steps, its downloads per source, its log.

import { $, el } from './dom.js?v=BUILD';
import { post } from './engine.js?v=BUILD';
import { on } from './state.js?v=BUILD';

/** How much downloaded data this browser keeps, with a way to delete it. */
export function showStored(bytes, reading) {
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
    bytes && !reading ? el('button', { type: 'button', class: 'link', text: 'Delete', onclick: () => post({ type: 'clear' }) }) : '',
  );
}

// The stages of a run, in the order they happen, and the progress lines
// of the analysis module that announce each of them.
const STEPS = [
  { id: 'place', label: 'Placing alleles on their genes', now: 'Reading the reference sequence of each gene and aligning your alleles to it.', starts: /reference:|^aligning/ },
  { id: 'download', label: 'Fetching population data', now: 'Fetching the regions of your genes from each dataset.', starts: /^downloading|^progress:/ },
  { id: 'read', label: 'Reading frequencies and haplotypes', now: 'Looking up every difference in the population data.', starts: /^reading / },
  { id: 'judge', label: 'Judging the alleles', now: 'Scoring every difference and writing the results.', starts: /^score:/ },
];

export const run = { step: -1, started: 0, timer: null, sources: new Map(), titles: new Map(), store: '', kept: '' };

export function startProgress() {
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
  [/^reading (.+?): \S+$/, (m) => `${m[1][0].toUpperCase()}${m[1].slice(1)} read`],
  [/^(\d+) founders, (\d+) haplotypes/, (m) => `${Number(m[2]).toLocaleString('en')} chromosomes of ${Number(m[1]).toLocaleString('en')} unrelated people`],
  [/^(\d+) regions, (\d+) from cache/, (m) => `${m[1]} regions`],
  [/^score: calibrated from (.*)/, (m) => `Score calibrated from ${m[1].replace(/^this run: /, '')}`],
];

export function logLine(message) {
  // What is to come: so many regions of a source, at an address. A
  // source read from several files (one per chromosome) is one line.
  const planned = message.match(/^plan: (\d+) regions of (.+) at (\S+)$/);
  if (planned) {
    const name = `${planned[2][0].toUpperCase()}${planned[2].slice(1)}`;
    run.titles.set(planned[3], name);
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
    const source = counted[5] && run.sources.get(run.titles.get(counted[5]));
    if (source) source[counted[4]] += 1;
    if (run.step === 1) {
      setBar(done / total);
      $('progress-now').textContent = `${done} of ${total} regions ready: ${done - kept} downloaded, ${kept} already in this browser.`;
      renderSources();
    }
    return;
  }
  if (run.step !== 1) $('sources').hidden = true;
  const reading = text.match(/^reading (.+?): \S+$/);
  if (reading && run.step === 2) $('progress-now').textContent = `Reading ${reading[1]}.`;
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

// The clock of a run stops when the page moves on.
on('section', (section) => {
  if (section !== 'progress') clearInterval(run.timer);
});

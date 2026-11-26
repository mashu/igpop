// The page before a run: what the file holds, which data will judge it,
// and how its names match the gene list.

import { $, el, hint, plural, selectTab } from './dom.js?v=BUILD';
import { post } from './engine.js?v=BUILD';
import { LINKS, decisionSelect, linkOf } from './mapping.js?v=BUILD';
import { showNameEditor, updateNameEditor } from './name-editor.js?v=BUILD';
import { BY_SEQUENCE, LEFT_OUT, aliasText, askCheck, currentFasta, decision, geneListText, keepNames, names, short } from './names.js?v=BUILD';
import { startProgress } from './progress.js?v=BUILD';
import { show, state } from './state.js?v=BUILD';
import { sortGenes } from './words.js?v=BUILD';

export function checkFailed(message) {
  $('check-summary').textContent = 'The gene list could not be used.';
  $('check-error').textContent = message;
  $('check-error').hidden = false;
  $('check-run').disabled = true;
}

/** What the run will do with the file, at a glance: how many alleles and
 *  names, and which data judge how many of the alleles. */
function renderOverview(groups, genes, sources, alleles) {
  const used = (gene) => groups.filter((g) => g.gene === gene.gene && g.state !== 'out');
  const allelesOn = (list) => list.reduce((n, gene) => n + used(gene).reduce((k, g) => k + g.alleles.length, 0), 0);
  const builds = [...new Set(sources.map((s) => s.build))];
  const open = groups.filter((g) => g.state === 'none');
  const openAlleles = open.reduce((n, g) => n + g.alleles.length, 0);
  const apart = genes.filter((g) => g.build && !g.evidence && used(g).length);
  const parts = [
    ...builds.map((build, i) => ({ key: i ? 'older' : 'main', count: allelesOn(genes.filter((g) => g.build === build && g.evidence)),
      name: i ? `judged with the older data of ${build}` : `judged with the data of ${build}` })),
    { key: 'none', count: allelesOn(apart), name: 'no population data' },
    { key: 'open', count: openAlleles, name: 'name not matched to a gene' },
  ].filter((part) => part.count > 0);
  const total = parts.reduce((n, part) => n + part.count, 0) || 1;
  const namesCount = new Set(alleles.map((a) => a.named)).size;
  const fact = (count, label, sub, cls = '') => el('div', { class: `fact ${cls}` },
    el('b', { text: count }), el('span', { text: ` ${label}` }), sub ? el('small', { text: sub }) : null);
  const notes = [];
  if (open.length) {
    notes.push(el('p', { class: 'check-note warn' }, el('span', { class: 'plan-warn', 'aria-hidden': 'true', text: '!' }),
      ` ${plural(new Set(open.map((g) => g.named)).size, 'gene name')} (${plural(openAlleles, 'allele')}) not matched to a gene. Unmatched alleles are reported as "not placed". `,
      el('button', { type: 'button', class: 'link', text: 'Review the names', onclick: () => selectTab('tab-check-names') })));
  }
  if (apart.length) {
    notes.push(el('p', { class: 'check-note' }, el('b', { text: `${apart.map((g) => g.gene).join(', ')}: ` }), 'no population data ',
      hint('Why there are no data', el('span', { class: 'help-text' }, el('span', { text: `The reference has ${apart.length === 1 ? 'this gene' : 'these genes'} only on a piece of sequence that is not part of a chromosome. The alleles are compared with it, but there are no dependable population data there: an allele that differs from it will be "not assessed".` })))));
  }
  if (builds.length > 1 && parts.some((part) => part.key === 'older')) {
    notes.push(el('p', { class: 'check-note' }, el('b', { text: `${genes.filter((g) => g.build === builds[1] && g.evidence && used(g).length).map((g) => g.gene).join(', ')}: ` }), `older data of ${builds[1]} `,
      hint('Why older data', el('span', { class: 'help-text' }, el('span', { text: `These genes are missing from ${builds[0]}. They are judged with the older data of ${builds[1]}, so read their verdicts with more caution.` })))));
  }
  $('check-overview').replaceChildren(
    el('div', { class: 'facts-row' },
      fact(alleles.length.toLocaleString('en'), alleles.length === 1 ? 'allele' : 'alleles'),
      fact(String(namesCount), namesCount === 1 ? 'gene name' : 'gene names', open.length ? `${new Set(open.map((g) => g.named)).size} not matched` : 'all matched to a gene', open.length ? 'warn' : 'good'),
      fact(String(sources.filter((source) => source.build === builds[0]).length), 'datasets')),
    el('div', { class: 'judged' },
      el('h3', { class: 'plan-title', text: 'Which data will judge your alleles' }),
      el('div', { class: 'judged-bar', role: 'img', 'aria-label': parts.map((part) => `${part.count} ${part.name}`).join(', ') },
        parts.map((part) => el('i', { class: `j-${part.key}`, style: `flex:${Math.max(part.count, total / 60)}`, title: `${part.count} ${part.name}` }))),
      el('p', { class: 'legend' }, parts.map((part) => el('span', { class: 'key' }, el('i', { class: `swatch j-${part.key}` }), el('b', { text: ` ${part.count}` }), ` ${part.name}`)))),
    ...notes);
  $('tab-check-names').textContent = open.length ? `Gene names (${openAlleles} not matched)` : 'Gene names';
  if (open.length && !renderOverview.shown) selectTab('tab-check-names');
  renderOverview.shown = true;
}

/** The datasets the run will read, per reference build: who is in each
 *  and how many people, as bars on one scale. */
function renderEvidencePlan(groups, genes, sources) {
  const used = (gene) => groups.filter((g) => g.gene === gene.gene && g.state !== 'out');
  const builds = [...new Set(sources.map((s) => s.build))];
  // People on a scale of ten, from 1,000 to 100,000.
  const width = (people) => `${Math.max(4, Math.min(100, 50 * (Math.log10(people) - 3)))}%`;
  const tables = builds.map((build) => {
    const members = genes.filter((g) => g.build === build && g.evidence && used(g).length);
    const mine = sources.filter((s) => s.build === build);
    if (!members.length) return null;
    return el('table', { class: 'plan-table' },
      el('caption', {}, el('i', { class: `swatch j-${build === builds[0] ? 'main' : 'older'}` }), ' ', el('b', { text: `Genes on ${build}` }), ` ${plural(members.length, 'gene')}`),
      el('thead', { class: 'sr' }, el('tr', {}, ...['Dataset', 'Gives', 'People'].map((h) => el('th', { scope: 'col', text: h })))),
      el('tbody', {}, mine.map((s) => el('tr', {},
        el('th', { scope: 'row' }, s.title, ' ', hint(`About ${s.title}`, el('span', { class: 'help-text' }, el('span', { text: `${s.about[0].toUpperCase()}${s.about.slice(1)}.` })))),
        el('td', {}, el('span', { class: `plan-kind k-${s.gives}`, text: s.gives === 'haplotypes' ? 'whole sequences' : 'frequencies' })),
        el('td', { class: 'plan-size' }, s.people
          ? [el('span', { class: 'plan-track' }, el('i', { class: `k-${s.gives}`, style: `width:${width(s.people)}` })), el('span', { class: 'plan-people', text: `${s.people.toLocaleString('en')} people` })]
          : el('span', { class: 'plan-people', text: 'differs by position' }))))),
      mine.some((s) => s.gives === 'haplotypes') ? null : el('tfoot', {}, el('tr', {}, el('td', { colspan: '3', class: 'muted', text: 'The haplotype check is switched off.' }))));
  }).filter(Boolean);
  $('check-evidence').replaceChildren(tables.length ? el('section', { class: 'plan', 'aria-label': 'The datasets the run will read' },
    tables,
    el('p', { class: 'muted plan-foot', text: 'Only the regions of your genes are downloaded. Bars compare the number of people, from 1,000 to 100,000 in steps of ten.' })) : '');
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

/** The check of the file: how its names connect to the gene list. Called
 *  with a new answer of the analysis module, or without to draw again. */
export function renderCheck(result) {
  if (result) Object.assign(names, { alleles: result.alleles, genes: result.genes, sources: result.sources });
  const { alleles, genes, sources = [] } = names;
  $('check-error').hidden = true;
  const groups = mappingGroups(alleles, genes);
  renderEvidencePlan(groups, genes, sources);
  renderOverview(groups, genes, sources, alleles);
  const of = (state) => groups.filter((g) => g.state === state);
  const allelesOf = (list) => list.reduce((n, g) => n + g.alleles.length, 0);
  const open = of('none');
  $('check-summary').textContent = '';

  const tally = ['none', 'sequence', 'out', 'yours', 'rule', 'same'].filter((state) => of(state).length)
    .map((state) => el('span', { class: `map-tally t-${state}` }, el('b', { text: String(allelesOf(of(state))) }), ` ${LINKS[state].label}`));
  const head = el('thead', {}, el('tr', {}, ...['Name in your file', '', 'Gene of the list', 'Change'].map((h) => el('th', { scope: 'col', text: h }))));
  const shown = groups.filter((g) => g.state !== 'same');
  const same = of('same');
  $('check-map').replaceChildren(el('section', { class: 'map', 'aria-label': 'How the names of your file match the gene list' },
    el('p', { class: 'map-tallies' }, 'Alleles: ', tally),
    of('rule').length ? el('p', { class: 'map-note' }, el('b', { text: 'Some names were redirected for you. ' }),
      `${[...new Set(of('rule').map((g) => g.named))].join(', ')}: the gene list has no gene of ${new Set(of('rule').map((g) => g.named)).size === 1 ? 'this name' : 'these names'}, but it knows the alleles under other gene names, shown below. Check that this is what you mean; you can redirect a line or correct the names.`) : null,
    shown.length ? el('div', { class: 'scroll' }, el('table', { class: 'map-table' }, head, el('tbody', {}, shown.flatMap((g) => mappingRow(g, sources))))) : null,
    open.length ? el('p', { class: 'map-open' }, `${plural(allelesOf(open), 'allele')} not matched: without a choice they are reported as "not placed". `,
      el('button', { type: 'button', class: 'ghost', text: 'Find all of them by sequence', onclick: () => { for (const g of open) decide(g, BY_SEQUENCE); } })) : null,
    same.length ? el('details', { class: 'fold' }, el('summary', { text: `${plural(same.length, 'name')} found in the gene list as written (${plural(allelesOf(same), 'allele')})` }),
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
  post({
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

/** A new allele file: show the check page and ask what its names match. */
export function startCheck(fasta, name) {
  state.prefix = name.replace(/\.[^.]+$/, '');
  names.pending = { fasta, name };
  names.renames = {};
  renderOverview.shown = false;
  selectTab('tab-check-data');
  $('check-name').textContent = name;
  $('check-summary').textContent = 'Reading the gene names…';
  for (const id of ['check-evidence', 'check-map']) $(id).replaceChildren();
  $('check-run').disabled = true;
  show('check');
  askCheck();
}

export function wireCheck() {
  $('check-run').addEventListener('click', runChecked);
  for (const id of ['opt-build', 'opt-hap', 'opt-gnomad']) $(id).addEventListener('change', askCheck);
}

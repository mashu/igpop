// The names of the sequences, corrected one by one or all at once.

import { $, el, plural, selectTab } from './dom.js?v=BUILD';
import { decisionSelect, linkOf } from './mapping.js?v=BUILD';
import { askCheck, keepNames, names, original } from './names.js?v=BUILD';

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
export function updateNameEditor() {
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
export function showNameEditor(named) {
  selectTab('tab-check-names');
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

export function wireNameEditor() {
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


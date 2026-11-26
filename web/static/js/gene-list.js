// The gene list, shown and edited.

import { $, el, plural } from './dom.js?v=BUILD';
import { askCheck, geneListText, geneRows, keepNames, names, parseGeneList } from './names.js?v=BUILD';
import { sortGenes } from './words.js?v=BUILD';

export function renderGeneList() {
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

export function wireGeneList() {
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


// The gene names of the file in hand and the gene list they are matched
// to, with what this browser remembers of both.

import { $ } from './dom.js?v=BUILD';
import { post } from './engine.js?v=BUILD';

// The gene list and the names of this browser. `rows` is the user's own
// gene list (null = the built-in one); `mapping` says what to do with a
// gene name the list does not have: another gene's name, 'sequence' (find
// the gene by sequence) or 'skip' (leave its alleles out).
export const names = { builtin: { coords: '', aliases: '' }, rows: null, mapping: {}, renames: {}, open: new Set(), pending: null, alleles: [], genes: [], asked: 0 };

const KEPT = 'igpop-gene-names';
try {
  const kept = JSON.parse(localStorage.getItem(KEPT) || '{}');
  if (Array.isArray(kept.rows)) names.rows = kept.rows;
  if (kept.mapping && typeof kept.mapping === 'object') names.mapping = kept.mapping;
} catch { /* nothing kept, or storage is off */ }

export function keepNames() {
  try { localStorage.setItem(KEPT, JSON.stringify({ rows: names.rows, mapping: names.mapping })); } catch { /* storage is off */ }
}

const COLUMNS = ['gene', 'build', 'contig', 'start', 'end', 'strand'];

/** Rows of a gene list in text form; throws with the line that is wrong. */
export function parseGeneList(text) {
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

export const geneRows = () => names.rows || (names.builtin.coords ? parseGeneList(names.builtin.coords) : []);

export const geneListText = (rows) => [COLUMNS.join('\t'), ...rows.map((r) => COLUMNS.map((c) => r[c]).join('\t'))].join('\n');

// What the user can decide for a name: a gene of the list, BY_SEQUENCE
// (find the gene by sequence during the run) or LEFT_OUT.
export const BY_SEQUENCE = '?';

export const LEFT_OUT = 'skip';

/** The decision that applies to an allele: its own, else that of its name. */
export const decision = (allele) => names.mapping[allele.name] ?? names.mapping[allele.named];

/** The alias table of a run: the built-in one and the user's decisions. */
export const aliasText = () => [names.builtin.aliases, ...Object.entries(names.mapping)
  .filter(([, to]) => to !== LEFT_OUT).map(([from, to]) => `${from}\t${to}`)].join('\n');

/** The FASTA text with the names the user corrected. */
export function currentFasta() {
  const { fasta } = names.pending;
  if (!Object.keys(names.renames).length) return fasta;
  return fasta.replace(/^>(\S+)/gm, (line, name) => (names.renames[name] ? `>${names.renames[name]}` : line));
}

export function askCheck() {
  if (!names.pending) return;
  names.asked += 1;
  post({
    type: 'check', id: names.asked, fasta: currentFasta(), name: names.pending.name, build: $('opt-build').value,
    haplotypes: $('opt-hap').checked,
    gnomad: $('opt-gnomad').checked,
    coords: names.rows ? geneListText(names.rows) : undefined, aliases: aliasText(),
  });
}

export const short = (allele) => allele.name.slice(allele.named.length) || allele.name;

export const original = (name) => Object.keys(names.renames).find((from) => names.renames[from] === name);


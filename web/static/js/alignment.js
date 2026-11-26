// The alignment of a gene: every allele under the reference, with the
// regions of the gene and the verdict of each difference.

import { el, hint, plural } from './dom.js?v=BUILD';
import { frequencies, populations } from './populations.js?v=BUILD';
import { state } from './state.js?v=BUILD';
import { given, parseBlock } from './tsv.js?v=BUILD';
import { scoreWords } from './variants.js?v=BUILD';
import { BASE_WORDS, VERDICTS, badge, percent } from './words.js?v=BUILD';

/** The alignment of one gene: every allele under the reference, one column
 *  per base. Pointing at a column highlights it and says what is there. */
const CODONS = (() => {
  const bases = 'TCAG';
  const acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const table = {};
  [...bases].forEach((a, i) => [...bases].forEach((b, j) => [...bases].forEach((c, k) => { table[a + b + c] = acids[16 * i + 4 * j + k]; })));
  return table;
})();

/** The framework and complementarity-determining regions of a heavy-chain
 *  V gene, found from its three conserved residues (the two cysteines and
 *  the tryptophan), as [name, first column, last column]. Approximate: the
 *  exact borders are defined on an alignment of all V genes, which igpop
 *  does not make. Null if the residues are not where a V gene has them. */
function regionsOf(reference, first, last) {
  const dna = reference.slice(first, last + 1).toUpperCase();
  let protein = '';
  for (let i = 0; i + 3 <= dna.length; i += 3) protein += CODONS[dna.slice(i, i + 3)] || 'X';
  // The first cysteine is codon 22 of the gene in every human IGHV gene;
  // a reference allele can have lost it.
  const cys1 = protein.slice(17, 28).includes('C') ? protein.slice(17, 28).indexOf('C') + 17 : 21;
  const found = protein.slice(cys1 + 5).search(/W[A-Z]RQ/);
  const trp = found >= 0 ? cys1 + 5 + found : protein.indexOf('W', cys1 + 5);
  const cys2 = protein.lastIndexOf('C');
  if (trp < 0 || cys2 < trp + 45) return null;
  const parts = [['FR1', 0, cys1 + 3], ['CDR1', cys1 + 4, trp - 3], ['FR2', trp - 2, trp + 14], ['CDR2', trp + 15, cys2 - 38], ['FR3', cys2 - 37, cys2]];
  if (parts.some(([, from, to]) => to < from)) return null;
  const regions = parts.map(([name, from, to]) => [name, first + 3 * from, first + 3 * to + 2]);
  if (first + 3 * cys2 + 3 <= last) regions.push(['CDR3', first + 3 * cys2 + 3, last]);
  return regions;
}

export function alignmentPanel(block, alleleName, opened = new Map()) {
  const aln = parseBlock(block);
  const columns = aln.positions.length;
  const names = ['reference', ...(aln.rows.has(alleleName) ? [alleleName] : []), ...[...aln.rows.keys()].filter((n) => n !== 'reference' && n !== alleleName)];
  // Where the gene itself begins and ends among the columns: the block
  // shows a few bases on either side of it.
  const span = aln.title.match(/:(\d+)-(\d+)\s+strand\s+(\S)/);
  const ends = span ? (span[3] === '-' ? [Number(span[2]), Number(span[1])] : [Number(span[1]), Number(span[2])]) : null;
  const geneFirst = ends ? Math.max(0, aln.positions.indexOf(ends[0])) : 0;
  const geneLast = ends && aln.positions.indexOf(ends[1]) >= 0 ? aln.positions.indexOf(ends[1]) : columns - 1;
  const inGene = (c) => c - geneFirst + 1;
  const numbers = [];
  for (let n = 10; geneFirst + n - 1 <= geneLast; n += 10) numbers.push(el('span', { class: 'aln-number', style: `left:${geneFirst + n - 0.5}ch`, text: String(n) }));
  const regions = /^IGHV/.test(aln.title) ? regionsOf(aln.rows.get('reference') || '', geneFirst, geneLast) : null;
  const regionAt = (c) => (regions || []).find(([, from, to]) => c >= from && c <= to);
  const regionLine = regions ? [
    el('span', { style: `display:inline-block;width:${regions[0][1]}ch` }),
    ...regions.map(([name, from, to]) => el('span', { class: `aln-region ${name.startsWith('CDR') ? 'cdr' : 'fr'}`, title: name === 'CDR3' ? 'The first bases of CDR3; the rest comes from the D and J genes' : name, style: `width:${to - from + 1}ch` }, to - from + 1 >= name.length + 1 ? el('small', { text: name }) : null))] : null;
  // The verdict of each difference, by allele and column, so that the
  // ones that are not supported stand out.
  const verdictAt = (name) => {
    const marks = new Map();
    for (const v of state.variants.get(name) || []) {
      if (v.verdict === 'supported' || v.type === 'INS') continue;
      const pos = Number(v.genome_pos);
      const covered = v.type === 'SNV' ? [pos] : Array.from({ length: Math.max(0, v.vcf_ref.length - 1) }, (_, k) => pos + 1 + k);
      for (const p of covered) {
        const c = aln.positions.indexOf(p);
        if (c >= 0) marks.set(c, v.verdict);
      }
    }
    return marks;
  };
  const line = (name, cls, content) => el('div', { class: `aln-line ${cls}`, 'data-name': name }, content);
  const coloured = (seq, name) => {
    const marks = verdictAt(name);
    const parts = [];
    let plain = '';
    [...seq].forEach((piece, c) => {
      if (!/^[A-Za-z-]$/.test(piece)) { plain += piece === '.' ? '·' : piece; return; }
      if (plain) { parts.push(plain); plain = ''; }
      parts.push(el('b', { class: `nt nt-${piece === '-' ? 'gap' : piece.toUpperCase()}${marks.has(c) ? ` nv-${marks.get(c)}` : ''}`, text: piece }));
    });
    if (plain) parts.push(plain);
    return parts;
  };
  const marker = el('div', { class: 'aln-column', hidden: true });
  // A band over the columns of one region: shown while the pointer is on
  // its label, and kept when the label is clicked.
  const band = el('div', { class: 'aln-band', hidden: true });
  let pinned = null;
  const showBand = (region) => {
    band.hidden = !region;
    if (region) Object.assign(band.style, { left: `${region[1]}ch`, width: `${region[2] - region[1] + 1}ch` });
  };
  for (const [i, label] of (regionLine || []).slice(1).entries()) {
    const region = regions[i];
    Object.assign(label, { tabIndex: 0 });
    label.setAttribute('role', 'button');
    label.setAttribute('aria-pressed', 'false');
    label.addEventListener('mouseenter', () => showBand(region));
    label.addEventListener('mouseleave', () => showBand(pinned));
    const toggle = (event) => {
      event.stopPropagation();
      pinned = pinned === region ? null : region;
      for (const other of regionLine.slice(1)) other.setAttribute('aria-pressed', String(other === label && pinned === region));
      showBand(pinned);
    };
    label.addEventListener('click', toggle);
    label.addEventListener('keydown', (event) => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); toggle(event); } });
  }
  const seqs = el('div', { class: 'aln-seqs', style: `width:${columns}ch` }, marker, band,
    line('', 'aln-ruler aln-numbers', numbers),
    el('div', { class: 'aln-line aln-ruler aln-ticks', style: `--from:${geneFirst}ch;--length:${geneLast - geneFirst + 1}ch` }),
    regionLine ? line('', 'aln-ruler aln-regions', regionLine) : null,
    names.map((name) => line(name, name === 'reference' ? 'aln-ref' : name === alleleName ? 'current' : '',
      name === 'reference' ? aln.rows.get(name) : coloured(aln.rows.get(name) || '', name))));
  const labels = el('div', { class: 'aln-names' },
    line('', 'aln-ruler aln-numbers', 'position'), line('', 'aln-ruler aln-ticks', ' '),
    regionLine ? line('', 'aln-ruler aln-regions', 'region') : null,
    names.map((name) => {
      // Each allele carries its verdict as a small mark before its name.
      const allele = state.alleles.find((a) => a.allele_name === name);
      const verdict = allele ? VERDICTS[allele.verdict] : null;
      return line(name, name === 'reference' ? 'aln-ref' : name === alleleName ? 'current' : '', [
        el('i', { class: `aln-verdict${allele ? ` v-${allele.verdict}` : ''}`, title: verdict ? verdict.label : '', 'aria-hidden': 'true' }),
        name === 'reference' ? 'reference genome' : name,
        verdict ? el('span', { class: 'sr', text: `, ${verdict.label}` }) : null]);
    }));
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
      el('b', { text: col < geneFirst ? `${geneFirst - col} base${geneFirst - col > 1 ? 's' : ''} before the gene` : col > geneLast ? `${col - geneLast} base${col - geneLast > 1 ? 's' : ''} after the gene` : `Base ${inGene(col)} of the gene${regionAt(col) ? `, in ${regionAt(col)[0] === 'CDR3' ? 'the first bases of CDR3' : regionAt(col)[0]}` : ''}` }),
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
  const own = aln.notes.get(alleleName) || '';
  const identity = own.match(/identity ([\d.]+)/);
  const where = aln.title.match(/^(\S+)\s+(\S+)\s+(\S+):(\d+)-(\d+)\s+strand\s+(\S)/);
  const key = (cls, words) => el('span', { class: 'aln-key' }, el('b', { class: `nt nt-A ${cls}`, text: 'A' }), words);
  const help = el('span', {},
    el('b', { text: 'How to read the alignment' }),
    el('span', { text: 'Every allele of the gene under the reference genome, read 5′ to 3′. A dot is the same base as the reference, a coloured letter is a difference, a dash is a missing base, and an empty place is not covered by the allele.' }),
    el('span', { text: 'A difference that is not supported is marked: underlined when it is inconclusive or has no population data, framed when it needs confirmation. The mark before a name is the verdict of the whole allele.' }),
    el('span', { text: 'Point at a column to see what is there; click a difference of the opened allele for its evidence. With the keyboard, the arrow keys jump between the differences of the opened allele (hold Shift to move one base).' }),
    regions ? el('span', { text: 'FR and CDR are the framework and complementarity-determining regions. Their borders are placed from the conserved cysteines and tryptophan of a V gene, not from an IMGT alignment, and can be off by a codon or two. A V gene holds only the first bases of CDR3. Point at a region to see its columns; click it to keep them highlighted.' }) : null,
    aln.notes.size ? el('span', { text: 'Inserted bases have no column: they are listed in the table of differences.' }) : null);
  return el('details', { class: 'aln' }, el('summary', { text: 'Show the alignment of this gene' }),
    el('div', { class: 'aln-head' },
      el('span', { class: 'aln-title' }, el('b', { text: where ? where[1] : aln.title }),
        where ? el('span', { class: 'muted', text: `${where[2]} ${where[3]}:${Number(where[4]).toLocaleString('en')}–${Number(where[5]).toLocaleString('en')}, ${where[6] === '-' ? 'reverse' : 'forward'} strand, ${plural(names.length - 1, 'allele')}` }) : null,
        identity ? el('span', { class: 'muted', text: `this allele ${(100 * Number(identity[1])).toFixed(1)}% identical` }) : null),
      el('span', { class: 'aln-legend' }, key('', 'supported'), key('nv-inconclusive', 'inconclusive'), key('nv-needs_confirmation', 'needs confirmation'),
        hint('How to read the alignment', help))),
    grid, tip);
}


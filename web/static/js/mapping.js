// How a name of the file is connected to a gene of the list, and the
// control that changes it.

import { el } from './dom.js?v=BUILD';
import { BY_SEQUENCE, LEFT_OUT, decision, geneRows } from './names.js?v=BUILD';
import { sortGenes } from './words.js?v=BUILD';

// How a name of the file is connected to a gene of the list.
export const LINKS = {
  same: { label: 'same name', arrow: 'solid' },
  rule: { label: 'built-in rule', arrow: 'solid', help: 'The gene list knows these alleles under another gene name.' },
  yours: { label: 'your match', arrow: 'solid' },
  sequence: { label: 'by sequence', arrow: 'dashed', help: 'Each allele will be placed on the gene its sequence fits best. The results say which.' },
  out: { label: 'left out', arrow: 'broken', help: 'These alleles are not analysed.' },
  none: { label: 'not matched', arrow: 'broken', help: 'No gene of the list has this name. Choose a gene, let igpop find it by sequence, correct the name, or add the gene to the list.' },
};

/** How an allele is connected to the gene list: a key of LINKS. */
export function linkOf(allele, place) {
  const chosen = decision(allele);
  return chosen === LEFT_OUT ? 'out'
    : allele.gene === BY_SEQUENCE ? 'sequence'
      : !(place && place.build) ? 'none'
        : allele.gene === allele.named ? 'same'
          : chosen === allele.gene ? 'yours' : 'rule';
}

/** The choices for a group or an allele, as a select. The list of all
 *  genes is long, and there can be hundreds of these selects: the genes
 *  are added when the select is first used. */
export function decisionSelect(label, current, similarTo, onChange) {
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


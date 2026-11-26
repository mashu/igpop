// The haplotype check of an opened allele.

import { el } from './dom.js?v=BUILD';
import { groupTitle } from './populations.js?v=BUILD';
import { state } from './state.js?v=BUILD';
import { given, number } from './tsv.js?v=BUILD';
import { changeWords, detailSection } from './variants.js?v=BUILD';
import { HAPLOTYPES } from './words.js?v=BUILD';

// The phased panel is that of 1000 Genomes: its chromosomes are counted
// in the five groups of that project, named as in datasets.tsv.
const PANEL_GROUPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'];
const PANEL = { GRCh38: '1kGP_30x', GRCh37: '1kGP_phase3' };

export function haplotypePanel(allele) {
  const info = HAPLOTYPES[allele.hap_status];
  if (!info || allele.hap_status === 'not_checked') return null;
  const parts = [el('p', { text: info.help })];
  const which = given(allele.hap_exact) ? 'exact' : given(allele.hap_carrying) ? 'carrying' : null;
  // Which differences the phased panel has as sites, by name.
  const absent = new Set([...(allele.hap_untested_str || '').matchAll(/,(\d+)\)/g)].map((m) => m[1]));
  const label = (v) => (given(v.rsid) ? v.rsid : `the ${changeWords(v)} at position ${v.aln_pos}`);
  const mine = state.variants.get(allele.allele_name) || [];
  const [inPanel, notInPanel] = [mine.filter((v) => !absent.has(v.genome_pos)).map(label), mine.filter((v) => absent.has(v.genome_pos)).map(label)];
  const listed = (names) => (names.length > 1 ? `${names.slice(0, -1).join(', ')} and ${names[names.length - 1]}` : names[0] || '');
  if (notInPanel.length) {
    parts.push(el('p', { class: 'muted', text: `The phased panel has ${inPanel.length ? listed(inPanel) : 'none of the differences'} as ${inPanel.length === 1 ? 'a site' : 'sites'}, but not ${listed(notInPanel)}. The panel records only variants found among its own people and kept by its quality filters, so a difference can be seen in other data and still be missing here; how often each one is seen is in the table above.` }));
  }
  if (which) {
    const total = number(allele[`hap_${which}`]);
    const donors = given(allele.hap_exact_donors) ? `, in ${allele.hap_exact_donors} ${allele.hap_exact_donors === '1' ? 'person' : 'people'}` : '';
    // What the bars count: the whole sequence, or only the part that
    // could be looked up.
    const title = which === 'exact' ? 'Phased chromosomes with exactly this sequence'
      : inPanel.length === 1 ? `Phased chromosomes that carry ${inPanel[0]}`
        : `Phased chromosomes that carry ${listed(inPanel)} together`;
    const caption = which === 'exact'
      ? `${total.toLocaleString('en')} of ${Number(allele.hap_total).toLocaleString('en')} chromosomes of unrelated people${donors}.`
      : `${total.toLocaleString('en')} of ${Number(allele.hap_total).toLocaleString('en')} chromosomes of unrelated people.`;
    const shares = PANEL_GROUPS.map((pop) => {
      const [count, size] = (allele[`hap_${which}_${pop}`] || '0/0').split('/').map(Number);
      return { pop, count, size, share: size ? count / size : 0 };
    });
    const top = Math.max(...shares.map((s) => s.share), 1e-9);
    parts.push(el('div', { class: 'hap-chart' },
      el('h5', { text: title }),
      el('p', { text: caption }),
      el('table', { class: 'bars' }, el('caption', { class: 'sr', text: `${title}, per population group` }),
        el('tbody', {}, shares.map((s) => el('tr', {},
          el('th', { scope: 'row', text: groupTitle(allele.build, PANEL[allele.build], s.pop) }),
          el('td', {}, el('span', { class: 'bar-fill', style: `--share:${(100 * s.share) / top}%` })),
          el('td', { class: 'num', text: `${s.count.toLocaleString('en')} of ${s.size.toLocaleString('en')}` }),
          el('td', { class: 'num', text: `${(100 * s.share).toFixed(s.share < 0.01 ? 2 : 1)}%` }))))),
      which === 'carrying' ? el('p', { class: 'muted', text: notInPanel.length
        ? `This counts ${inPanel.length === 1 ? `${inPanel[0]} alone` : 'these differences only'}, not the allele: ${listed(notInPanel)} ${notInPanel.length === 1 ? 'is' : 'are'} not in the panel, so the check cannot say whether a chromosome carries all differences of the allele together.`
        : 'Chromosomes that carry these differences also carry further variants inside the gene, so none has exactly the sequence of this allele.' }) : null));
    const coincidence = number(allele.hap_coincidence);
    if (which === 'exact' && total > 0 && coincidence >= 0) {
      const share = coincidence < 0.0005 ? 'fewer than 1 in 2000' : `${(100 * coincidence).toFixed(coincidence < 0.01 ? 2 : 1)}% of`;
      parts.push(el('p', { class: 'muted', text: `For comparison: ${share} arbitrary single-base changes of a real chromosome are found on this many chromosomes.` }));
    }
    const p = number(allele.hap_exact_pop_heterogeneity_p);
    if (which === 'exact' && p < 0.001) parts.push(el('p', { class: 'muted', text: 'The frequency differs clearly between populations: this haplotype is population-specific.' }));
  }
  return detailSection('Haplotype check', 'Does the whole sequence, with all its differences together, exist on single chromosomes of real people?', ...parts);
}


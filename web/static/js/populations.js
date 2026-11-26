// Who carries a difference: the populations that looked at it, with
// their counts.

import { el } from './dom.js?v=BUILD';
import { state } from './state.js?v=BUILD';
import { given, number, parseTSV } from './tsv.js?v=BUILD';

export const POPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'];

export const POP_NAMES = { AFR: 'African', AMR: 'American', EAS: 'East Asian', EUR: 'European', SAS: 'South Asian' };

/** Every frequency reported for a variant, as {kind, name, af}. */
export function frequencies(v) {
  const out = [];
  for (const pop of POPS) {
    const af = number(v[`kgp_AF_${pop}`]);
    if (af > 0) out.push({ kind: 'kgp', name: `1000 Genomes ${POP_NAMES[pop]}`, af });
  }
  const gnomad = number(v.gnomad_max_group_AF);
  if (gnomad > 0) out.push({ kind: 'gnomad', name: `gnomAD ${v.gnomad_max_group}`, af: gnomad });
  if (given(v.more_groups)) {
    for (const entry of v.more_groups.split(';')) {
      const [source, group, counts] = entry.split(':');
      const [ac, an] = counts.split('/').map(Number);
      if (ac > 0) out.push({ kind: 'dbsnp', name: `${source.replace('_', ' ')} ${group}`, af: ac / an });
    }
  }
  if (given(v.dbsnp_freq)) {
    for (const entry of v.dbsnp_freq.split(';')) {
      const [name, af] = entry.split(':');
      if (Number(af) > 0) out.push({ kind: 'dbsnp', name, af: Number(af) });
    }
  }
  return out;
}

const GNOMAD_NAMES = { afr: 'African / African American', ami: 'Amish', amr: 'Admixed American', asj: 'Ashkenazi Jewish', eas: 'East Asian', fin: 'Finnish', mid: 'Middle Eastern', nfe: 'European (non-Finnish)', sas: 'South Asian' };

export const sourceTitle = (name) => ({ '1kGP_30x': '1000 Genomes, 30x', '1kGP_phase3': '1000 Genomes, phase 3', dbSNP: 'dbSNP', ALFA: 'ALFA by population', GenomeAsia: 'the GenomeAsia pilot', ToMMo_61k: 'ToMMo 61KJPN', WBBC: 'the Westlake BioBank for Chinese' }[name] || name.replace('gnomAD_v', 'gnomAD v'));

/** A stratum such as `1kGP_30x:AFR` in words. */
export function stratumWords(stratum) {
  const [source, population] = stratum.split(':');
  if (source.startsWith('1kGP')) return `${POP_NAMES[population] || population} people of ${sourceTitle(source)}`;
  if (source.startsWith('gnomAD')) return `the ${GNOMAD_NAMES[population] || population} group of ${sourceTitle(source)}`;
  if (MORE_SOURCES[source]) return `${MORE_SOURCES[source].group(population)} of ${MORE_SOURCES[source].title}`;
  return `the study ${population}`;
}

/** The further datasets read beside 1000 Genomes and gnomAD, by the name
 *  the tables use. */
const MORE_SOURCES = {
  ALFA: { title: 'ALFA (NCBI)', group: (g) => `the ${g} group` },
  GenomeAsia: { title: 'the GenomeAsia pilot', group: (g) => `the ${g} group` },
  ToMMo_61k: { title: 'ToMMo 61KJPN', group: (g) => `${g} people` },
  WBBC: { title: 'the Westlake BioBank for Chinese', group: (g) => `people from ${g}` },
};

/** The continental group a population belongs to, which gives its bar a
 *  colour; mixed and unplaced populations have none. */
const ANCESTRIES = [['eur', 'European'], ['afr', 'African'], ['eas', 'East Asian'], ['mid', 'Middle Eastern'], ['sas', 'South Asian'], ['amr', 'American']];

const ANCESTRY_OF = [
  [/^(Other|Other Asian|Oceanian)$/, ''],
  [/African/, 'afr'],
  [/American/, 'amr'],
  [/South Asian/, 'sas'],
  [/East Asian|Northeast Asian|Southeast Asian|Japanese|China|Lingnan|^(KOREAN|Korea\w*|TOMMO|Vietnamese)$/, 'eas'],
  [/Middle Eastern|^Qatari$/, 'mid'],
  [/Europ|Finnish|Amish|Ashkenazi|West Eurasian|^(GoNL|TWINSUK|ALSPAC|Estonian|NorthernSweden|FINRISK|GENOME_DK)$/, 'eur'],
];

export const ancestry = (name) => (ANCESTRY_OF.find(([pattern]) => pattern.test(name)) || [0, ''])[1];

export function ancestryLegend() {
  return el('p', { class: 'legend ancestry-key' },
    ...ANCESTRIES.map(([key, name]) => el('span', { class: 'key' }, el('i', { class: `swatch a-${key}` }), ` ${name}`)),
    el('span', { class: 'key' }, el('i', { class: 'swatch' }), ' mixed or other'));
}

/** Every population that looked at a variant, with its counts:
 *  [{ source, total: {ac, an}, note, rows: [{ name, ac, an, estimated }] }]. */
export function populations(v) {
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
  if (given(v.more_groups)) {
    const bySource = new Map();
    for (const entry of v.more_groups.split(';')) {
      const [source, group, counts] = entry.split(':');
      const [ac, an] = counts.split('/').map(Number);
      if (!bySource.has(source)) bySource.set(source, []);
      bySource.get(source).push({ name: group, ac, an });
    }
    for (const [source, rows] of bySource) {
      const title = MORE_SOURCES[source] ? MORE_SOURCES[source].title.replace(/^the /, '') : source;
      out.push({
        source: title[0].toUpperCase() + title.slice(1),
        total: { ac: rows.reduce((sum, r) => sum + r.ac, 0), an: rows.reduce((sum, r) => sum + r.an, 0) },
        note: '',
        rows,
      });
    }
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


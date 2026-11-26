// Who carries a difference: the populations that looked at it, with
// their counts. Everything here comes from three tables of the run:
// populations.tsv (the counts), site_quality.tsv (how far a dataset's
// counts can be trusted at a position) and datasets.tsv (who the
// populations are). The page knows no dataset by name.

import { el } from './dom.js?v=BUILD';
import { state } from './state.js?v=BUILD';
import { given, number, parseTSV } from './tsv.js?v=BUILD';

// The continental groups, in the order of their key. The names of the
// groups are those of datasets.tsv; a population of none is drawn grey.
const ANCESTRIES = [['european', 'European'], ['african', 'African'], ['east_asian', 'East Asian'], ['middle_eastern', 'Middle Eastern'], ['south_asian', 'South Asian'], ['american', 'American']];

const keyOf = (row) => [row.build, row.contig, row.genome_pos, row.reference_bases, row.query_variant].join('|');

function grouped(text) {
  const out = new Map();
  for (const row of text ? parseTSV(text) : []) {
    const key = keyOf(row);
    if (!out.has(key)) out.set(key, []);
    out.get(key).push(row);
  }
  return out;
}

/** Read the population tables of a run into `state`. */
export function loadPopulations(files) {
  state.populations = grouped(files['populations.tsv']);
  state.quality = grouped(files['site_quality.tsv']);
  // Per build: the datasets in the order of the table, each with its
  // title and its groups by code.
  state.datasets = new Map();
  for (const row of files['datasets.tsv'] ? parseTSV(files['datasets.tsv']) : []) {
    if (!state.datasets.has(row.build)) state.datasets.set(row.build, new Map());
    const datasets = state.datasets.get(row.build);
    if (!datasets.has(row.dataset)) datasets.set(row.dataset, { name: row.dataset, title: row.dataset_title, estimated: row.counts === 'estimated', people: number(row.people), groups: new Map() });
    datasets.get(row.dataset).groups.set(row.group, { title: row.group_title, ancestry: row.ancestry, chromosomes: number(row.chromosomes), used: row.used === 'yes' });
  }
}

/** The datasets of a build, as [{ name, title, estimated, people, groups }]. */
export const datasetsOf = (build) => [...(state.datasets.get(build) || new Map()).values()];

/** A dataset of a build by the name the tables use. */
const dataset = (build, name) => (state.datasets.get(build) || new Map()).get(name);

/** The title of a group of a dataset, such as "African" for AFR. */
export const groupTitle = (build, name, group) => dataset(build, name)?.groups.get(group)?.title || group;

/** A stratum such as `1kGP_30x:AFR`, in words. */
export function stratumWords(build, stratum) {
  const [name, group] = stratum.split(':');
  const from = dataset(build, name);
  if (!from) return stratum;
  return from.estimated ? `the study ${group}, relayed by dbSNP` : `the ${groupTitle(build, name, group)} group of ${from.title}`;
}

/** Every population that looked at a difference, with its counts:
 *  [{ source, total: {ac, an}, note, rows: [{ name, ancestry, ac, an, estimated }] }],
 *  one entry per dataset, in the order of datasets.tsv. */
export function populations(v) {
  const rows = state.populations.get(keyOf(v)) || [];
  const quality = new Map((state.quality.get(keyOf(v)) || []).map((row) => [row.dataset, row]));
  return datasetsOf(v.build).map((from) => {
    const mine = rows.filter((row) => row.dataset === from.name);
    if (!mine.length) return null;
    const site = quality.get(from.name);
    const [one, two] = site ? [number(site.people_one_copy), number(site.people_two_copies)] : [NaN, NaN];
    const sum = (field) => mine.reduce((total, row) => total + Number(row[field]), 0);
    return {
      source: from.title,
      // A dataset gives its own total where its groups do not add up to
      // everyone; studies relayed by dbSNP have no common total.
      total: from.estimated ? null : site ? { ac: Number(site.AC), an: Number(site.AN) } : { ac: sum('AC'), an: sum('AN') },
      note: from.estimated ? 'dbSNP gives frequencies only; the counts are worked out from the estimated size of each study'
        : one > 0 || two > 0 ? [one >= 0 ? `${one} ${one === 1 ? 'person' : 'people'} with one copy` : '', two >= 0 ? `${two} with two` : ''].filter(Boolean).join(', ')
          : '',
      rows: mine.map((row) => {
        const group = from.groups.get(row.group) || { title: row.group, ancestry: 'other' };
        return { name: group.title, ancestry: group.ancestry, ac: Number(row.AC), an: Number(row.AN), estimated: from.estimated };
      }),
    };
  }).filter(Boolean);
}

/** Every population that has seen a difference, as { source, name, ancestry, af }. */
export function frequencies(v) {
  return populations(v).flatMap((from) => from.rows.filter((row) => row.ac > 0).map((row) => ({ source: from.source, name: `${from.source}: ${row.name}`, ancestry: row.ancestry, af: row.ac / row.an })));
}

/** The site criteria a difference's position fails, as [{ source, failed: [keyword] }]. */
export function failures(v) {
  return (state.quality.get(keyOf(v)) || []).filter((row) => given(row.failed) && row.failed !== 'PASS')
    .map((row) => ({ source: dataset(v.build, row.dataset)?.title || row.dataset, failed: row.failed.split(',') }));
}

/** The key to the colours of the population bars. */
export function ancestryLegend() {
  return el('p', { class: 'legend ancestry-key' },
    ...ANCESTRIES.map(([key, name]) => el('span', { class: 'key' }, el('i', { class: `swatch a-${key}` }), ` ${name}`)),
    el('span', { class: 'key' }, el('i', { class: 'swatch a-other' }), ' mixed or other'));
}

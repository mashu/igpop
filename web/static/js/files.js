// Files given to the page: an allele file starts a run; the output files
// of an earlier run are put on display as they are.

import { startCheck } from './check.js?v=BUILD';
import { display } from './results.js?v=BUILD';
import { fail, state } from './state.js?v=BUILD';
import { OUTPUTS } from './words.js?v=BUILD';

export async function handleFiles(fileList) {
  const files = [...fileList];
  if (!files.length) return;
  const outputs = files.filter((f) => OUTPUTS.some((suffix) => f.name.endsWith(`.${suffix}`)));
  if (outputs.length) {
    const loaded = {};
    for (const file of outputs) {
      const suffix = OUTPUTS.find((s) => file.name.endsWith(`.${s}`));
      loaded[suffix] = await file.text();
      state.prefix = file.name.slice(0, -suffix.length - 1);
    }
    if (!loaded['alleles.tsv']) return fail('To browse earlier results, include the .alleles.tsv file (and .variants.tsv and .populations.tsv for the evidence).');
    return display(loaded, []);
  }
  const file = files[0];
  const fasta = await file.text();
  if (!fasta.trimStart().startsWith('>')) return fail(`${file.name} does not look like a FASTA file (it should start with ">").`);
  startCheck(fasta, file.name);
}

/** The example of the start page, read like a file of the user's. */
export async function loadExample() {
  try {
    const response = await fetch(new URL(`../example_KIARVA_IGHV.fasta${new URL(import.meta.url).search}`, import.meta.url));
    if (!response.ok) throw new Error(`status ${response.status}`);
    await handleFiles([new File([await response.text()], 'KIARVA_IGHV.fasta', { type: 'text/plain' })]);
  } catch (error) {
    fail(`The example could not be loaded (${error.message}).`);
  }
}

// igpop web app: the page. The analysis itself runs in a worker
// (worker.js); the modules here start it and turn the tables it returns
// into the views of the page. The tables are exactly the files the
// command-line program writes, so the page can also display those.
//
// The modules, each using only the ones before it:
//   dom, tsv, words      building blocks: elements, table text, fixed words
//   state                what is on display, and a channel for "this changed"
//   engine               the worker: requests and messages
//   populations, score   the data behind a difference: who carries it, its score
//   names, mapping       the gene names of a file and what they are matched to
//   progress             a run in progress
//   variants, haplotype, alignment   the parts of an opened allele
//   limit, pages         the verdict limit; the other tabs of the results
//   results              the results page
//   name-editor, gene-list, check    the page before a run
//   files                files given to the page
//   main                 this file: wiring

import { checkFailed, renderCheck, wireCheck } from './check.js?v=BUILD';
import { $, wireTabs } from './dom.js?v=BUILD';
import { listen, startWorker } from './engine.js?v=BUILD';
import { handleFiles, loadExample } from './files.js?v=BUILD';
import { renderGeneList, wireGeneList } from './gene-list.js?v=BUILD';
import { wireLimit } from './limit.js?v=BUILD';
import { wireNameEditor } from './name-editor.js?v=BUILD';
import { names } from './names.js?v=BUILD';
import { logLine, showStored } from './progress.js?v=BUILD';
import { display, wireResults } from './results.js?v=BUILD';
import { fail, show } from './state.js?v=BUILD';

/** A file is chosen by a click, the keyboard or a drop anywhere. */
function wireFiles() {
  const drop = $('drop');
  const input = $('file');
  drop.addEventListener('click', () => input.click());
  drop.addEventListener('keydown', (event) => {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      input.click();
    }
  });
  input.addEventListener('change', () => handleFiles(input.files).then(() => { input.value = ''; }));
  for (const type of ['dragenter', 'dragover']) {
    document.addEventListener(type, (event) => {
      event.preventDefault();
      drop.classList.add('over');
    });
  }
  for (const type of ['dragleave', 'drop']) document.addEventListener(type, () => drop.classList.remove('over'));
  document.addEventListener('drop', (event) => {
    event.preventDefault();
    if (!$('progress').hidden) return;
    handleFiles(event.dataTransfer.files);
  });
  $('example').addEventListener('click', loadExample);
  for (const button of document.querySelectorAll('[data-restart]')) button.addEventListener('click', () => show('start'));
}

// What each message of the worker leads to.
listen('ready', (message) => {
  if (typeof message.geneList !== 'string') {
    fail('This page and its analysis module are from different versions. Reload the page (Ctrl+Shift+R, or Cmd+Shift+R on a Mac) and try again.');
    return;
  }
  $('version').textContent = `igpop ${message.version}`;
  names.builtin = { coords: message.geneList, aliases: message.aliases };
  renderGeneList();
});
listen('checked', (message) => message.id === names.asked && renderCheck(message.genes));
listen('check-failed', (message) => message.id === names.asked && checkFailed(message.message));
listen('*', (message) => 'stored' in message && showStored(message.stored, message.reading));
listen('log', (message) => logLine(message.message));
listen('error', (message) => fail(message.message));
listen('done', (message) => display(message.output.files, message.output.warnings, message.output.summary));

wireFiles();
wireTabs();
wireResults();
wireCheck();
wireGeneList();
wireLimit();
wireNameEditor();
// The worker is started with the page, so the analysis module is loaded
// and the stored data counted before the first file is dropped.
startWorker();

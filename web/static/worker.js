// Runs the analysis off the main thread. The WebAssembly module first
// downloads the regions it needs, several at a time, and then analyses
// without pauses, reading with synchronous requests, which browsers allow
// only in workers. Downloaded pieces are kept in the
// browser (store.js) so that later runs need little or no downloading.
import init, { check, gene_aliases, gene_list, run, use_store, version } from './pkg/igpop_web.js';
import { Store } from './store.js';

let store = null;
const ready = Promise.all([init(), Store.open()]).then(([, opened]) => {
  store = opened;
  use_store(store);
  postMessage({ type: 'ready', version: version(), stored: store ? store.size : null, geneList: gene_list(), aliases: gene_aliases() });
});

self.onmessage = async (event) => {
  const job = event.data;
  await ready;
  if (job.type === 'clear') {
    if (store) store.clear();
    postMessage({ type: 'stored', stored: store ? store.size : null });
    return;
  }
  if (job.type === 'check') {
    // Which genes the file names and whether the gene list has them:
    // answered at once, nothing is downloaded.
    try {
      postMessage({ type: 'checked', id: job.id, genes: check(job) });
    } catch (error) {
      postMessage({ type: 'check-failed', id: job.id, message: String(error && error.message ? error.message : error) });
    }
    return;
  }
  try {
    const output = await run(job, (message) => postMessage({ type: 'log', message }));
    postMessage({ type: 'done', output, stored: store ? store.size : null });
  } catch (error) {
    postMessage({ type: 'error', message: String(error && error.message ? error.message : error) });
  } finally {
    if (store) store.flush();
  }
};

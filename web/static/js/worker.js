// Runs the analysis off the main thread. The WebAssembly module first
// downloads the regions it needs, several at a time, and then analyses
// without pauses, reading with synchronous requests, which browsers allow
// only in workers. Downloaded pieces are kept in the
// browser (store.js) so that later runs need little or no downloading.
// The page starts this worker with its build stamp (worker.js?v=...). The
// module, its WebAssembly file and the store are loaded with the same
// stamp, so that all parts of the app come from one build.
const stamp = self.location.search;
const at = (path) => new URL(`${path}${stamp}`, import.meta.url);

let store = null;
let igpop = null;
const ready = Promise.all([import(at('../pkg/igpop_web.js')), import(at('./store.js'))]).then(async ([module, { Store }]) => {
  await module.default({ module_or_path: at('../pkg/igpop_web_bg.wasm') });
  igpop = module;
  store = await Store.open();
  igpop.use_store(store);
  postMessage({ type: 'ready', version: igpop.version(), stored: store ? store.size : null, reading: Boolean(store && !store.writer), geneList: igpop.gene_list(), aliases: igpop.gene_aliases() });
});
ready.catch((error) => postMessage({ type: 'error', message: `The analysis module could not be loaded (${error && error.message ? error.message : error}). Reload the page and try again.` }));

self.onmessage = async (event) => {
  const job = event.data;
  await ready;
  if (job.type === 'clear') {
    if (store) store.clear();
    postMessage({ type: 'stored', stored: store ? store.size : null, reading: Boolean(store && !store.writer) });
    return;
  }
  if (job.type === 'check') {
    // Which genes the file names and whether the gene list has them:
    // answered at once, nothing is downloaded.
    try {
      postMessage({ type: 'checked', id: job.id, genes: igpop.check(job) });
    } catch (error) {
      postMessage({ type: 'check-failed', id: job.id, message: String(error && error.message ? error.message : error) });
    }
    return;
  }
  try {
    const output = await igpop.run(job, (message) => postMessage({ type: 'log', message }));
    postMessage({ type: 'done', output, stored: store ? store.size : null, reading: Boolean(store && !store.writer) });
  } catch (error) {
    postMessage({ type: 'error', message: String(error && error.message ? error.message : error) });
  } finally {
    if (store) store.flush();
  }
};

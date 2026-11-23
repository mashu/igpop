// Runs the analysis off the main thread. The WebAssembly module blocks
// until it is done and reads remote data with synchronous requests, which
// browsers allow only in workers.
import init, { run, version } from './pkg/igpop_web.js';

const ready = init().then(() => postMessage({ type: 'ready', version: version() }));

self.onmessage = async (event) => {
  const job = event.data;
  try {
    await ready;
    const output = run(
      job.fasta, job.name, job.build, job.marker, job.haplotypes, job.paralog,
      (message) => postMessage({ type: 'log', message }),
    );
    postMessage({ type: 'done', output });
  } catch (error) {
    postMessage({ type: 'error', message: String(error && error.message ? error.message : error) });
  }
};

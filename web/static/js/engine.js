// The analysis runs in a worker (worker.js). This module starts it, sends
// it requests and hands its messages to whoever asked for them by type.

import { fail } from './state.js?v=BUILD';

let worker = null;
const handlers = {};

/** Call `handler` with every message of this type; '*' = every message. */
export function listen(type, handler) {
  (handlers[type] ||= []).push(handler);
}

/** Send a request to the worker; false if it has not been started. */
export function post(message) {
  if (!worker) return false;
  worker.postMessage(message);
  return true;
}

export function startWorker() {
  // Every file of the app carries the same build stamp in its address, so
  // that a browser never combines a new page with an older, cached worker
  // or analysis module.
  const stamp = new URL(import.meta.url).search;
  worker = new Worker(new URL(`./worker.js${stamp}`, import.meta.url), { type: 'module' });
  worker.onmessage = (event) => {
    for (const type of ['*', event.data.type]) for (const handler of handlers[type] || []) handler(event.data);
  };
  worker.onerror = (event) => fail(`The analysis module could not be started (${event.message || 'unknown error'}). A recent Chrome, Edge, Firefox or Safari is needed.`);
}

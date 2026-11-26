// What is on display, and a channel on which modules say what changed.

import { $ } from './dom.js?v=BUILD';

// What the page holds: the tables of the run on display and the choices
// made about them.
export const state = { limit: 10, prefix: 'igpop', marker: '_S', build: '', files: {}, alleles: [], variants: new Map(), alignments: new Map(), verdict: '', doubt: new Set() };

export function show(section) {
  for (const id of ['start', 'check', 'progress', 'error', 'results']) $(id).hidden = id !== section;
  // The steps of a run, with the present one marked.
  const order = ['start', 'check', 'progress', 'results'];
  $('stages').hidden = !order.includes(section) || section === 'start';
  for (const item of $('stages').children) {
    const at = order.indexOf(item.dataset.step) - order.indexOf(section);
    item.className = at < 0 ? 'done' : at === 0 ? 'now' : '';
    if (at === 0) item.setAttribute('aria-current', 'step'); else item.removeAttribute('aria-current');
  }
  emit('section', section);
}

export function fail(message) {
  $('error-text').textContent = message;
  show('error');
  $('error').querySelector('button').focus();
}

// Candidates are named by the marker, as in the command-line program.
export const isCandidate = (allele) => allele.allele_name.includes(state.marker);

// A module that changes something others show says so here, instead of
// calling them: 'section' (the part of the page on display changed),
// 'verdicts' (the verdicts were worked out again with another limit).
const listeners = {};
export function on(event, handler) {
  (listeners[event] ||= []).push(handler);
}
export function emit(event, ...args) {
  for (const handler of listeners[event] || []) handler(...args);
}

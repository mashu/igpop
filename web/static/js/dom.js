// Small tools for building the page: elements, hints, tabs.

export const $ = (id) => document.getElementById(id);

export const el = (tag, attrs = {}, ...children) => {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value === false || value == null) continue;
    if (key === 'class') node.className = value;
    else if (key === 'text') node.textContent = value;
    else if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? '' : value);
  }
  for (const child of children.flat()) {
    if (child != null) node.append(child);
  }
  return node;
};

export const svg = (tag, attrs = {}, ...children) => {
  const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
  node.append(...children);
  return node;
};

export const plural = (n, word) => `${n} ${word}${n === 1 ? '' : 's'}`;

/** Shows `content` next to `trigger` on hover, keyboard focus or tap.
 *  The explanation is placed on the page, not inside the table, so that a
 *  scrolling table cannot cut it off. */
export function explain(trigger, content) {
  explain.count = (explain.count || 0) + 1;
  const tip = el('span', { class: 'hint-text', role: 'tooltip', id: `hint-${explain.count}` }, content);
  trigger.setAttribute('aria-describedby', tip.id);
  const show = () => {
    tip.classList.add('open');
    const at = trigger.getBoundingClientRect();
    const width = tip.offsetWidth;
    const left = Math.min(Math.max(8, at.left + at.width / 2 - width / 2), window.innerWidth - width - 8);
    const below = at.bottom + 8 + tip.offsetHeight < window.innerHeight;
    tip.style.left = `${left}px`;
    tip.style.top = below ? `${at.bottom + 8}px` : `${Math.max(8, at.top - 8 - tip.offsetHeight)}px`;
  };
  const hide = () => tip.classList.remove('open');
  for (const type of ['mouseenter', 'focus']) trigger.addEventListener(type, show);
  for (const type of ['mouseleave', 'blur']) trigger.addEventListener(type, hide);
  trigger.addEventListener('keydown', (event) => { if (event.key === 'Escape') hide(); });
  return el('span', { class: 'hint' }, trigger, tip);
}

/** A small "?" with an explanation. */
export function hint(label, content) {
  return explain(el('button', { type: 'button', class: 'hint-button', 'aria-label': label, text: '?' }), content);
}

export function selectTab(id) {
  for (const tab of $(id).closest('[role=tablist]').querySelectorAll('[role=tab]')) {
    const selected = tab.id === id;
    tab.setAttribute('aria-selected', String(selected));
    tab.tabIndex = selected ? 0 : -1;
    $(tab.getAttribute('aria-controls')).hidden = !selected;
  }
}

/** Tabs answer to a click and to the arrow keys. */
export function wireTabs() {
  for (const list of document.querySelectorAll('[role=tablist]')) {
    const tabs = [...list.querySelectorAll('[role=tab]')];
    tabs.forEach((tab, i) => {
      tab.addEventListener('click', () => selectTab(tab.id));
      tab.addEventListener('keydown', (event) => {
        const step = { ArrowRight: 1, ArrowLeft: -1 }[event.key];
        if (!step) return;
        const next = tabs[(i + step + tabs.length) % tabs.length];
        selectTab(next.id);
        next.focus();
      });
    });
  }
}

// Downloaded pieces of the public data files, kept in the browser between
// visits (Origin Private File System). Used from the worker only: its
// synchronous file handles are what lets the analysis, which runs without
// pauses, read and write the store while it works.
//
// Layout: `pieces.bin` holds the pieces one after another; `pieces.idx` has
// one line per piece: offset, length, time stored, key.

const MAX_BYTES = 1.5e9;
const MAX_AGE_MS = 180 * 24 * 3600 * 1000;

export class Store {
  /** Opens the store; null where the browser does not offer one, or
   *  another tab of this page holds it. */
  static async open() {
    try {
      const root = await navigator.storage.getDirectory();
      const dir = await root.getDirectoryHandle('igpop', { create: true });
      const open = async (name) => (await dir.getFileHandle(name, { create: true })).createSyncAccessHandle();
      const store = new Store(await open('pieces.bin'), await open('pieces.idx'));
      store.load();
      return store;
    } catch (error) {
      return null;
    }
  }

  constructor(data, index) {
    this.data = data;
    this.index = index;
    this.pieces = new Map();
    this.dataEnd = 0;
    this.indexEnd = 0;
    this.unflushed = 0;
  }

  load() {
    const text = new Uint8Array(this.index.getSize());
    this.index.read(text, { at: 0 });
    let oldest = Date.now();
    for (const line of new TextDecoder().decode(text).split('\n')) {
      const [offset, length, time, key] = line.split('\t');
      if (key === undefined) continue;
      this.pieces.set(decodeURIComponent(key), [Number(offset), Number(length)]);
      oldest = Math.min(oldest, Number(time));
    }
    this.dataEnd = this.data.getSize();
    this.indexEnd = text.length;
    // The store only grows; when it is large or old it starts afresh.
    if (this.dataEnd > MAX_BYTES || Date.now() - oldest > MAX_AGE_MS) this.clear();
  }

  /** The piece stored under `key`, or undefined. */
  get(key) {
    const piece = this.pieces.get(key);
    if (!piece) return undefined;
    const bytes = new Uint8Array(piece[1]);
    return this.data.read(bytes, { at: piece[0] }) === piece[1] ? bytes : undefined;
  }

  put(key, bytes) {
    if (this.pieces.has(key)) return;
    this.data.write(bytes, { at: this.dataEnd });
    const line = new TextEncoder().encode(`${this.dataEnd}\t${bytes.length}\t${Date.now()}\t${encodeURIComponent(key)}\n`);
    this.index.write(line, { at: this.indexEnd });
    this.pieces.set(key, [this.dataEnd, bytes.length]);
    this.dataEnd += bytes.length;
    this.indexEnd += line.length;
    this.unflushed += 1;
    if (this.unflushed >= 50) this.flush();
  }

  flush() {
    this.data.flush();
    this.index.flush();
    this.unflushed = 0;
  }

  clear() {
    this.data.truncate(0);
    this.index.truncate(0);
    this.flush();
    this.pieces.clear();
    this.dataEnd = 0;
    this.indexEnd = 0;
  }

  /** Bytes stored. */
  get size() {
    return this.dataEnd;
  }
}

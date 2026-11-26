// Downloaded pieces of the public data files, kept in the browser between
// visits (Origin Private File System). Used from the worker only: its
// synchronous file handles are what lets the analysis, which runs without
// pauses, read and write the store while it works.
//
// Several tabs can use the store at once: one adds to it, the others read
// it (see `Store.open`).
//
// Layout: `pieces.bin` holds the pieces one after another; `pieces.idx` has
// one line per piece: offset, length, time stored, key.

const MAX_BYTES = 1.5e9;
const MAX_AGE_MS = 180 * 24 * 3600 * 1000;

const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

/** True if this worker may add to the store: only one tab at a time does,
 *  and it stays the one until it is closed or reloaded. */
function claimWriting() {
  if (!navigator.locks) return Promise.resolve(true);
  return new Promise((resolve) => {
    navigator.locks.request('igpop-store-writer', { ifAvailable: true }, (lock) => {
      resolve(Boolean(lock));
      // Held for the life of this worker.
      return lock ? new Promise(() => {}) : undefined;
    }).catch(() => resolve(true));
  });
}

export class Store {
  /** Opens the store; null where the browser does not offer one.
   *
   *  Several tabs can have it open: one of them adds to it (`writer`), the
   *  others read what is there, see what the writer adds, and keep their
   *  own downloads in memory. Where the browser only allows one tab to
   *  open the files, a tab that finds them taken tries again for a few
   *  seconds, which is what a reload needs while the old page lets go. */
  static async open() {
    try {
      const root = await navigator.storage.getDirectory();
      const dir = await root.getDirectoryHandle('igpop', { create: true });
      const files = [await dir.getFileHandle('pieces.bin', { create: true }), await dir.getFileHandle('pieces.idx', { create: true })];
      const both = async (options) => {
        const data = await files[0].createSyncAccessHandle(options);
        try {
          return [data, await files[1].createSyncAccessHandle(options)];
        } catch (error) {
          data.close();
          throw error;
        }
      };
      let handles = null;
      let shared = true;
      for (let attempt = 0; attempt < 12 && !handles; attempt += 1) {
        try {
          handles = await both({ mode: 'readwrite-unsafe' });
        } catch {
          // No shared mode in this browser, or another tab holds the files
          // alone: open them alone, waiting a little for a page that is
          // being reloaded to let go.
          try {
            handles = await both();
            shared = false;
          } catch {
            await wait(300);
          }
        }
      }
      if (!handles) return null;
      const store = new Store(handles[0], handles[1]);
      store.writer = shared ? await claimWriting() : true;
      store.load();
      return store;
    } catch {
      return null;
    }
  }

  constructor(data, index) {
    this.data = data;
    this.index = index;
    this.pieces = new Map();
    // Pieces of a tab that reads only: kept until the tab is closed.
    this.own = new Map();
    this.writer = true;
    this.dataEnd = 0;
    this.indexEnd = 0;
    this.unflushed = 0;
  }

  /** Read the lines of the index from `indexEnd` on. */
  read() {
    const size = this.index.getSize();
    if (size <= this.indexEnd) return Date.now();
    const bytes = new Uint8Array(size - this.indexEnd);
    this.index.read(bytes, { at: this.indexEnd });
    // Only whole lines: the writer may be in the middle of one.
    const whole = bytes.lastIndexOf(10) + 1;
    let oldest = Date.now();
    for (const line of new TextDecoder().decode(bytes.subarray(0, whole)).split('\n')) {
      const [offset, length, time, key] = line.split('\t');
      if (key === undefined) continue;
      this.pieces.set(decodeURIComponent(key), [Number(offset), Number(length)]);
      oldest = Math.min(oldest, Number(time));
    }
    this.indexEnd += whole;
    return oldest;
  }

  load() {
    const oldest = this.read();
    this.dataEnd = this.data.getSize();
    // The store only grows; when it is large or old it starts afresh.
    if (this.writer && (this.dataEnd > MAX_BYTES || Date.now() - oldest > MAX_AGE_MS)) this.clear();
  }

  /** The piece stored under `key`, or undefined. */
  get(key) {
    const mine = this.own.get(key);
    if (mine) return mine;
    let piece = this.pieces.get(key);
    if (!piece && !this.writer) {
      // Another tab may have added it since.
      this.read();
      piece = this.pieces.get(key);
    }
    if (!piece) return undefined;
    const bytes = new Uint8Array(piece[1]);
    return this.data.read(bytes, { at: piece[0] }) === piece[1] ? bytes : undefined;
  }

  put(key, bytes) {
    if (this.pieces.has(key) || this.own.has(key)) return;
    if (!this.writer) {
      this.own.set(key, bytes.slice());
      return;
    }
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

  /** Empties the store. A tab that only reads leaves it alone. */
  clear() {
    this.own.clear();
    if (!this.writer) return;
    this.data.truncate(0);
    this.index.truncate(0);
    this.flush();
    this.pieces.clear();
    this.dataEnd = 0;
    this.indexEnd = 0;
  }

  /** Bytes stored. */
  get size() {
    return this.writer ? this.dataEnd : this.data.getSize();
  }
}

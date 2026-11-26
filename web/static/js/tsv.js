// The tables of a run as text: reading them into rows.

export function parseTSV(text) {
  const lines = text.split('\n').filter((line) => line.length);
  const header = lines.shift().split('\t');
  return lines.map((line) => {
    const cells = line.split('\t');
    return Object.fromEntries(header.map((name, i) => [name, cells[i] ?? '.']));
  });
}

export const given = (cell) => cell != null && cell !== '.' && cell !== '-' && cell !== '';

export const number = (cell) => (given(cell) ? Number(cell) : NaN);

/** Blocks of alignments.txt by gene name. */
export function parseAlignments(text) {
  const blocks = new Map();
  for (const block of (text || '').split(/\n(?=>)/)) {
    const match = block.match(/^>(\S+)/);
    if (match) blocks.set(match[1], block.trimEnd());
  }
  return blocks;
}

/** One gene block of alignments.txt as columns: the positions on the
 *  genome, and one string per row (reference first). */
export function parseBlock(block) {
  const lines = block.split('\n');
  const where = lines[0].match(/^>(\S+)\s+(\S+)\s+(\S+):(\d+)-(\d+)/) || [];
  const out = { title: lines[0].slice(1), build: where[2] || '', contig: where[3] || '', positions: [], rows: new Map(), notes: new Map() };
  let [width, first] = [0, 0];
  for (const line of lines.slice(1)) {
    const chunk = line.match(/^\s+(\d+) \.\. (\d+)\s*$/);
    if (chunk) {
      const [from, to] = [Number(chunk[1]), Number(chunk[2])];
      width = Math.abs(to - from) + 1;
      for (let k = 0; k < width; k += 1) out.positions.push(from + (to >= from ? k : -k));
      continue;
    }
    const note = line.match(/^ {2}(\S+)\s+(identity .*)$/);
    if (note) {
      out.notes.set(note[1], note[2]);
      continue;
    }
    const row = line.match(/^(\S+)(\s+)/);
    if (!row || !width || line.startsWith(' ')) continue;
    // The reference row comes first and has no blanks: it shows in which
    // column of the text the sequences begin.
    if (row[1] === 'reference') first = row[0].length;
    const before = out.positions.length - width;
    const seq = (out.rows.get(row[1]) || '').padEnd(before, ' ') + line.slice(first).padEnd(width, ' ').slice(0, width);
    out.rows.set(row[1], seq);
  }
  return out;
}


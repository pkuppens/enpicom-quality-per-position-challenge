# Reflection and Lessons Learned

## Original Steps (Initial Attempt)

1. Explore the assignment
2. Explore the repo — created uv venv with Python 3.12.11
3. Start on the assignment — parser first, support plain and gzip
4. Read per 4 lines (LINES_PER_READ)
5. Extract ID, quality string — convert to bytes for histogram efficiency
6. Abandoned histogram; rewrote with sorted approach

---

## Core Mistake: Per Sequence vs Per Position

**The implementation computed one boxplot per sequence (read), not per position.**

- The spec says: *"one boxplot data structure per read position, summarizing the distribution of quality scores across all reads at that position"*
- The implementation used `sequence_id` as `Boxplot.name` and produced 250 boxplots (one per record)
- The correct result: one boxplot per position index 0..max_position, so 493 boxplots for the fixture (reads vary from 250 to 493 bases)

**Lesson:** Read the spec carefully; "per position" means transpose scores across reads by position index, not one statistic set per read.

---

## What Worked

- **Streaming:** Reading FASTQ in chunks of 4 lines avoids loading the whole file into memory.
- **Sanger encoding:** `ord(char) - 33` for Phred quality.
- **Percentile interpolation:** Virtual index with linear interpolation between adjacent sorted values; straightforward and correct.

---

## What Did Not Work

- **Histogram-based O(n):** Attempted for efficiency; complexity outweighed benefit for typical sizes; reverted to sorted approach.
- **Per-sequence boxplot:** Misread the requirement; should have been per-position from the start.
- **`Boxplot.name`:** Used `sequence_id`; spec requires `str(position_index)`.

---

## Variable-Length Reads

The fixture has reads of length 250–493. The transpose step must:

- Group scores by position index across all reads
- Only include a read at position `i` if that read has length > `i`
- Use an extensible structure (`defaultdict(list)`) — later reads can be longer than earlier ones

---

## Defensive Programming

- Validate quality string length vs sequence length; raise `ValueError` on mismatch
- Validate quality characters (ASCII ≥ 33 for Sanger)
- Handle empty file → return `[]`
- Never close stdin when reading from `"-"`; use `contextlib.nullcontext`

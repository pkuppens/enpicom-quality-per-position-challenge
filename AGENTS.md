# Agent Skills and Learnings

Project-level guidance for AI agents working on this codebase. May be migrated to personal skills later.

---

## Quality Per Position Challenge

### What Works

- **Streaming FASTQ reads** in chunks of 4 lines (`LINES_PER_READ`); avoid loading the entire file.
- **Sanger encoding:** `ord(char) - 33` for Phred quality scores.
- **Transpose by position:** `defaultdict(list)` mapping `position -> [scores]` from all reads; extensible for variable-length reads.
- **Percentile interpolation:** Virtual index with linear interpolation between adjacent sorted values; always return float.
- **stdio support:** Pass `"-"` as filename; use `contextlib.nullcontext(sys.stdin)` so stdin is never closed.

### What Failed (Initial Attempt)

- Computing boxplot **per sequence** instead of **per position**.
- Using `sequence_id` as `Boxplot.name`; spec requires `str(position_index)`.

### Conventions

- README mentions `.fastq`; fixtures use `.fq`. Both are supported; `open_fastq_stream` checks `.gz` suffix for gzip.
- `Boxplot.name` must be `str(0-based position)`.
- Support stdin via `"-"`; never close stdio.
- Extensible position handling: later reads may be longer than earlier reads; do not preallocate.

### File Extensions (ADR 0001)

- Both `.fq` and `.fastq` (and gzipped variants) are valid.
- `Boxplot` TypedDict includes `name: str`; upstream stub omitted it.

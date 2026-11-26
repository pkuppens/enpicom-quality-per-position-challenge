# Quality Per Position — Algorithm Specification (Pseudo-code)

This document defines the canonical algorithm for `quality_per_position_boxplot_data`. **Python, Rust, and Mojo implementations must be identical** in structure and logic. All parsing, decoding, and transposition functions must be transcribed to each language.

**Native execution:** Rust and Mojo run as standalone binaries (not Python extensions), so each can use language-specific optimizations (SIMD, zero-cost abstractions) without FFI overhead.

## Constants

```
LINES_PER_READ = 4
SEQUENCE_ID = 0
SEQUENCE = 1
OPTIONAL_SEQUENCE_ID = 2
QUALITY_STRING = 3
SANGER_ASCII_QUALITY_SCORE_OFFSET = 33
```

## Data Types

```
CenterLine = { min: float, max: float }
Box = { q1: float, q3: float }
Boxplot = { centerLine: CenterLine, box: Box, median: float, average: float, name: str }
```

**Standalone binaries (Rust, Mojo):** Read FASTQ path from argv[1]; write a JSON array of Boxplot to stdout. Example:

```json
[{"centerLine":{"min":28.0,"max":42.0},"box":{"q1":35.0,"q3":40.0},"median":38.0,"average":37.2,"name":"0"},...]
```

---

## Function 1: open_fastq_stream(filename: str) -> Stream

Open a FASTQ input stream. Never close stdin.

```
IF filename == "-":
    RETURN stdin (no-close context)
IF filename ends with ".gz":
    RETURN gzip.open(filename, "rt", encoding="utf-8")
ELSE:
    RETURN open(filename, "r", encoding="utf-8")
```

---

## Function 2: read_fastq_records(stream: Stream) -> Iterator[Record]

Yield FASTQ records as lists of LINES_PER_READ lines. Each record = [id_line, sequence_line, optional_id_line, quality_line].

```
LOOP:
    record = [stream.readline() for i in 1..LINES_PER_READ]
    IF len(record) < LINES_PER_READ OR record[0] is empty:
        RETURN (stop iteration)
    YIELD record
```

---

## Function 3: parse_quality_scores(quality_string: str) -> List[int]

Decode quality string to integer scores using Sanger/Phred+33 encoding.

**Input:** quality_string — raw quality line (after strip)  
**Output:** list of integer Phred scores, one per character

```
RETURN [ord(c) - SANGER_ASCII_QUALITY_SCORE_OFFSET for c in quality_string]
```

Assumes valid input: all chars have ASCII >= 33. Caller validates.

---

## Function 4: transpose_by_position(records: Iterator[Record]) -> Map[int, List[int]]

Transpose quality scores from per-read to per-position. Key = 0-based position index; value = list of scores from all reads at that position.

```
position_to_scores = empty map (int -> list of int)

FOR each record IN records:
    sequence = record[SEQUENCE].strip()
    quality_string = record[QUALITY_STRING].strip()
    scores = parse_quality_scores(quality_string)
    FOR i, score IN enumerate(scores):
        position_to_scores[i].append(score)

RETURN position_to_scores
```

---

## Function 5: compute_position_stats(scores_at_position: List[int], position: int) -> Boxplot

Compute boxplot statistics for a single position using histogram-based percentiles.

**Input:** scores_at_position — quality scores at this position; position — 0-based index  
**Output:** Boxplot with centerLine, box, median, average, name  
**Raises:** if scores_at_position is empty

### 5.1 Build histogram

```
n = len(scores_at_position)
IF n == 0: RAISE

hist = empty map (score -> count)
FOR s IN scores_at_position:
    hist[s] += 1

min_score = min(hist.keys())
max_score = max(hist.keys())
average = sum(scores_at_position) / n
```

### 5.2 Cumulative counts (sorted order)

```
scores_sorted = sorted(hist.keys())
cum_before = empty map (score -> int)
cum = 0
FOR s IN scores_sorted:
    cum_before[s] = cum
    cum += hist[s]
```

### 5.3 Score at sorted index k (0-based)

```
FUNCTION score_at_sorted_index(k: int) -> float:
    k = clamp(k, 0, n - 1)
    FOR s IN scores_sorted (in order):
        IF cum_before[s] <= k < cum_before[s] + hist[s]:
            RETURN float(s)
    RETURN float(max_score)
```

### 5.4 Percentile interpolation

Virtual indices: vi = p * (n - 1) for p in [0.25, 0.5, 0.75]. Linear interpolation between adjacent sorted elements.

```
FUNCTION percentile(p: float) -> float:
    vi = p * (n - 1)
    lower_idx = floor(vi)
    upper_idx = lower_idx + 1

    IF upper_idx >= n:
        RETURN score_at_sorted_index(lower_idx)
    ELSE:
        score_lower = score_at_sorted_index(lower_idx)
        score_upper = score_at_sorted_index(upper_idx)
        frac = vi - lower_idx
        RETURN (1 - frac) * score_lower + frac * score_upper

q1 = percentile(0.25)
median = percentile(0.5)
q3 = percentile(0.75)
```

### 5.5 Return Boxplot

```
RETURN Boxplot(
    centerLine = { min: min_score, max: max_score },
    box = { q1: q1, q3: q3 },
    median = median,
    average = average,
    name = str(position)
)
```

---

## Function 6: quality_per_position_boxplot_data(fastq_filename: str) -> List[Boxplot]

Main entry point. Orchestrates parsing, transposition, and per-position stats.

```
position_to_scores = empty map (int -> list of int)

WITH open_fastq_stream(fastq_filename) AS stream:
    FOR record IN read_fastq_records(stream):
        sequence = record[SEQUENCE].strip()
        quality_string = record[QUALITY_STRING].strip()
        scores = parse_quality_scores(quality_string)
        FOR i, score IN enumerate(scores):
            position_to_scores[i].append(score)

IF position_to_scores is empty:
    RETURN []

max_position = max(position_to_scores.keys())   # 0-based, so positions 0..max_position

RETURN [
    compute_position_stats(position_to_scores[i], i)
    FOR i IN 0..max_position (inclusive)
]
```

---

## Function Mapping (Python ↔ Rust ↔ Mojo)

| Pseudo-code | Python | Rust | Mojo |
|-------------|--------|------|------|
| open_fastq_stream | `open_fastq_stream` | `open_fastq_stream` | `open_fastq_stream` |
| read_fastq_records | `read_fastq_records` | `read_fastq_records` | `read_fastq_records` |
| parse_quality_scores | `_parse_quality_scores` | `parse_quality_scores` | `parse_quality_scores` |
| transpose_by_position | inline in main | `transpose_by_position` or inline | `transpose_by_position` or inline |
| compute_position_stats | `_compute_position_stats` | `compute_position_stats` | `compute_position_stats` |
| quality_per_position_boxplot_data | `quality_per_position_boxplot_data` | `quality_per_position_boxplot_data` | `quality_per_position_boxplot_data` |

All three implementations must expose the same functions and produce bit-identical (or float-tolerance-identical) results for the same input.

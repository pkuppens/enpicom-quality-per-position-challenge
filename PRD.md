# Product Requirements Document: Quality Per Position

## Reference

- **Original problem:** [README.md](README.md) — FASTQ format and `quality_per_position_boxplot_data` specification
- **Challenge:** [CHALLENGE.md](CHALLENGE.md) — Implementation expectations
- **Upstream:** [ENPICOM/quality-per-position-challenge](https://github.com/ENPICOM/quality-per-position-challenge)

---

## Functional Requirements

| ID | Requirement | Traceability |
|----|-------------|--------------|
| FR-1 | Parse FASTQ format (plain text or gzip-compressed) | README: "FASTQ files may be plain text (.fastq) or gzip-compressed (.fastq.gz)" |
| FR-2 | Read records as exactly four lines (ID, sequence, separator, quality) | README: FASTQ format description |
| FR-3 | Convert quality string to integer scores using Sanger/Phred+33 encoding | README: "subtracting 33 from the ASCII value gives the quality score" |
| FR-4 | Transpose quality scores by position index across all reads | README step 2: "groups scores from all reads by their position index" |
| FR-5 | Handle variable-length reads; position i includes only reads with length > i | README: "Reads can differ in length" |
| FR-6 | Compute per-position statistics: min, max, q1, median, q3, average | README step 3: table of fields |
| FR-7 | Use linear interpolation between adjacent sorted values for percentiles | README: "Percentiles are computed with linear interpolation" |
| FR-8 | Return `list[Boxplot]` with one entry per position, ordered 0..max_position | README step 4 |
| FR-9 | Each `Boxplot.name` is the string representation of the 0-based position index | README: "stored in the `name` field" |

---

## Internal Requirements (Ambiguity Resolution)

| ID | Requirement | Rationale |
|----|-------------|-----------|
| IR-1 | Support both `.fq` and `.fastq` (and `.fq.gz`, `.fastq.gz`) | README mentions `.fastq`; fixtures use `.fq`; both are common |
| IR-2 | Support stdin via `"-"` as filename; never close stdio | Enables piping; closing stdin would break shell pipelines |
| IR-3 | Extensible position handling: later reads may be longer than earlier reads | Do not preallocate; extend position map as new positions appear |
| IR-4 | Raise `ValueError` on malformed records (e.g. quality length ≠ sequence length) | No silent skip; fail fast for bad data |
| IR-5 | Raise `ValueError` on invalid quality characters (e.g. ASCII < 33) | Defensive input validation |
| IR-6 | Return `[]` for empty file or empty stdin | Defensible empty result |

---

## Non-Functional Requirements / Best Practices

| ID | Requirement | Rationale |
|----|-------------|-----------|
| NF-1 | Single responsibility per module/function | SOLID, maintainability |
| NF-2 | Defensive programming; validate inputs before processing | Robustness |
| NF-3 | Stream FASTQ records; avoid loading entire file into memory | Scalability for large files |
| NF-4 | Use `defaultdict(list)` or equivalent for position→scores; never preallocate fixed size | Extensible for variable-length reads |
| NF-5 | Percentile computation via virtual index with linear interpolation; return float for clients | Consistency with spec; interoperability |

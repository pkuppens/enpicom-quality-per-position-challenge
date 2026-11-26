# ADR 0001: File extensions and TypedDict

## Status

Accepted

## Context

The README specifies `.fastq` and `.fastq.gz`, while the upstream fixtures and current code use `test.fq` and `test.fq.gz`. The upstream `quality_per_position.py` stub omits the `name` field from the `Boxplot` TypedDict; the README explicitly requires it.

## Decision

| Topic | README | Implementation | Decision |
|-------|--------|----------------|----------|
| Extensions | `.fastq`, `.fastq.gz` | Fixtures: `test.fq`, `test.fq.gz` | Both `.fq` and `.fastq` are valid; `open_fastq_file` checks `.gz` only. Document that both are supported; fixtures stay `.fq`. |
| Boxplot `name` | Required, `str(position_index)` | Upstream stub omitted `name` | Add `name: str` to TypedDict; implementation sets `name=str(i)`. |
| Test docstring | - | `test_open_gzipped_fastq` says ".fastq.gz" but file is `test.fq.gz` | Fix docstring to ".fq.gz" or "gzipped FASTQ" for clarity. |

## Consequences

- `open_fastq_file` detects gzip by `.gz` suffix; file extension (.fq vs .fastq) does not affect behavior.
- `Boxplot` TypedDict includes `name: str` as the 0-based position index string.
- Test docstrings are updated for consistency with actual fixture names.

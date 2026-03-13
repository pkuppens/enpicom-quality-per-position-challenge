# [FEAT]: Python vs Rust vs Mojo Benchmark for Quality Per Position

## Goal

Demonstrate speed gains of the `quality_per_position_boxplot_data` algorithm by implementing Mojo and Rust ports, with reproducible benchmarks, synthetic data, and a three-way comparison table (Python vs Rust vs Mojo).

## Background

The reference implementation in [`quality_per_position.py`](../quality_per_position.py) reads FASTQ files, transposes quality scores by position, and computes per-position boxplot statistics (min, max, q1, median, q3, average) using histogram-based percentiles. The existing fixture `fixtures/test.fq` has ~250 reads and ~493 positions—too small to show meaningful performance differences. Mojo and Rust implementations will enable a fair comparison and document setup for users new to Mojo.

**Implementation constraint:** Python, Rust, and Mojo implementations must be **identical in pseudo-code**. All parsing, decoding, and transposition functions must be transcribed to each language. See [docs/quality-per-position-algorithm-spec.md](docs/quality-per-position-algorithm-spec.md) for the canonical specification.

**Native execution:** Rust and Mojo run as **standalone native binaries**, not Python extensions. This allows each to use language-specific optimizations (SIMD, zero-cost abstractions, etc.) without FFI overhead. The benchmark harness invokes them via subprocess.

## Mojo Environment Setup (for Users New to Mojo)

**Note:** Mojo has no Windows wheels. On Windows, use WSL2 or Docker (see [docs/benchmark-mojo-rust-setup.md](docs/benchmark-mojo-rust-setup.md)).

### Windows: WSL2 (recommended)

1. Install WSL2 + Ubuntu 24.04: `wsl --install -d Ubuntu-24.04`
2. Open project in Cursor via WSL path: `\\wsl$\Ubuntu-24.04\home\<user>\...\enpicom-quality-per-position-challenge`
3. In WSL terminal: create Linux venv (do not reuse Windows `.venv`), then:
   ```bash
   uv venv && source .venv/bin/activate
   uv pip install mojo --extra-index-url https://modular.gateway.scarf.sh/simple/
   mojo --version
   ```

### Linux / macOS: Add Mojo (uv)

| Step | Action |
|------|--------|
| 1 | Activate venv: `source .venv/bin/activate` |
| 2 | Install Mojo: `uv pip install mojo --extra-index-url https://modular.gateway.scarf.sh/simple/` |
| 3 | Verify: `mojo --version` |

### Linux / macOS: Add Mojo (Pixi)

| Step | Action |
|------|--------|
| 1 | `pixi init -c https://conda.modular.com/max-nightly/ -c conda-forge` |
| 2 | `pixi add mojo` |
| 3 | `pixi shell` → `mojo --version` |

### IDE

Install "Mojo" extension from Modular (VS Code / Cursor).

### System requirements

From [Modular install docs](https://docs.modular.com/mojo/manual/install):

- **macOS**: Apple Silicon (M1–M5), Xcode 16+, Python 3.10–3.14
- **Linux**: Ubuntu 24.04, x86-64 (SSE4.2) or Graviton, g++/clang++, Python 3.10–3.14
- **Windows**: WSL2 with Ubuntu 24.04

### Quick verification

```bash
# After pixi shell or uv venv activation
echo 'def main() raises: print("Hello, Mojo!")' > hello.mojo
mojo hello.mojo
# Expected: Hello, Mojo!
```

## Tasks

- [ ] Add `tmp/` to `.gitignore`
- [ ] Create `scripts/generate_benchmark_fastq.py` (synthetic data generator)
- [ ] Create `scripts/benchmark_quality_per_position.py` (benchmark harness)
- [ ] Implement `quality_per_position.mojo` — standalone binary; transcribe all spec functions; use Mojo-native optimizations; read path from argv, write JSON to stdout
- [ ] Implement Rust crate `quality_per_position` — standalone binary; transcribe all spec functions; use Rust-native optimizations; read path from argv, write JSON to stdout
- [ ] Add tests: `test_benchmark_data_generator.py`, `test_mojo_quality_per_position.py`, `test_rust_quality_per_position.py`
- [ ] Create `docs/benchmark-mojo-rust-setup.md` with full Mojo setup (extracted from this issue)
- [ ] Run benchmarks and document comparison table in README or docs

## Synthetic Benchmark Data

The generator creates FASTQ files on first run, stored in `tmp/` (gitignored). Data is never committed.

| Tier | Reads | Read length | Total bases | Purpose |
|------|-------|-------------|-------------|---------|
| S | 1,000 | 250 | 250k | Sanity check |
| M | 50,000 | 250 | 12.5M | Typical workload |
| L | 500,000 | 300 | 150M | Stress test |
| XL | 2,000,000 | 350 | 700M | Where Python struggles |

**Generator script** (`scripts/generate_benchmark_fastq.py`):

- Accept tier (S/M/L/XL) or custom `--reads N --length L`
- Output to `tmp/benchmark_{tier}.fq`
- Use Sanger encoding; quality chars ASCII 33–126 (Phred 0–93)
- Deterministic (fixed seed) for reproducibility
- Create `tmp/` if missing

## Implementation Plan

All implementations must follow [docs/quality-per-position-algorithm-spec.md](docs/quality-per-position-algorithm-spec.md). Transcribe these functions identically:

| Function | Purpose |
|----------|---------|
| `open_fastq_stream` | Open file, gzip, or stdin |
| `read_fastq_records` | Parse 4-line FASTQ records |
| `parse_quality_scores` | Decode Sanger/Phred+33 quality string → `List[int]` |
| (inline transposition) | Map position → scores across reads |
| `compute_position_stats` | Histogram-based Q1/median/Q3, min, max, average |
| `quality_per_position_boxplot_data` | Main entry: read → transpose → compute |

### Mojo

- **Standalone binary** (no Python FFI): `mojo run quality_per_position.mojo <path>` or `mojo build` → `./quality_per_position <path>`
- Transcribe all functions from the spec; use Mojo-native optimizations (SIMD, etc.)
- Reads FASTQ path from argv; writes boxplot JSON to stdout
- Invoked by benchmark harness via subprocess

### Rust

- **Standalone binary** (no Python FFI): `cargo build --release` → `./target/release/quality_per_position <path>`
- Transcribe all functions from the spec; use Rust-native optimizations (zero-cost abstractions, etc.)
- Reads FASTQ path from argv; writes boxplot JSON to stdout
- Invoked by benchmark harness via subprocess

### Benchmark harness

- Script: `scripts/benchmark_quality_per_position.py`
- **Python**: in-process, direct import
- **Rust**: subprocess `./target/release/quality_per_position <path>`
- **Mojo**: subprocess `mojo run quality_per_position.mojo <path>` (or built binary)
- Flow: (1) Ensure synthetic data exists; (2) Run each impl N times, record median; (3) Output comparison table

## Comparison Table (target output format)

| Tier | Reads | Python (s) | Rust (s) | Mojo (s) | Rust speedup | Mojo speedup |
|------|-------|------------|----------|----------|--------------|--------------|
| S | 1k | — | — | — | — | — |
| M | 50k | — | — | — | — | — |
| L | 500k | — | — | — | — | — |
| XL | 2M | — | — | — | — | — |

*(Values to be filled by benchmark runs; illustrative expectations: Rust 5–9x, Mojo 10–19x over Python on larger tiers.)*

## Acceptance Criteria

- [ ] **Generator**: Run `python scripts/generate_benchmark_fastq.py S`; `tmp/benchmark_S.fq` exists, has 1000 reads × 250 bases; `python -c "from quality_per_position import quality_per_position_boxplot_data; r=quality_per_position_boxplot_data('tmp/benchmark_S.fq'); assert len(r)==250 and r[0]['name']=='0'"` exits 0
- [ ] **Benchmark harness**: Run `python scripts/benchmark_quality_per_position.py S`; prints table with Python timing; exits 0
- [ ] **Mojo correctness**: `mojo run quality_per_position.mojo fixtures/test.fq` outputs JSON; parse and compare to Python result (structure match; float tolerance). Mojo must implement all spec functions as native Mojo.
- [ ] **Mojo uses GPU acceleration**: Mojo supports it; and it is supposed to speed up significantly
- [ ] **Rust correctness**: `./target/release/quality_per_position fixtures/test.fq` outputs JSON; parse and compare to Python result. Rust must implement all spec functions as native Rust.
- [ ] **Regression**: `uv run pytest tests/` exits 0; existing tests unchanged
- [ ] **Docs**: `docs/benchmark-mojo-rust-setup.md` exists with full Mojo setup instructions

## Tests

- `tests/test_benchmark_data_generator.py` — generator produces valid FASTQ, correct structure, deterministic output
- `tests/test_mojo_quality_per_position.py` — run Mojo binary via subprocess, parse JSON stdout, assert structure and values match Python (skip if Mojo not installed)
- `tests/test_rust_quality_per_position.py` — run Rust binary via subprocess, parse JSON stdout, assert structure and values match Python (skip if binary not built)

## Out of Scope

- Changing the Python reference implementation
- Production deployment of Mojo/Rust versions
- Supporting Windows native (Mojo requires WSL on Windows)
- Benchmarking gzip vs plain FASTQ separately (use plain for simplicity)

## Estimate

**Size: L** — Multiple language implementations, benchmark harness, synthetic data generator, tests, and documentation.

## Metadata

- Labels: `enhancement`, `benchmark`, `size:L`
- Milestone: (none)

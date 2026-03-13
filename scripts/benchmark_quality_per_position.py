"""Benchmark Python vs Rust vs Mojo implementations of quality_per_position_boxplot_data.

Ensures synthetic data exists, runs each implementation N times, and prints
a comparison table. Rust and Mojo run as native binaries (subprocess); Python
runs in-process. Rust and Mojo are optional; skipped if not available.

Usage:
    python scripts/benchmark_quality_per_position.py [S|M|L|XL]
    python scripts/benchmark_quality_per_position.py --path tmp/benchmark_S.fq
"""

import argparse
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

# Add project root for imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

DEFAULT_RUNS = 5
TIERS = ("S", "M", "L", "XL")

# Paths for native binaries
RUST_BINARY = PROJECT_ROOT / "target" / "release" / "quality_per_position"
MOJO_SCRIPT = PROJECT_ROOT / "quality_per_position.mojo"


def ensure_benchmark_data(tier: str) -> Path:
    """Generate benchmark FASTQ if it does not exist."""
    out_path = PROJECT_ROOT / "tmp" / f"benchmark_{tier}.fq"
    if not out_path.exists():
        subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "generate_benchmark_fastq.py"), tier],
            cwd=PROJECT_ROOT,
            check=True,
        )
    return out_path


def benchmark_python(path: Path, runs: int) -> float:
    """Run Python implementation in-process; return median time in seconds."""
    times: list[float] = []
    for _ in range(runs):
        t0 = time.perf_counter()
        from quality_per_position import quality_per_position_boxplot_data

        quality_per_position_boxplot_data(str(path))
        times.append(time.perf_counter() - t0)
    return statistics.median(times)  # 'averages' over runs


def benchmark_mojo(path: Path, runs: int) -> float | None:
    """Run Mojo implementation via subprocess; return median time or None."""
    if not MOJO_SCRIPT.exists() or not shutil.which("mojo"):
        return None

    times: list[float] = []
    for _ in range(runs):
        t0 = time.perf_counter()
        subprocess.run(
            ["mojo", "run", str(MOJO_SCRIPT), str(path)],
            cwd=PROJECT_ROOT,
            capture_output=True,
            check=True,
        )
        times.append(time.perf_counter() - t0)
    return statistics.median(times)


def benchmark_rust(path: Path, runs: int) -> float | None:
    """Run Rust implementation via subprocess; return median time or None."""
    if not RUST_BINARY.exists():
        return None

    times: list[float] = []
    for _ in range(runs):
        t0 = time.perf_counter()
        subprocess.run(
            [str(RUST_BINARY), str(path)],
            cwd=PROJECT_ROOT,
            capture_output=True,
            check=True,
        )
        times.append(time.perf_counter() - t0)
    return statistics.median(times)


def print_table(
    tier: str,
    reads: int | str,
    py_time: float,
    rust_time: float | None,
    mojo_time: float | None,
) -> None:
    """Print comparison table row."""
    rust_speedup = f"{py_time / rust_time:.1f}x" if rust_time and rust_time > 0 else "—"
    mojo_speedup = f"{py_time / mojo_time:.1f}x" if mojo_time and mojo_time > 0 else "—"
    rust_s = f"{rust_time:.3f}" if rust_time is not None else "—"
    mojo_s = f"{mojo_time:.3f}" if mojo_time is not None else "—"
    reads_str = str(reads)
    print(f"| {tier:3} | {reads_str:>6} | {py_time:.3f} | {rust_s:>6} | {mojo_s:>6} | {rust_speedup:>10} | {mojo_speedup:>10} |")


def main() -> None:
    """Run benchmarks and print comparison table."""
    parser = argparse.ArgumentParser(
        description="Benchmark quality_per_position implementations."
    )
    parser.add_argument(
        "tier",
        nargs="?",
        default="S",
        choices=TIERS,
        help="Tier: S, M, L, XL (default: S)",
    )
    parser.add_argument(
        "--path",
        type=Path,
        help="Use specific FASTQ file instead of tier",
    )
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS, help="Runs per impl")
    args = parser.parse_args()

    if args.path:
        path = Path(args.path).resolve()
        if not path.exists():
            print(f"Error: {path} not found", file=sys.stderr)
            sys.exit(1)
        tier = "custom"
        reads = "?"
    else:
        path = ensure_benchmark_data(args.tier)
        tier = args.tier
        reads = {"S": 1_000, "M": 50_000, "L": 500_000, "XL": 2_000_000}[tier]

    print(f"Benchmarking {path} ({tier}, {reads} reads) — {args.runs} runs each\n")

    py_time = benchmark_python(path, args.runs)
    rust_time = benchmark_rust(path, args.runs)
    mojo_time = benchmark_mojo(path, args.runs)

    print("| Tier | Reads | Python (s) | Rust (s) | Mojo (s) | Rust speedup | Mojo speedup |")
    print("|------|-------|------------|----------|----------|--------------|--------------|")
    print_table(tier, reads, py_time, rust_time, mojo_time)


if __name__ == "__main__":
    main()

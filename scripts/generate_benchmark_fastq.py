"""Generate synthetic FASTQ files for benchmark comparisons.

Creates deterministic FASTQ files in tmp/ for benchmarking Python vs Rust vs Mojo
implementations of quality_per_position_boxplot_data. Data is never committed.

Usage:
    python scripts/generate_benchmark_fastq.py S          # tier S (1k reads, 250)
    python scripts/generate_benchmark_fastq.py M           # tier M (50k reads)
    python scripts/generate_benchmark_fastq.py --reads 10000 --length 200  # custom
"""

import argparse
import random
import sys
from pathlib import Path

# Sanger/Phred+33: ASCII 33 = score 0, ASCII 126 = score 93
SANGER_ASCII_OFFSET = 33
QUALITY_RANGE = (0, 93)  # Phred scores
BASES = "ACGTN"

TIERS = {
    "S": {"reads": 1_000, "length": 250},
    "M": {"reads": 50_000, "length": 250},
    "L": {"reads": 500_000, "length": 300},
    "XL": {"reads": 2_000_000, "length": 350},
}

DEFAULT_SEED = 42


def _score_to_char(score: int) -> str:
    """Convert Phred score to Sanger ASCII character."""
    return chr(score + SANGER_ASCII_OFFSET)


def _generate_random_quality(length: int, rng: random.Random) -> str:
    """Generate quality string with scores in typical range (28–42)."""
    # Typical Illumina scores: 28–42 (ASCII ; to I)
    low, high = 28, 42
    return "".join(_score_to_char(rng.randint(low, high)) for _ in range(length))


def _generate_sequence(length: int, rng: random.Random) -> str:
    """Generate random nucleotide sequence."""
    return "".join(rng.choices(BASES, k=length))


def generate_fastq(
    output_path: Path,
    num_reads: int,
    read_length: int,
    seed: int = DEFAULT_SEED,
) -> None:
    """Write synthetic FASTQ file to output_path.

    Args:
        output_path: Path to write .fq file.
        num_reads: Number of reads to generate.
        read_length: Length of each read (bases).
        seed: Random seed for reproducibility.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)

    with open(output_path, "w", encoding="utf-8") as f:
        for i in range(num_reads):
            seq = _generate_sequence(read_length, rng)
            qual = _generate_random_quality(read_length, rng)
            f.write(f"@read_{i}/1\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")


def main() -> None:
    """Parse args and generate FASTQ."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic FASTQ for benchmark comparisons."
    )
    parser.add_argument(
        "tier",
        nargs="?",
        choices=list(TIERS),
        help="Tier: S, M, L, XL",
    )
    parser.add_argument("--reads", type=int, help="Custom number of reads")
    parser.add_argument("--length", type=int, help="Custom read length")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: tmp/benchmark_{tier}.fq)",
    )
    args = parser.parse_args()

    if args.tier:
        tier = TIERS[args.tier]
        num_reads = args.reads or tier["reads"]
        read_length = args.length or tier["length"]
        out_name = args.output or Path("tmp") / f"benchmark_{args.tier}.fq"
    elif args.reads and args.length:
        num_reads = args.reads
        read_length = args.length
        out_name = args.output or Path("tmp") / "benchmark_custom.fq"
    else:
        parser.error("Provide tier (S/M/L/XL) or both --reads and --length")
        sys.exit(1)

    out_path = Path(out_name)
    generate_fastq(out_path, num_reads, read_length, seed=args.seed)
    print(f"Generated {out_path} ({num_reads} reads × {read_length} bases)")


if __name__ == "__main__":
    main()

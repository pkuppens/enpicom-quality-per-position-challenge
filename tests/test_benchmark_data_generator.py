"""Tests for the synthetic FASTQ benchmark data generator."""

import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GENERATOR = PROJECT_ROOT / "scripts" / "generate_benchmark_fastq.py"
TMP = PROJECT_ROOT / "tmp"


def test_generator_produces_valid_fastq(tmp_path):
    """Generator produces valid FASTQ with correct structure."""
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "S", "--output", str(tmp_path / "bench_S.fq")],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    out = tmp_path / "bench_S.fq"
    assert out.exists()

    lines = out.read_text().strip().split("\n")
    assert len(lines) >= 4  # At least one read
    assert lines[0].startswith("@")
    assert lines[1]  # sequence
    assert lines[2] == "+"
    assert len(lines[3]) == len(lines[1])  # quality matches sequence length


def test_generator_tier_s_has_expected_reads_and_length(tmp_path):
    """Tier S produces 1000 reads of length 250."""
    out = tmp_path / "bench_S.fq"
    subprocess.run(
        [sys.executable, str(GENERATOR), "S", "--output", str(out)],
        cwd=PROJECT_ROOT,
        capture_output=True,
        check=True,
    )
    lines = out.read_text().strip().split("\n")
    # Count records (4 lines each); avoid counting @ in quality string (Phred 31 = '@')
    num_reads = len(lines) // 4
    assert num_reads == 1000
    # First read sequence length
    first_seq = lines[1]
    assert len(first_seq) == 250


def test_generator_deterministic(tmp_path):
    """Same seed produces identical output."""
    out1 = tmp_path / "a.fq"
    out2 = tmp_path / "b.fq"
    for out in (out1, out2):
        subprocess.run(
            [sys.executable, str(GENERATOR), "S", "--output", str(out), "--seed", "42"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            check=True,
        )
    assert out1.read_text() == out2.read_text()


def test_generated_fastq_parsable_by_quality_per_position(tmp_path):
    """Generated file is parsable and yields correct boxplot structure."""
    out = tmp_path / "bench_S.fq"
    subprocess.run(
        [sys.executable, str(GENERATOR), "S", "--output", str(out)],
        cwd=PROJECT_ROOT,
        capture_output=True,
        check=True,
    )

    from quality_per_position import quality_per_position_boxplot_data

    result = quality_per_position_boxplot_data(str(out))
    assert len(result) == 250
    assert result[0]["name"] == "0"
    assert "centerLine" in result[0]
    assert "box" in result[0]
    assert "median" in result[0]
    assert "average" in result[0]

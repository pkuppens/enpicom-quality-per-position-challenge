"""Tests for Rust implementation of quality_per_position_boxplot_data.

Runs Rust as a native binary via subprocess. Skips if quality_per_position binary
is not built. Rust binary must output JSON (list of Boxplot) to stdout.
"""

import json
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUST_BINARY = PROJECT_ROOT / "target" / "release" / "quality_per_position"
FASTQ_PLAIN = "fixtures/test.fq"


def _rust_available() -> bool:
    """Check if Rust binary exists."""
    return RUST_BINARY.exists()


def _run_rust(fastq_path: str) -> list:
    """Run Rust binary, capture JSON stdout, parse to list of boxplots."""
    result = subprocess.run(
        [str(RUST_BINARY), str(PROJECT_ROOT / fastq_path)],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


@pytest.mark.skipif(not _rust_available(), reason="Rust quality_per_position binary not built")
def test_rust_result_structure_matches_python():
    """Rust result has same structure as Python (list of boxplots)."""
    from quality_per_position import quality_per_position_boxplot_data

    py_result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    rust_result = _run_rust(FASTQ_PLAIN)

    assert isinstance(rust_result, list)
    assert len(rust_result) == len(py_result)


@pytest.mark.skipif(not _rust_available(), reason="Rust quality_per_position binary not built")
def test_rust_result_values_match_python():
    """Rust boxplot values match Python within float tolerance."""
    from quality_per_position import quality_per_position_boxplot_data

    py_result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    rust_result = _run_rust(FASTQ_PLAIN)

    assert len(rust_result) == len(py_result)
    for i, (py_bp, rust_bp) in enumerate(zip(py_result, rust_result)):
        assert rust_bp["name"] == py_bp["name"] == str(i)
        assert abs(rust_bp["centerLine"]["min"] - py_bp["centerLine"]["min"]) < 0.01
        assert abs(rust_bp["centerLine"]["max"] - py_bp["centerLine"]["max"]) < 0.01
        assert abs(rust_bp["box"]["q1"] - py_bp["box"]["q1"]) < 0.01
        assert abs(rust_bp["box"]["q3"] - py_bp["box"]["q3"]) < 0.01
        assert abs(rust_bp["median"] - py_bp["median"]) < 0.01
        assert abs(rust_bp["average"] - py_bp["average"]) < 0.01

"""Tests for Mojo implementation of quality_per_position_boxplot_data.

Runs Mojo as a native binary via subprocess. Skips if Mojo or quality_per_position.mojo
is not available. Mojo binary must output JSON (list of Boxplot) to stdout.
"""

import json
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MOJO_SCRIPT = PROJECT_ROOT / "quality_per_position.mojo"
FASTQ_PLAIN = "fixtures/test.fq"


def _mojo_available() -> bool:
    """Check if Mojo binary/script exists and mojo CLI is installed."""
    import shutil

    return MOJO_SCRIPT.exists() and shutil.which("mojo") is not None


def _run_mojo(fastq_path: str) -> list:
    """Run Mojo binary, capture JSON stdout, parse to list of boxplots."""
    result = subprocess.run(
        ["mojo", "run", str(MOJO_SCRIPT), str(PROJECT_ROOT / fastq_path)],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


@pytest.mark.skipif(not _mojo_available(), reason="Mojo or quality_per_position.mojo not available")
def test_mojo_result_structure_matches_python():
    """Mojo result has same structure as Python (list of boxplots)."""
    from quality_per_position import quality_per_position_boxplot_data

    py_result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    mojo_result = _run_mojo(FASTQ_PLAIN)

    assert isinstance(mojo_result, list)
    assert len(mojo_result) == len(py_result)


@pytest.mark.skipif(not _mojo_available(), reason="Mojo or quality_per_position.mojo not available")
def test_mojo_result_values_match_python():
    """Mojo boxplot values match Python within float tolerance."""
    from quality_per_position import quality_per_position_boxplot_data

    py_result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    mojo_result = _run_mojo(FASTQ_PLAIN)

    assert len(mojo_result) == len(py_result)
    for i, (py_bp, mojo_bp) in enumerate(zip(py_result, mojo_result)):
        assert mojo_bp["name"] == py_bp["name"] == str(i)
        assert abs(mojo_bp["centerLine"]["min"] - py_bp["centerLine"]["min"]) < 0.01
        assert abs(mojo_bp["centerLine"]["max"] - py_bp["centerLine"]["max"]) < 0.01
        assert abs(mojo_bp["box"]["q1"] - py_bp["box"]["q1"]) < 0.01
        assert abs(mojo_bp["box"]["q3"] - py_bp["box"]["q3"]) < 0.01
        assert abs(mojo_bp["median"] - py_bp["median"]) < 0.01
        assert abs(mojo_bp["average"] - py_bp["average"]) < 0.01

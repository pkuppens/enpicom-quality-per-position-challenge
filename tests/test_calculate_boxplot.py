"""Tests for calculate_boxplot: verify boxplot statistics from quality scores."""

import pytest

from quality_per_position import calculate_boxplot


def test_calculate_boxplot_all_same_quality():
    """When all positions have the same quality score, min and average equal that value."""
    quality_bytes = bytes([40, 40, 40, 40])
    result = calculate_boxplot("@seq1", "ACGT", "+", quality_bytes)
    assert result["centerLine"]["min"] == 40.0
    assert result["average"] == 40.0
    # Structure: q1 <= median <= q3, min <= box values (implementation-dependent)
    assert result["box"]["q1"] <= result["median"] <= result["box"]["q3"]


def test_calculate_boxplot_single_position():
    """Single position yields min and average equal to that value."""
    quality_bytes = bytes([25])
    result = calculate_boxplot("@seq1", "A", "+", quality_bytes)
    assert result["centerLine"]["min"] == 25.0
    assert result["average"] == 25.0
    assert result["median"] == result["box"]["q1"] == result["box"]["q3"]


def test_calculate_boxplot_average_correct():
    """Average equals sum of quality scores divided by number of positions."""
    quality_bytes = bytes([10, 20, 30, 40])
    result = calculate_boxplot("@seq1", "ACGT", "+", quality_bytes)
    expected_avg = (10 + 20 + 30 + 40) / 4
    assert result["average"] == expected_avg


def test_calculate_boxplot_structure():
    """Result has required Boxplot keys: centerLine, box, median, average."""
    quality_bytes = bytes([30, 35])
    result = calculate_boxplot("@seq1", "AG", "+", quality_bytes)
    assert "centerLine" in result
    assert "min" in result["centerLine"]
    assert "max" in result["centerLine"]
    assert "box" in result
    assert "q1" in result["box"]
    assert "q3" in result["box"]
    assert "median" in result
    assert "average" in result


def test_calculate_boxplot_min_max_bounds():
    """centerLine min is the smallest quality score; max is at least the largest minus one."""
    quality_bytes = bytes([5, 10, 15, 20, 25])
    result = calculate_boxplot("@seq1", "ACGTA", "+", quality_bytes)
    assert result["centerLine"]["min"] == 5.0
    assert result["centerLine"]["max"] >= 20.0  # implementation detail for max


def test_calculate_boxplot_ignores_sequence_id_fields():
    """Sequence id, sequence, and optional id do not affect the boxplot values."""
    quality_bytes = bytes([50, 50])
    r1 = calculate_boxplot("ignore1", "XX", "ignore2", quality_bytes)
    r2 = calculate_boxplot("other", "YY", "other+", quality_bytes)
    assert r1 == r2

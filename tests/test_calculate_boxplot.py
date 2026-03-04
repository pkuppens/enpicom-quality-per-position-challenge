"""Tests for calculate_boxplot: verify boxplot statistics from quality scores."""

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
    assert result["box"]["q1"] == 10.0
    assert result["median"] == 15.0
    assert result["box"]["q3"] == 20.0
    assert result["centerLine"]["max"] == 25.0
    assert result["average"] == 15.0


def test_calculate_boxplot_linear_interpolation():
    """Check linear interpolation between adjacent sorted values."""
    quality_bytes = bytes([0, 100])
    result = calculate_boxplot("@linear_interpolation", "ACGTA", "+", quality_bytes)
    assert result["centerLine"]["min"] == 0.0
    assert result["box"]["q1"] == 25.0
    assert result["median"] == 50.0
    assert result["box"]["q3"] == 75.0
    assert result["centerLine"]["max"] == 100.0
    assert result["average"] == 50.0

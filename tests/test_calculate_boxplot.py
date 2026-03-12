"""Tests for _compute_position_stats: verify boxplot statistics from scores at a position."""

from quality_per_position import _compute_position_stats


def test_compute_position_stats_all_same_quality():
    """When all scores are equal, min and average equal that value."""
    scores = [40.0, 40.0, 40.0, 40.0]
    result = _compute_position_stats(scores, 0)
    assert result["centerLine"]["min"] == 40.0
    assert result["average"] == 40.0
    assert result["name"] == "0"
    assert result["box"]["q1"] <= result["median"] <= result["box"]["q3"]


def test_compute_position_stats_single_score():
    """Single score yields min and average equal to that value."""
    scores = [25.0]
    result = _compute_position_stats(scores, 5)
    assert result["centerLine"]["min"] == 25.0
    assert result["average"] == 25.0
    assert result["name"] == "5"
    assert result["median"] == result["box"]["q1"] == result["box"]["q3"]


def test_compute_position_stats_average_correct():
    """Average equals sum of scores divided by count."""
    scores = [10.0, 20.0, 30.0, 40.0]
    result = _compute_position_stats(scores, 0)
    expected_avg = (10 + 20 + 30 + 40) / 4
    assert result["average"] == expected_avg


def test_compute_position_stats_structure():
    """Result has required Boxplot keys: centerLine, box, median, average, name."""
    scores = [30.0, 35.0]
    result = _compute_position_stats(scores, 1)
    assert "centerLine" in result
    assert "min" in result["centerLine"]
    assert "max" in result["centerLine"]
    assert "box" in result
    assert "q1" in result["box"]
    assert "q3" in result["box"]
    assert "median" in result
    assert "average" in result
    assert "name" in result
    assert result["name"] == "1"


def test_compute_position_stats_min_max_bounds():
    """centerLine min is smallest; max is largest; percentiles in order."""
    scores = [5.0, 10.0, 15.0, 20.0, 25.0]
    result = _compute_position_stats(scores, 0)
    assert result["centerLine"]["min"] == 5.0
    assert result["box"]["q1"] == 10.0
    assert result["median"] == 15.0
    assert result["box"]["q3"] == 20.0
    assert result["centerLine"]["max"] == 25.0
    assert result["average"] == 15.0


def test_compute_position_stats_linear_interpolation():
    """Check linear interpolation between adjacent sorted values."""
    scores = [0.0, 100.0]
    result = _compute_position_stats(scores, 0)
    assert result["centerLine"]["min"] == 0.0
    assert result["box"]["q1"] == 25.0
    assert result["median"] == 50.0
    assert result["box"]["q3"] == 75.0
    assert result["centerLine"]["max"] == 100.0
    assert result["average"] == 50.0

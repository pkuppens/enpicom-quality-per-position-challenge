"""Tests for quality_per_position_boxplot_data: end-to-end integration."""

from quality_per_position import quality_per_position_boxplot_data

# Fixture path relative to project root
FASTQ_PLAIN = "fixtures/test.fq"
FASTQ_GZ = "fixtures/test.fq.gz"


def test_quality_per_position_boxplot_data_returns_list():
    """Result is a list of boxplot dictionaries."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    assert isinstance(result, list)


def test_quality_per_position_boxplot_data_record_count():
    """Number of boxplots equals number of FASTQ records in the file."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    # test.fq has 1000 lines = 250 records (4 lines per read)
    assert len(result) == 250


def test_quality_per_position_boxplot_data_each_item_has_required_keys():
    """Each boxplot has centerLine, box, median, and average."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    for bp in result:
        assert "centerLine" in bp
        assert "box" in bp
        assert "median" in bp
        assert "average" in bp
        assert "min" in bp["centerLine"]
        assert "max" in bp["centerLine"]
        assert "q1" in bp["box"]
        assert "q3" in bp["box"]


def test_quality_per_position_boxplot_data_numeric_types():
    """All boxplot values are numeric (int or float)."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    for bp in result:
        assert isinstance(bp["centerLine"]["min"], (int, float))
        assert isinstance(bp["centerLine"]["max"], (int, float))
        assert isinstance(bp["box"]["q1"], (int, float))
        assert isinstance(bp["box"]["q3"], (int, float))
        assert isinstance(bp["median"], (int, float))
        assert isinstance(bp["average"], (int, float))


def test_quality_per_position_boxplot_data_plain_and_gz_match():
    """Plain .fq and gzipped .fq.gz produce identical boxplot results."""
    plain_result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    gz_result = quality_per_position_boxplot_data(FASTQ_GZ)
    assert plain_result == gz_result

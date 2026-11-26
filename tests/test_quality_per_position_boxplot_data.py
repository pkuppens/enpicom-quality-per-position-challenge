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
    """Number of boxplots equals max read length (one per position)."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    # test.fq has reads of length 250-493; max position is 492, so 493 boxplots
    assert len(result) == 493


def test_quality_per_position_boxplot_data_each_item_has_required_keys():
    """Each boxplot has centerLine, box, median, average, and name."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    for bp in result:
        assert "centerLine" in bp
        assert "box" in bp
        assert "median" in bp
        assert "average" in bp
        assert "name" in bp
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


def test_quality_per_position_boxplot_data_name_is_position():
    """Each boxplot name is the string representation of its 0-based position index."""
    result = quality_per_position_boxplot_data(FASTQ_PLAIN)
    for i, bp in enumerate(result):
        assert bp["name"] == str(i)


def test_empty_file_returns_empty_list():
    """Empty FASTQ file yields empty list of boxplots."""
    import io
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fq", delete=False) as f:
        f.write("")
        empty_path = f.name
    try:
        result = quality_per_position_boxplot_data(empty_path)
        assert result == []
    finally:
        import os
        os.unlink(empty_path)


def test_single_read_single_position():
    """Single read with one position yields one boxplot with name '0'."""
    import io
    content = "@id\nA\n+\nI\n"  # I = ASCII 73, score 40
    with io.StringIO(content) as s:
        # We need a file path - use temp file
        import tempfile
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fq", delete=False) as f:
            f.write(content)
            path = f.name
    try:
        result = quality_per_position_boxplot_data(path)
        assert len(result) == 1
        assert result[0]["name"] == "0"
        assert result[0]["centerLine"]["min"] == 40.0
        assert result[0]["average"] == 40.0
    finally:
        import os
        os.unlink(path)


def test_stdin_via_dash(monkeypatch):
    """Reading from '-' uses stdin; verify result from mocked stdin."""
    import io
    content = "@id\nAC\n+\nII\n"  # 2 positions, score 40 each
    monkeypatch.setattr("sys.stdin", io.StringIO(content))
    result = quality_per_position_boxplot_data("-")
    assert len(result) == 2
    assert result[0]["name"] == "0"
    assert result[1]["name"] == "1"


def test_quality_length_mismatch_raises():
    """Malformed record (quality length != sequence length) raises ValueError."""
    import tempfile
    import os
    content = "@id\nAC\n+\nI\n"  # sequence len 2, quality len 1
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fq", delete=False) as f:
        f.write(content)
        path = f.name
    try:
        import pytest
        with pytest.raises(ValueError, match="does not match"):
            quality_per_position_boxplot_data(path)
    finally:
        os.unlink(path)


def test_invalid_quality_char_raises():
    """Invalid quality character (ASCII < 33) raises ValueError."""
    import tempfile
    import os
    content = "@id\nAAA\n+\nA A\n"  # space (32) in middle; gives score -1
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fq", delete=False) as f:
        f.write(content)
        path = f.name
    try:
        import pytest
        with pytest.raises(ValueError, match="Invalid quality"):
            quality_per_position_boxplot_data(path)
    finally:
        os.unlink(path)


def test_later_read_longer_than_earlier():
    """Reads with increasing length extend the position map correctly."""
    import tempfile
    import os
    # Read 1: 2 bases. Read 2: 3 bases. Read 3: 4 bases.
    content = "@r1\nAG\n+\nII\n@r2\nACG\n+\nIII\n@r3\nACGT\n+\nIIII\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fq", delete=False) as f:
        f.write(content)
        path = f.name
    try:
        result = quality_per_position_boxplot_data(path)
        assert len(result) == 4  # positions 0,1,2,3
        assert result[0]["name"] == "0"
        assert result[3]["name"] == "3"
        # Position 3 has only 1 score (from read 3)
        assert result[3]["centerLine"]["min"] == 40.0
        assert result[3]["average"] == 40.0
    finally:
        os.unlink(path)

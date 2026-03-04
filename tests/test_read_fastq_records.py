"""Tests for read_fastq_records: verify FASTQ records are yielded as 4-line chunks."""

import io

import pytest

from quality_per_position import LINES_PER_READ, read_fastq_records


def test_read_fastq_records_empty_stream():
    """Empty stream yields no records."""
    stream = io.StringIO("")
    records = list(read_fastq_records(stream))
    assert records == []


def test_read_fastq_records_incomplete_record_yielded_with_empty_lines():
    """Stream ending mid-record yields one record; missing lines are empty strings."""
    stream = io.StringIO("@id1\nACGT\n")
    records = list(read_fastq_records(stream))
    assert len(records) == 1
    assert records[0][0].strip() == "@id1"
    assert records[0][1].strip() == "ACGT"
    assert records[0][2] == ""
    assert records[0][3] == ""


def test_read_fastq_records_single_record():
    """Stream with one complete record yields one 4-line record."""
    stream = io.StringIO("@seq1\nACGT\n+\n!!!!\n")
    records = list(read_fastq_records(stream))
    assert len(records) == 1
    assert len(records[0]) == LINES_PER_READ
    assert records[0][0].strip() == "@seq1"
    assert records[0][1].strip() == "ACGT"
    assert records[0][2].strip() == "+"
    assert records[0][3].strip() == "!!!!"


def test_read_fastq_records_multiple_records():
    """Stream with multiple complete records yields all records."""
    content = "@seq1\nACGT\n+\n!!!!\n@seq2\nTGCA\n+\nIIII\n"
    stream = io.StringIO(content)
    records = list(read_fastq_records(stream))
    assert len(records) == 2
    assert records[0][1].strip() == "ACGT"
    assert records[1][1].strip() == "TGCA"


def test_read_fastq_records_trailing_newline_preserved():
    """Lines include trailing newline when present in source."""
    stream = io.StringIO("@seq1\nACGT\n+\n!!!!\n")
    records = list(read_fastq_records(stream))
    assert records[0][0] == "@seq1\n"
    assert records[0][1] == "ACGT\n"

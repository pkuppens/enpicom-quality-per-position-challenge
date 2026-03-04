"""Tests for open_fastq_file: verify plain and gzipped FASTQ files open and yield identical content."""

import pytest

from quality_per_position import open_fastq_file

# Fixture paths relative to project root
FASTQ_PLAIN = "fixtures/test.fq"
FASTQ_GZ = "fixtures/test.fq.gz"


def read_first_four_lines(path: str) -> list[str]:
    """Read the first four lines from a FASTQ file path."""
    with open_fastq_file(path) as f:
        return [f.readline() for _ in range(4)]


def test_open_plain_fastq_reads_four_lines():
    """Plain .fq file opens and yields at least four readable lines."""
    lines = read_first_four_lines(FASTQ_PLAIN)
    assert len(lines) == 4
    assert lines[0].startswith("@")
    assert lines[1]  # sequence
    assert lines[2].startswith("+")
    assert len(lines[3].strip()) == len(lines[1].strip())  # quality same length as sequence


def test_open_gzipped_fastq_reads_four_lines():
    """Gzipped .fastq.gz file opens and yields at least four readable lines."""
    lines = read_first_four_lines(FASTQ_GZ)
    assert len(lines) == 4
    assert lines[0].startswith("@")
    assert lines[1]  # sequence
    assert lines[2].startswith("+")
    assert len(lines[3].strip()) == len(lines[1].strip())


def test_plain_and_gzipped_content_correspond():
    """First four lines from plain and gzipped versions match."""
    plain_lines = read_first_four_lines(FASTQ_PLAIN)
    gz_lines = read_first_four_lines(FASTQ_GZ)
    assert plain_lines == gz_lines

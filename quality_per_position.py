"""Quality per position statistics for FASTQ files.

Reads a FASTQ file and returns boxplot data per read position, summarizing
the distribution of quality scores across all reads at each position.
"""

import contextlib
import gzip
import math
import sys
from collections import defaultdict
from typing import Iterator, TextIO, TypedDict


class CenterLine(TypedDict):
    min: float
    max: float


class Box(TypedDict):
    q1: float
    q3: float


class Boxplot(TypedDict):
    centerLine: CenterLine
    box: Box
    median: float
    average: float
    name: str


LINES_PER_READ = 4
SEQUENCE_ID = 0
SEQUENCE = 1
OPTIONAL_SEQUENCE_ID = 2
QUALITY_STRING = 3
SANGER_ASCII_QUALITY_SCORE_OFFSET = 33


def open_fastq_stream(filename: str):
    """Open FASTQ stream: file path or '-' for stdin. Never closes stdin."""
    if filename == "-":
        return contextlib.nullcontext(sys.stdin)
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt", encoding="utf-8")
    return open(filename, "r", encoding="utf-8")


def open_fastq_file(filename: str):
    """Open FASTQ file or stdin ('-'). Returns context manager; never closes stdin."""
    return open_fastq_stream(filename)


def read_fastq_records(stream: TextIO) -> Iterator[list[str]]:
    """Yield FASTQ records as lists of LINES_PER_READ lines."""
    while True:
        record = [stream.readline() for _ in range(LINES_PER_READ)]
        if len(record) < LINES_PER_READ or not record[0]:
            return
        yield record


def _parse_quality_scores(quality_string: str, sequence_length: int) -> list[float]:
    """Convert quality string to list of scores (Sanger/Phred+33).

    Raises:
        ValueError: If quality length != sequence length or if any char has ASCII < 33.
    """
    if len(quality_string) != sequence_length:
        raise ValueError(
            f"Quality string length ({len(quality_string)}) does not match "
            f"sequence length ({sequence_length})"
        )
    scores = []
    for char in quality_string:
        q = ord(char) - SANGER_ASCII_QUALITY_SCORE_OFFSET
        if q < 0:
            raise ValueError(
                f"Invalid quality character (ASCII {ord(char)}) at position {len(scores)}; "
                "Sanger encoding requires ASCII >= 33"
            )
        scores.append(float(q))
    return scores


def _compute_position_stats(scores_at_position: list[float], position: int) -> Boxplot:
    """Compute boxplot statistics for a single position across all reads."""
    n = len(scores_at_position)
    if n == 0:
        raise ValueError(f"Position {position} has no scores")
    sorted_scores = sorted(scores_at_position)
    average = sum(scores_at_position) / n
    center_line_min = sorted_scores[0]
    center_line_max = sorted_scores[-1]
    virtual_indices = [p * (n - 1) for p in [0.25, 0.5, 0.75]]
    percentiles = []
    for vi in virtual_indices:
        lower = int(math.floor(vi))
        upper = lower + 1
        if upper >= n:
            percentiles.append(sorted_scores[lower])
        else:
            frac = vi - lower
            percentiles.append(
                (1.0 - frac) * sorted_scores[lower] + frac * sorted_scores[upper]
            )
    q1, median, q3 = percentiles[0], percentiles[1], percentiles[2]
    return Boxplot(
        centerLine={"min": center_line_min, "max": center_line_max},
        box={"q1": q1, "q3": q3},
        median=median,
        average=average,
        name=str(position),
    )


def quality_per_position_boxplot_data(fastq_filename: str) -> list[Boxplot]:
    """Read FASTQ file or stdin ('-') and return one boxplot per read position."""
    position_to_scores: dict[int, list[float]] = defaultdict(list)
    with open_fastq_stream(fastq_filename) as stream:
        for record in read_fastq_records(stream):
            sequence = record[SEQUENCE].strip()
            quality_string = record[QUALITY_STRING].strip()
            scores = _parse_quality_scores(quality_string, len(sequence))
            for i, score in enumerate(scores):
                position_to_scores[i].append(score)
    if not position_to_scores:
        return []
    max_position = max(position_to_scores.keys())
    return [
        _compute_position_stats(position_to_scores[i], i)
        for i in range(max_position + 1)
    ]

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


def _parse_quality_scores(quality_string: str) -> list[int]:
    """Convert quality string to list of integer scores (Sanger/Phred+33).

    Uses a single allocation (list comprehension) instead of repeated append.
    No per-element function calls, so it is faster than an append loop or map.

    Assumes the quality string is valid: length matches the sequence line and all
    characters have ASCII >= 33 (Sanger encoding). Caller must validate input.

    Args:
        quality_string: Raw quality line from FASTQ (after strip).

    Returns:
        list[int]: Integer quality scores (one per base, no conversion to float).
    """
    return [ord(c) - SANGER_ASCII_QUALITY_SCORE_OFFSET for c in quality_string]


def _compute_position_stats(scores_at_position: list[int], position: int) -> Boxplot:
    """Compute boxplot statistics for a single position using histogram-based percentiles.

    Uses a histogram instead of sorting to compute Q1, median, and Q3 efficiently.
    Scores are discrete (typically Phred+33, range 0..94) and lists can be long,
    so O(n) histogram + O(buckets) lookup is preferable to O(n log n) sorting.

    Percentile interpolation uses virtual indices vi = p * (n-1) and linear
    interpolation between the scores at floor(vi) and floor(vi)+1 in sorted order.
    The histogram maps each score to its count; cumulative counts give sorted-order
    positions. We resolve which bucket contains each discrete index, then interpolate.

    Edge cases (see inline comments in _percentile_from_histogram):

    1. Same bucket: lower_idx and upper_idx refer to elements in the same histogram
       bucket (same score). The interpolation formula yields that score as a float
       implicitly, we we'll return a float anyway (desired behavior).

    2. Adjacent buckets: lower and upper fall in consecutive score buckets. Linear
       interpolation between the two bucket scores.

    3. Non-adjacent buckets (empty buckets in between): Lower and upper fall in
       buckets with score difference > 1. We still interpolate linearly between
       the two bucket scores; the result estimates the percentile within the gap.

    4. Percentiles outside [0, n-1]: Clamp indices so we never look beyond the data
       range. For standard p in [0.25, 0.5, 0.75] this only matters when n=1.

    Args:
        scores_at_position: Quality scores at this position (Phred+33, typically 0..94).
        position: 0-based position index (for name and error messages).

    Returns:
        Boxplot with centerLine, box, median, average, name.

    Raises:
        ValueError: If scores_at_position is empty.
    """
    n = len(scores_at_position)
    if n == 0:
        raise ValueError(f"Position {position} has no scores")

    # Build histogram: score -> count. Scores are discrete (Phred+33, typically 0..94).
    hist: dict[int, int] = defaultdict(int)
    for s in scores_at_position:
        hist[int(s)] += 1

    min_score = min(hist.keys())
    max_score = max(hist.keys())
    average = sum(scores_at_position) / n

    # Cumulative count before each score. Only non-empty buckets: fewer iterations
    # in lookup and less memory when scores have gaps (e.g. distinct values 28–42).
    # Elements at sorted indices [cum_before[s], cum_before[s]+hist[s]) have score s.
    scores_sorted = sorted(hist.keys())
    cum_before: dict[int, int] = {}
    cum = 0
    for s in scores_sorted:
        cum_before[s] = cum
        cum += hist[s]

    # Lookup iterates over buckets only (O(b), b = distinct scores), not n.
    # We assume k values increase across calls; resume from last bucket to avoid
    # re-scanning. For n>=3 the six k values are monotonic; for n=2 we may reset.
    _bucket_idx = [0]

    def _score_at_sorted_index(k: int) -> float:
        """Return the score of the element at sorted index k (0-based)."""
        k = max(0, min(k, n - 1))
        # If k fell before current bucket (edge case: n=2), reset
        if _bucket_idx[0] > 0 and k < cum_before[scores_sorted[_bucket_idx[0]]]:
            _bucket_idx[0] = 0
        for i in range(_bucket_idx[0], len(scores_sorted)):
            s = scores_sorted[i]
            if cum_before[s] <= k < cum_before[s] + hist[s]:
                _bucket_idx[0] = i
                return float(s)
        return float(max_score)

    # Percentile interpolation via virtual indices in sorted-order space.
    #
    # Edge cases for lower_idx = floor(vi) and upper_idx = lower_idx + 1:
    #
    # (1) Same bucket: Both indices point to elements in the same histogram bucket
    #     (e.g. many scores equal 40). score_lower == score_upper, so interpolation
    #     yields that discrete score.
    #
    # (2) Adjacent buckets: lower and upper fall in consecutive score buckets
    #     (e.g. 10 and 11). Linear interpolation produces a value between them.
    #
    # (3) Non-adjacent buckets (empty buckets in between): lower in bucket 10,
    #     upper in bucket 25, buckets 11..24 empty. We still interpolate between
    #     10 and 25; the result estimates the percentile within that gap.
    #
    # (4) Percentiles outside [0, n-1]: When upper_idx >= n (e.g. n=1 or vi at
    #     boundary), we use only the lower element. _score_at_sorted_index clamps
    #     its argument to [0, n-1] for safety.
    virtual_indices = [p * (n - 1) for p in [0.25, 0.5, 0.75]]
    percentiles = []
    for vi in virtual_indices:
        lower_idx = int(math.floor(vi))
        upper_idx = lower_idx + 1

        if upper_idx >= n:
            percentiles.append(_score_at_sorted_index(lower_idx))
        else:
            score_lower = _score_at_sorted_index(lower_idx)
            score_upper = _score_at_sorted_index(upper_idx)
            frac = vi - lower_idx
            percentiles.append(
                (1.0 - frac) * score_lower + frac * score_upper
            )

    q1, median, q3 = percentiles[0], percentiles[1], percentiles[2]
    return Boxplot(
        centerLine={"min": float(min_score), "max": float(max_score)},
        box={"q1": q1, "q3": q3},
        median=median,
        average=average,
        name=str(position),
    )


def quality_per_position_boxplot_data(fastq_filename: str) -> list[Boxplot]:
    """Read FASTQ file or stdin ('-') and return one boxplot per read position."""
    position_to_scores: dict[int, list[int]] = defaultdict(list)
    with open_fastq_stream(fastq_filename) as stream:
        for record in read_fastq_records(stream):
            sequence = record[SEQUENCE].strip()
            quality_string = record[QUALITY_STRING].strip()
            scores = _parse_quality_scores(quality_string)
            for i, score in enumerate(scores):
                position_to_scores[i].append(score)
    if not position_to_scores:
        return []
    max_position = max(position_to_scores.keys())
    return [
        _compute_position_stats(position_to_scores[i], i)
        for i in range(max_position + 1)
    ]

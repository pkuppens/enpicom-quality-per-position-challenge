import gzip
import math
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


# Reads in fastq files come in pairs of four
LINES_PER_READ = 4
# First line contains a leading @ followed by the sequency id
SEQUENCE_ID = 0
# Second line contains the sequence
SEQUENCE = 1
# Third line contains a leading + followed by an optional sequence id equal to the one in the first line
OPTIONAL_SEQUENCE_ID = 2
# Fourth line contains a quality string equal in length to the sequence
QUALITY_STRING = 3
"""
A character in the quality string corresponds to a code in the ascii table.
The character which has quality score zero has ascii code 33.
Going from character to ascii code to quality score requires an offset.
"""
SANGER_ASCII_QUALITY_SCORE_OFFSET = 33


def open_fastq_file(filename: str) -> TextIO:
    """Open a FASTQ file, handling plain (.fastq) and gzip-compressed (.fastq.gz) formats.

    Args:
        filename: Path to the FASTQ file (.fastq or .fastq.gz).

    Returns:
        TextIO: A text stream of the file contents (decompressed if gzipped).

    Error handling: Explicitly left out of scope for this assignment.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt", encoding="utf-8")
    return open(filename, "r", encoding="utf-8")


def read_fastq_records(stream: TextIO) -> Iterator[list[str]]:
    """Yield FASTQ records as lists of LINES_PER_READ lines.

    Args:
        stream: Text stream of a FASTQ file.

    Yields:
        list[str]: Each record as 4 lines [sequence_id, sequence, optional_id, quality].
    """
    while True:
        record = [stream.readline() for _ in range(LINES_PER_READ)]
        if len(record) < LINES_PER_READ or not record[0]:
            return
        yield record


def quality_per_position_boxplot_data(fastq_filename: str) -> list[Boxplot]:
    boxplots = []

    # Step 1. Open the fastq file
    # This was explicitly specified as:
    # FASTQ files may be plain text (.fq) or gzip-compressed (.fq.gz), by example. No explicit requirement on this.
    # So we'll create and call a function to open the file, depending on the extension, and return the TextIO stream of the gunzipped content
    with open_fastq_file(fastq_filename) as fastq_file:
        # Since the assignment is in the context of Petabytes scale company, and the main loop includes a performance and memory timing,
        # we'll do some small optimizations. The first will be to stream the FASTQ file, instead of reading it all into memory.
        # Sample file is 1000 lines, 250 records, which is not a lot of data, but still, we'll do it.

        # The elegant way of processing is to create a generator function that yields the lines of the FASTQ file, in blocks of 4 lines.
        # Then parses them into records, and yields the records for further processing (the boxplot calculations)
        # This was explicitly specified as: Reads in fastq files come in pairs of four
        # So we'll read the file in chunks of 4 lines (use LINES_PER_READ const)
        for record in read_fastq_records(fastq_file):
            # This is a sequence identifier in the first line of the record
            sequence_id = record[SEQUENCE_ID].strip()
            # The sequence is in the second line of the record
            sequence = record[SEQUENCE].strip()
            # The optional sequence id is in the third line of the record
            optional_sequence_id = record[OPTIONAL_SEQUENCE_ID].strip()
            # The quality string is in the fourth line of the record
            quality_string = record[QUALITY_STRING].strip()

            # Specification says: Loads quality scores — parses the file and converts every quality string into a list of integer scores, producing one list per read.
            # I'll read as bytes (memory efficient, iterable over values, as we'll see later)
            quality_score_bytes = bytes(
                ord(char) - SANGER_ASCII_QUALITY_SCORE_OFFSET for char in quality_string
            )

            # Now we can calculate the boxplot for this record
            boxplot = calculate_boxplot(
                sequence_id, sequence, optional_sequence_id, quality_score_bytes
            )
            boxplots.append(boxplot)

    return boxplots


def calculate_boxplot(
    sequence_id: str,
    sequence: str,
    optional_sequence_id: str,
    quality_score_bytes: bytes,
) -> Boxplot:
    """Calculate the boxplot for a FASTQ record.

    Args:
        sequence_id (str): The sequence identifier (from the first line of the FASTQ record).
        sequence (str): The sequence string (second line of the FASTQ record).
        optional_sequence_id (str): The optional sequence identifier (third line of the FASTQ record, usually '+').
        quality_score_bytes (bytes): Bytes array of quality scores, where each value represents a Sanger-encoded Phred quality score for each base (fourth line of the FASTQ record, converted).

    Error handling: Explicitly left out of scope for this assignment. E.g. empty sequence, invalid quality scores, etc.
    """
    # We can ignore the first 3 parameters, as they are not used in the boxplot calculation
    # The quality scores are in the bytes iterable
    n_positions = len(quality_score_bytes)

    average = float(sum(quality_score_bytes)) / n_positions

    # Define virtual indices for the percentiles
    virtual_index = [p * (n_positions - 1) for p in [0.25, 0.5, 0.75]]

    quality_score_bytes_sorted = sorted(quality_score_bytes)

    centerLine_min = float(quality_score_bytes_sorted[0])
    centerLine_max = float(quality_score_bytes_sorted[-1])

    # iterate over the virtual indices, and determine the lower and upper index for the percentiles
    percentiles = []
    for p in virtual_index:
        lower_index = math.floor(p)
        upper_index = lower_index + 1

        if upper_index == n_positions:
            # theoretically on low samples or high percentile, we would exceed the array bounds
            percentile = float(quality_score_bytes_sorted[lower_index])
        else:
            # We need to interpolate between the two values
            fraction = p - lower_index
            percentile = (1.0 - fraction) * quality_score_bytes_sorted[
                lower_index
            ] + fraction * quality_score_bytes_sorted[upper_index]

        percentiles.append(percentile)

    # We're done with the centerLine, now we need to calculate the box
    q1 = percentiles[0]
    q3 = percentiles[2]
    median = percentiles[1]

    return Boxplot(
        centerLine={
            "min": centerLine_min,
            "max": centerLine_max,
        },
        box={
            "q1": q1,
            "q3": q3,
        },
        median=median,
        average=average,
        name=sequence_id,
    )

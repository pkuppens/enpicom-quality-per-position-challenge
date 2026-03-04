from typing import TypedDict


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


def quality_per_position_boxplot_data(fastq_filename: str) -> list[Boxplot]:
    boxplots = []

    return boxplots

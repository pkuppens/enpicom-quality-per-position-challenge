"""Measure the cost of FASTQ quality-string validation (none vs per-char loop vs cheap).

Usage: python scripts/benchmark_validation.py   (run from the repo root)
"""

import os
import random
import sys
import tempfile
import time

sys.path.insert(0, os.getcwd())
import quality_per_position as q  # noqa: E402

random.seed(1)
N, L = 100_000, 150
quals = ["".join(chr(random.randint(33, 73)) for _ in range(L)) for _ in range(2000)]
quals = [quals[i % 2000] for i in range(N)]
seqs = ["A" * L] * N


def fast(qs):
    return [ord(c) - 33 for c in qs]


def old(qs, n):
    if len(qs) != n:
        raise ValueError("len")
    out = []
    for c in qs:
        v = ord(c) - 33
        if v < 0:
            raise ValueError("chr")
        out.append(float(v))
    return out


def cheap(qs, n):
    if len(qs) != n:
        raise ValueError("len")
    if min(qs) < "!":
        raise ValueError("chr")
    return [ord(c) - 33 for c in qs]


def t(f):
    best = 1e9
    for _ in range(3):
        s = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - s)
    return best


a = t(lambda: [fast(x) for x in quals])
b = t(lambda: [old(x, L) for x in quals])
c = t(lambda: [cheap(x, L) for x in quals])
print(f"parse only, {N} reads x {L} bp")
print(f"  baseline, no validation: {a:.3f}s")
print(f"  old validation loop:     {b:.3f}s  ({b / a:.1f}x)")
print(
    f"  cheap validation:        {c:.3f}s  ({c / a:.2f}x, +{(c - a) / N * 1e6:.2f} us/read)"
)
# end to end on the generated file
path = os.path.join(tempfile.gettempdir(), "bench_val.fq")
with open(path, "w") as f:
    for i in range(N):
        f.write(f"@r{i}\n{seqs[i]}\n+\n{quals[i]}\n")
e2e = t(lambda: q.quality_per_position_boxplot_data(path))
print(f"end-to-end with the repo's validation: {e2e:.2f}s for {N} reads")
print(
    f"  parse (no validation) is {a / e2e * 100:.0f}% of that; old loop would add {(b - a) / e2e * 100:.0f}%"
)
os.unlink(path)

import sys
from pprint import pprint

from quality_per_position import quality_per_position_boxplot_data


def main():
    fastq_filename = "fixtures/test.fq" if len(sys.argv) == 1 else sys.argv[1]
    boxplots = quality_per_position_boxplot_data(fastq_filename)
    pprint(boxplots)


if __name__ == "__main__":
    import resource
    import sys
    import time

    t0 = time.perf_counter()

    main()

    elapsed = time.perf_counter() - t0
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # KiB on Linux

    print("\n--- profile ---", file=sys.stderr)
    print(f"runtime:      {elapsed:.3f}s", file=sys.stderr)
    print(f"peak memory:  {peak_kb / 1024:.2f} MiB", file=sys.stderr)

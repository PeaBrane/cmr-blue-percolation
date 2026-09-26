"""Reproduce the surgery-check counts of the diluted-model manuscript (optional; needs networkx).

The manuscript's appendix reports that a literal implementation of the surgery
map, check_surgery.py in this directory, found no failure on 6546 (+)-pivotal
instances in d = 2 (739 resolved in Step 1; 2285, 1892 and 1630 of types N, B
and S) and on 2197 in d = 3 (442, 793, 384 and 578). These totals are the sums
of the five runs below. The research logs did not record the seeds; each seed
was recovered by rerunning check_surgery.py and matching the first line of the
log, the cumulative counts after 10 samples.

    d  samples  seed  mode    L   research log
    2  200      13    path    10  surgery_d2_path_L10.out
    2  300      12    path    16  surgery_d2_path_L16.out
    2  100      11    random  10  surgery_d2_random_L10.out
    3  60       15    path    10  surgery_d3_path_L10.out
    3  100      14    path    12  surgery_d3_path_L12.out

By default the three d = 2 runs stop after 10 samples (10 to 20 CPU
seconds), and their counts are compared with the first line of the research
log.
--slow adds the first 10 samples of the two d = 3 runs, which take minutes
(about 15 CPU seconds for the first 2 samples at L = 10). --full completes
all five runs, compares their final counts with the last line of each log and
checks the manuscript's totals; it takes hours. Every run also passes all
assertions of check_surgery.py. This is a finite sanity check, not a proof
input.
"""

import argparse
import ast
from pathlib import Path
import subprocess
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
KEYS = ("step1", "N", "B", "S", "pivotal_edges")

# (d, samples, seed, mode, L, counts after 10 samples, final counts), counts as in KEYS.
RUNS = [
    (2, 200, 13, "path", 10, (9, 53, 7, 50, 119), (245, 911, 208, 723, 2087)),
    (2, 300, 12, "path", 16, (11, 55, 65, 28, 159), (471, 1115, 1673, 860, 4119)),
    (2, 100, 11, "random", 10, (0, 27, 4, 8, 39), (23, 259, 11, 47, 340)),
    (3, 60, 15, "path", 10, (25, 45, 6, 44, 120), (146, 303, 77, 230, 756)),
    (3, 100, 14, "path", 12, (36, 66, 37, 42, 181), (296, 490, 307, 348, 1441)),
]
SLOW = {(3, 10), (3, 12)}
# Totals displayed in the manuscript, as in KEYS.
TOTALS = {2: (739, 2285, 1892, 1630, 6546), 3: (442, 793, 384, 578, 2197)}


def run(d, samples, seed, mode, L):
    command = [sys.executable, str(HERE / "check_surgery.py"), str(d), str(samples), str(seed), mode, str(L)]
    out = subprocess.run(command, capture_output=True, text=True, check=True).stdout.splitlines()
    assert out[-1] == "ALL CHECKS PASSED", out[-3:]
    assert out[-2].startswith("FINAL "), out[-2]
    stats = ast.literal_eval(out[-2][len("FINAL "):])
    assert stats["samples"] == samples
    return tuple(stats[key] for key in KEYS)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--slow", action="store_true", help="also check the first 10 samples of the two d = 3 runs")
    group.add_argument("--full", action="store_true", help="complete every run and check the manuscript's totals")
    args = parser.parse_args()
    totals = {2: [0] * len(KEYS), 3: [0] * len(KEYS)}
    for d, samples, seed, mode, L, first, final in RUNS:
        if (d, L) in SLOW and not (args.slow or args.full):
            print(f"not run: d={d} {mode} L={L} seed={seed} (use --slow or --full)", flush=True)
            continue
        start = time.monotonic()
        n = samples if args.full else 10
        counts = run(d, n, seed, mode, L)
        assert counts == (final if args.full else first), (d, seed, mode, L, counts)
        totals[d] = [a + b for a, b in zip(totals[d], counts)]
        print(f"PASS d={d} {mode} L={L} seed={seed}, {n} samples: {counts[-1]} (+)-pivotal edges, "
              f"counts equal to the research log ({time.monotonic() - start:.1f}s)", flush=True)
    if args.full:
        assert {d: tuple(v) for d, v in totals.items()} == TOTALS
        print(f"PASS manuscript totals: {TOTALS[2][-1]} instances in d = 2 and {TOTALS[3][-1]} in d = 3")


if __name__ == "__main__":
    main()

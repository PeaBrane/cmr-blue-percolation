"""Full reproduction of the sparse-insulation window checks (d = 2 box, d = 3 cube, octahedral rule).

Runs, each in a child process of this interpreter, and compares every output line with expected_outputs/:
 1. the default light check check_window_general.py;
 2. oct_check.py for d = 3 at L = 6 (bracket R8), d = 4 at L = 5 (brackets R8 and old; about 45 s and
    3.8 GB each) and d = 5 at L = 5 (bracket R8, without the bulk clique; about 6.5 min, 0.4 GB): every
    rule design satisfies the design-lemma hypotheses with the costs of its case, the same-edge key pairs
    behave as the arc-reconstruction lemma requires, the realizable-key sum is at most K_d, and (d = 4) the
    exact bulk clique value is at most K_d on the stated bracket;
 3. audit_extra.py on both rule tables (needs networkx): own maximal cliques equal networkx.find_cliques,
    the anchor-variant diagnosis, type-N options contain 0, mutation tests (sanity, not a proof input);
 4. with --e2e, the end-to-end sanity runs of e2e_check.py whose parameters were recorded (not a proof
    input): d = 2 at L = 5, 6, 7 (600 random states each) and d = 3 at L = 5 (seed 21, first 180 states).

Needs Python 3.10+; step 3 needs networkx. usage: python reproduce_window_general.py [--jobs N] [--e2e]
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
EXP = HERE / "expected_outputs"
T3, T2 = "params/lp_3_cube_rho_256.json", "params/lp_2_box_rho.json"

# (label, script and arguments, expected output, mode): mode "all" compares the whole output, "prefix" compares
# the output up to the recorded last line (the recorded run was stopped early), ignoring the final summary line.
OCT = [
    ("oct_check d=3 L=6 R8", ["oct_check.py", 3, 6, 1, 1, "R8"], "oct_check_d3_L6_R8.out", "all"),
    ("oct_check d=4 L=5 R8", ["oct_check.py", 4, 5, 1, 1, "R8"], "oct_check_d4_L5_R8.out", "all"),
    ("oct_check d=4 L=5 old", ["oct_check.py", 4, 5, 1, 1, "old"], "oct_check_d4_L5_old.out", "all"),
    ("oct_check d=5 L=5 R8", ["oct_check.py", 5, 5, 1, 0, "R8"], "oct_check_d5_L5_R8.out", "all"),
]
AUDIT = [
    ("audit_extra T_3", ["audit_extra.py", T3], "audit_extra_d3.out", "all"),
    ("audit_extra T_2", ["audit_extra.py", T2], "audit_extra_d2.out", "all"),
]
E2E = [
    ("e2e d=2 L=6 seed 11", ["e2e_check.py", T2, 6, 600, 11, 0.40, 0.5], "e2e/e2e_d2_L6_s11.out", "all"),
    ("e2e d=2 L=7 seed 12", ["e2e_check.py", T2, 7, 600, 12, 0.45, 1.0], "e2e/e2e_d2_L7_s12.out", "all"),
    ("e2e d=2 L=5 seed 13", ["e2e_check.py", T2, 5, 600, 13, 0.35, 0.3], "e2e/e2e_d2_L5_s13.out", "all"),
    ("e2e d=3 L=5 seed 21", ["e2e_check.py", T3, 5, 180, 21, 0.20, 0.5], "e2e/e2e_d3_L5_s21.out", "prefix"),
]


def run(job):
    label, cmd, expected, mode = job
    start = time.monotonic()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1")
    proc = subprocess.run([sys.executable, *map(str, cmd)], cwd=HERE, env=env, capture_output=True, text=True)
    seconds = round(time.monotonic() - start, 1)
    if proc.returncode:
        return label, False, seconds, f"exit status {proc.returncode}: {proc.stderr[-2000:]}"
    got, want = proc.stdout, (EXP / expected).read_text()
    if mode == "prefix":
        got = "".join(line for line in got.splitlines(keepends=True) if not line.startswith("final "))
    ok = got == want
    return label, ok, seconds, "" if ok else f"output differs from expected_outputs/{expected}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--jobs", type=int, default=min(6, os.cpu_count() or 1),
                    help="parallel child processes (the two d = 4 runs need 3.8 GB each)")
    ap.add_argument("--e2e", action="store_true", help="also run the recorded end-to-end sanity runs")
    args = ap.parse_args()
    t0 = time.monotonic()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    light = subprocess.run([sys.executable, "check_window_general.py"], cwd=HERE, env=env, check=True,
                           capture_output=True, text=True).stdout
    print(light.split("\n{", 1)[0], flush=True)
    report = {"python": sys.version.split()[0], "jobs": args.jobs,
              "light_check": json.loads("{" + light.split("\n{", 1)[1])}
    jobs = list(OCT)
    if importlib.util.find_spec("networkx") is None:
        raise SystemExit("reproduce_window_general.py needs networkx for audit_extra.py")
    jobs += AUDIT
    if args.e2e:
        jobs += E2E
    # longest first
    order = {"oct_check d=5 L=5 R8": 0}
    jobs.sort(key=lambda j: order.get(j[0], 1))
    with ThreadPoolExecutor(args.jobs) as ex:
        results = list(ex.map(run, jobs))
    failed = [r for r in results if not r[1]]
    for label, ok, seconds, msg in results:
        print(f"{'PASS' if ok else 'FAIL'} {label} ({seconds}s){': ' + msg if msg else ''}", flush=True)
    report["runs"] = {label: {"ok": ok, "seconds": seconds} for label, ok, seconds, _ in results}
    report["total_seconds"] = round(time.monotonic() - t0, 1)
    if failed:
        raise SystemExit(f"FAILED: {len(failed)} of {len(results)} runs")
    print(f"PASS window_general full reproduction ({report['total_seconds']}s)")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

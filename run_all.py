"""Run the certificates; optionally add supplementary and optional-dependency checks.

The default run and --supplementary need only the Python standard library.
--optional-deps adds finite sanity checks and the code-base-A global
crosschecks, which need mpmath, SymPy, NumPy or networkx; a check whose
package cannot be imported is skipped with a message. Every script runs in a
child process of the interpreter that runs this file. The full reproductions
(macrostep/reproduce_engine.py, diluted_uniqueness/window_d2/reproduce_window_d2.py,
diluted_uniqueness/window_general/reproduce_window_general.py) are separate
commands; see README.md.
"""

import argparse
from importlib.util import find_spec
from pathlib import Path
import subprocess
import sys
import time

DEFAULT = [
    "verify_signed_star.py",
    "verify_signed_star_independent.py",
    "verify_physical_separation.py",
    "overlap_revealed/certify.py",
    "single_floor/verify_single_floor.py",
    "diluted_uniqueness/gd_constants.py",
    "diluted_uniqueness/window_d2/check_window_d2.py",
    "diluted_uniqueness/window_general/check_window_general.py",
    "macrostep/certify_macrostep.py",
    "overlap_revealed/noisy_cavity/certify_noisy_cavity.py",
]
SUPPLEMENTARY = [
    "checks/smallstar_audit.py",
    "checks/numerator_audit.py",
    "checks/check_susceptibility.py",
    "overlap_revealed/crosscheck/crosscheck_local.py",
    "single_floor/check_lemmas.py",
    "single_floor/check_side_results.py",
    "diluted_uniqueness/check_typeS_rule.py",
    "overlap_revealed/noisy_cavity/crosscheck/nesting_exact.py",
]
# (script and arguments, packages it imports)
OPTIONAL = [
    (["overlap_revealed/crosscheck/crosscheck_global.py"], ["mpmath"]),
    (["diluted_uniqueness/optional/check_classical_window.py"], ["mpmath"]),
    (["diluted_uniqueness/optional/check_decimation_exact.py"], ["sympy"]),
    (["diluted_uniqueness/optional/check_gibbs_reduction.py"], ["numpy"]),
    (["diluted_uniqueness/optional/surgery_counts.py"], ["networkx"]),
    (["overlap_revealed/noisy_cavity/crosscheck/crosscheck_global_A.py"], ["mpmath"]),
    (["overlap_revealed/noisy_cavity/crosscheck/verify_pair_bruteforce.py"], ["numpy"]),
    (["overlap_revealed/noisy_cavity/crosscheck/verify_holley_bruteforce.py"], ["numpy"]),
    (["overlap_revealed/noisy_cavity/crosscheck/lemmaN_exact.py"], ["numpy"]),
    (["overlap_revealed/noisy_cavity/crosscheck/pairing_exact.py"], ["numpy"]),
]


def main():
    if sys.flags.optimize:
        raise SystemExit("Verification requires assertions: run Python without -O.")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--supplementary", action="store_true",
                        help="also run the independent finite-model and crosscheck programs (standard library)")
    parser.add_argument("--optional-deps", action="store_true",
                        help="also run the checks that need mpmath, SymPy, NumPy or networkx, when importable")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    runs = [[script] for script in DEFAULT]
    if args.supplementary:
        runs += [[script] for script in SUPPLEMENTARY]
    skipped = []
    if args.optional_deps:
        for command, packages in OPTIONAL:
            missing = [name for name in packages if find_spec(name) is None]
            if missing:
                skipped.append((" ".join(command), missing))
            else:
                runs.append(command)
    for command in runs:
        print(f"\nRunning {' '.join(command)}", flush=True)
        start = time.monotonic()
        subprocess.run([sys.executable, str(root / command[0]), *command[1:]], cwd=root, check=True)
        print(f"(finished in {time.monotonic() - start:.1f}s)", flush=True)
    for command, missing in skipped:
        print(f"\nSKIPPED optional check {command}: cannot import {', '.join(missing)} with {sys.executable}")
    print(f"\nPASS: {len(runs)} verification scripts completed.")


if __name__ == "__main__":
    main()

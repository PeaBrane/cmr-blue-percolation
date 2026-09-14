"""Run the paper's certificates, optionally including finite-model checks."""

import argparse
from pathlib import Path
import subprocess
import sys


def main():
    if sys.flags.optimize:
        raise SystemExit("Verification requires assertions: run Python without -O.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplementary", action="store_true",
                        help="also run the independent finite-model checks")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    scripts = ["verify_signed_star.py", "verify_signed_star_independent.py",
               "verify_physical_separation.py"]
    if args.supplementary:
        scripts += ["checks/smallstar_audit.py", "checks/numerator_audit.py",
                    "checks/check_susceptibility.py"]
    for script in scripts:
        print(f"\nRunning {script}", flush=True)
        subprocess.run([sys.executable, str(root / script)], cwd=root, check=True)
    print(f"\nPASS: {len(scripts)} verification scripts completed.")


if __name__ == "__main__":
    main()

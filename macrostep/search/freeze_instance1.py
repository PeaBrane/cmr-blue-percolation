"""Re-derive the frozen instance-1 local inputs of ../params.json (optional; needs SciPy).

Nothing here is needed for verification: ../certify_macrostep.py reads only the frozen exact rationals. For the two
instance-1 points (d, t) = (9, 7/50) and (8, 3/20) this program repeats the floating-point choices of the second
implementation of Section 4.7 of the first manuscript (the search functions of
../../overlap_revealed/search/freeze_params.py):

* the bond-floor tangent points c_U: Nelder-Mead over the family
  c_U = kappa * omega(U)^lambda * (2 cosh^m beta)^(1 - lambda), each c_U rounded to 10 decimals;
* the Holley tangent weights w_k: the best point of {0} U {i/200 : 1 <= i <= 180}, rounded to 4 decimals.

The bond floor p, the line (g_K, g_h) and xbar are the certified inputs of Table 4 and are copied from
../params.json. Any positive tangent point and any weight in [0, 1] are admissible, so the search cannot affect
soundness.

The frozen values were produced by this procedure with SciPy 1.16.0 and NumPy 2.2.6 on Apple Silicon, and this
program reproduces them there (MATCH). With SciPy 1.18.1 and NumPy 2.5.3 on x86-64 Linux, Nelder-Mead stops at a
slightly different point, and the program reports MISMATCH: the tangent points differ from about the eighth
significant digit, while the weights agree.

    python freeze_instance1.py            # recompute and compare with ../params.json
    python freeze_instance1.py --print    # print the recomputed JSON blocks
"""
import argparse
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "overlap_revealed" / "search"))
import freeze_params as FP  # noqa: E402  (imports ../../overlap_revealed/certify.py as FP.exact)

PARAMS = HERE.parent / "params.json"


def derive(inst):
    d, t, p = inst["d"], Fr(inst["t"]), Fr(inst["p"])
    M = FP.Model(d, t)
    kappa, lam = FP.search_family(M, p)
    cU = FP.tangent_family(M, p, kappa, lam)
    weights = [FP.decimal(FP.choose_wk(M, k), 4) for k in range(M.m + 1)]
    return {"pB_tangent_points": {"form": "direct", "values": {str(U): FP.decimal(c, 10) for U, c in sorted(cU.items())}},
            "holley_weights": weights,
            "provenance": {"pB_tangent_points": f"Nelder-Mead family kappa={kappa!r}, lambda={lam!r}; c_U rounded to "
                                                f"10 decimals (search functions of overlap_revealed/search/freeze_params.py)",
                           "holley_weights": "best point of {0} U {i/200: 1<=i<=180} for the floating numerator bound"}}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--print", action="store_true")
    a = ap.parse_args()
    params = json.loads(PARAMS.read_text())
    same = True
    for key, inst in params["instances"].items():
        if inst["local_inputs"] != "tangent":
            continue
        new = derive(inst)
        if a.print:
            print(json.dumps({key: new}, indent=1))
        ok = (new["pB_tangent_points"] == inst.get("pB_tangent_points") and new["holley_weights"] == inst.get("holley_weights"))
        same &= ok
        print(f"{key}: recomputed tangent points and weights {'equal' if ok else 'DIFFER FROM'} the frozen ones")
    print("MATCH" if same else "MISMATCH (platform-dependent float search; any admissible value is valid)")


if __name__ == "__main__":
    main()

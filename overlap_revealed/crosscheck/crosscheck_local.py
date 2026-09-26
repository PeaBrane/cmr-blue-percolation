"""Code base A: independent exact recomputation of the local inputs (standard library only).

For every row of ../params.json this program uses the part-1 exact checkers
(certify_pB.py, certify_odds.py), which were written independently of
../certify.py, with the frozen tangent points and weights:

1. bond floor: every signed-residual vertex sum Gamma_j <= 0;
2. Holley line in all 2d+1 environments, with the balanced-maximum lemma for
   D_k re-checked;
3. K' <= |h| and 2dK' < 1 as g_K g_h <= 1 and g_K^(2d) < 7 < e^2, and the
   mean-field root test in its exact power form, at the four-decimal point
   xbar = ceil(mbar * 10^4) / 10^4 >= mbar;
4. the ordered local rationals (Delta_j, aligned witness, N_k, D_k) hash to
   local_vector_sha256 recorded by ../certify.py in
   ../../expected/overlap_revealed.json, so both code bases produce identical
   exact rationals.

It first reruns part 1's headline certificate (certify_local.py) as
published in the research notes, and confirms that its tangent points and
weights are the frozen headline inputs.
"""

from fractions import Fraction as Fr
from math import ceil
from pathlib import Path
import hashlib
import json
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import certify_local  # noqa: E402
from certify_odds import denominator, make, numerator  # noqa: E402
from certify_pB import aligned_witness, certify  # noqa: E402

PARAMS = HERE.parent / "params.json"
EXPECTED = HERE.parent.parent / "expected" / "overlap_revealed.json"


def tangent_points(row, w):
    spec = row["pB_tangent_points"]
    values = {int(U): Fr(v) for U, v in spec["values"].items()}
    if spec["form"] == "direct":
        return values
    assert spec["form"] == "rho_times_omega"
    return {U: v * (w ** ((1 + U) // 2) + w ** (-((1 + U) // 2))) for U, v in values.items()}


def check_row(row):
    d, t, p = row["d"], Fr(row["t"]), Fr(row["p"])
    gK, gh, mbar = Fr(row["g_K"]), Fr(row["g_h"]), Fr(row["mbar"])
    m, w, a, C, ch2 = make(d, t)
    ok, gam = certify(d, t, p, cU=tangent_points(row, w), verbose=False)
    assert ok and all(g <= 0 for g in gam)
    witness = aligned_witness(d, t)
    vector = list(gam) + [witness]
    for k in range(m + 1):
        N, _ = numerator(k, m, a, ch2, Fr(row["holley_weights"][k]))
        D, dv = denominator(k, m, a, ch2)
        assert D == dv[(m - k) // 2]
        S = 2 * k - m
        assert C ** S * N / D >= gh * gK ** S, (row["label"], k)
        vector += [N, D]
    assert gh < 1 and gK * gh <= 1 and gK ** m < 7
    xbar = Fr(ceil(mbar * 10 ** 4), 10 ** 4)
    zx = m * xbar
    P, Q = zx.numerator, zx.denominator
    assert (1 / gh) ** Q * gK ** P < ((1 + xbar) / (1 - xbar)) ** Q
    print(f"PASS A-local {row['label']}: max Gamma_j = {float(max(gam)):+.4e}, Holley line in {m + 1}"
          f" environments, root test at xbar = {xbar}", flush=True)
    return vector


def main():
    ok1, gam, _, _ = certify_local.pB_certificate()
    ok2, _ = certify_local.odds_certificate(verbose=False)
    ok3, rho = certify_local.density_certificate()
    assert ok1 and ok2 and ok3
    rows = json.loads(PARAMS.read_text())["rows"]
    head = rows[0]
    assert head["label"] == "d12-t3/25" and head["pB_tangent_points"]["form"] == "rho_times_omega"
    assert {int(U): Fr(v) for U, v in head["pB_tangent_points"]["values"].items()} == certify_local.RHO_U
    assert [Fr(w) for w in head["holley_weights"]] == [certify_local.W_K[k] for k in range(25)]
    print(f"PASS A-local headline as in the research log: max Gamma_j = {float(max(gam)):+.4e}, rho > {float(rho)};"
          f" its tangent points and weights are the frozen ones")
    vector = []
    for row in rows:
        vector += check_row(row)
    local_digest = hashlib.sha256("\n".join(str(v) for v in vector).encode()).hexdigest()
    expected = json.loads(EXPECTED.read_text())["local_vector_sha256"]
    assert local_digest == expected, (local_digest, expected)
    print(f"PASS code bases A and B agree on all {len(vector)} local rationals: sha256 {local_digest}")


if __name__ == "__main__":
    main()

"""First implementation: independent evaluation of the global criterion (optional; needs mpmath).

For every row of ../params.json this program runs the global criterion of this
directory (criterion.py, green.py, covariance.py, evaluate_part1.py), written
independently of ../certify.py. It encloses K' = log(g_K)/2 and
|h| = -log(g_h)/2 in outward-rounded 200-bit mpmath intervals, takes
tbar = K' (not tanh K'), finds its own mbar by interval bisection and sums the
resolvent to radius 60. It evaluates the basic Score, Score_c (the manuscript's
shared-edge lemma, Lemma 4.27) and, where 2dK' <= 0.35, Score'_c (Lemma 4.28,
which Table 2 does not use).

Table 2 of the manuscript shows the worse of the two implementations, so this
one must satisfy every displayed bound: rho_-, eta', Score, Score_c and theta_* from
../params.json, and the reference Score'_c bounds and their theta below. Assertions decide
every comparison in exact rational arithmetic on the interval endpoints.

    python crosscheck_global.py [LABEL ...]   # all rows, or only the rows with these labels
"""

from fractions import Fraction as Fr
from pathlib import Path
import json
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from covariance import cov_lower_sharp  # noqa: E402
from criterion import evaluate, iv, ivq, q_lo  # noqa: E402
from evaluate_part1 import enclose  # noqa: E402

PARAMS = HERE.parent / "params.json"


def dn(x, k):
    return f"{float(Fr((x.numerator * 10 ** k) // x.denominator, 10 ** k)):.{k}f}"


def up(x, k):
    return f"{float(Fr(-((-x.numerator * 10 ** k) // x.denominator), 10 ** k)):.{k}f}"

# Reference bounds Score'_c <= and the corresponding theta_* >= at the rows of Table 2 (Lemma 4.28; these
# values are not displayed in the manuscript, which does not use Lemma 4.28 for Table 2).
LEMMA_36 = {"d12-t3/25": ("0.892753", "0.04188"), "d12-t11/100": ("0.915188", "0.03630"),
            "d12-t23/200": ("0.900825", "0.04073"), "d12-t1/8": ("0.889209", "0.04078"),
            "d12-t13/100": (None, None), "d11-t13/100": ("0.934120", "0.02392"),
            "d11-t3/25": ("0.946226", "0.02185"), "d10-t27/200": ("0.987498", "0.00449"),
            "d10-t13/100": ("0.991486", "0.00325")}


def row_bounds(row):
    d, p, gK, gh = row["d"], Fr(row["p"]), Fr(row["g_K"]), Fr(row["g_h"])
    Klo, Khi = enclose(gK)
    hlo, hhi = enclose(gh)
    Hlo, Hhi = -hhi, -hlo
    assert Khi <= Hlo and 2 * d * Khi < 1
    o = evaluate(d, Khi, Hhi, p, verbose=False)
    rho, eta = o["rho"], o["eta"]
    ch = lambda x: (iv.exp(x) + iv.exp(-x)) / 2
    sh = lambda x: (iv.exp(x) - iv.exp(-x)) / 2
    den = (iv.exp(ivq(Khi)) * ch(2 * ivq(Hhi + (2 * d - 1) * Khi))
           + iv.exp(-ivq(Klo)) * ch(2 * ivq((2 * d - 1) * Khi)))
    c = q_lo(2 * sh(2 * ivq(Klo)) / den ** 2)
    assert rho ** 2 >= c / 4
    score_c = 1 / (d * (rho + c / (4 * rho)) * p) + o["term2"] + o["term3"]
    score_s = None
    if 2 * d * Khi <= Fr(35, 100):
        c2, _ = cov_lower_sharp(d, Klo, Khi, Hhi, o["mbar"], Khi)
        assert rho ** 2 >= c2 / 4
        score_s = 1 / (d * (rho + c2 / (4 * rho)) * p) + o["term2"] + o["term3"]
    theta = lambda s: Fr(99, 100) * rho * (1 - eta) * (1 - s) / (1 - Fr(1, d))
    return rho, eta, o["score"], score_c, score_s, theta


def main():
    rows = json.loads(PARAMS.read_text())["rows"]
    labels = sys.argv[1:]
    unknown = set(labels) - {row["label"] for row in rows}
    assert not unknown, f"unknown row labels: {sorted(unknown)}"
    for row in rows:
        if labels and row["label"] not in labels:
            continue
        start = time.monotonic()
        rho, eta, score, score_c, score_s, theta = row_bounds(row)
        used = score if row["criterion"] == "basic" else score_c
        assert eta < 1 and used < 1
        pub = {key: Fr(value) for key, value in row["published"].items()}
        assert rho >= pub["rho_minus_lower"] and eta <= pub["eta_prime_upper"]
        assert score <= pub["score_upper"] and score_c <= pub["score_c_upper"]
        assert theta(used) >= pub["theta_star_lower"]
        shown_s, shown_theta = LEMMA_36[row["label"]]
        assert (score_s is None) == (shown_s is None)
        if score_s is not None:
            assert score_s <= Fr(shown_s) and score_s < 1 and theta(score_s) >= Fr(shown_theta)
        extra = f", Score'_c <= {up(score_s, 6)}" if score_s is not None else ", Score'_c n/a"
        print(f"PASS first-implementation global {row['label']}: rho_- >= {dn(rho, 7)}, eta' <= {up(eta, 6)},"
              f" Score <= {up(score, 6)}, Score_c <= {up(score_c, 6)}{extra},"
              f" theta_* >= {dn(theta(used), 6)} ({time.monotonic() - start:.1f}s)", flush=True)
    scope = f"at the {len(labels)} selected rows" if labels else "at every row"
    print(f"PASS the first implementation satisfies every published global bound {scope}")


if __name__ == "__main__":
    main()

"""Code base A at the five noisy-cavity rows (optional; needs mpmath).

The published global numbers of the second d = 9 proof are the worse of two code bases. The default checker
../certify_noisy_cavity.py evaluates code base B (../../certify.py) and the second implementation of Lemma 3.6
(../lemma36.py). This program evaluates code base A (../../crosscheck/: 200-bit mpmath intervals, tbar = K', its own
mbar by interval bisection) at the same rows, with Lemma 3.6 where 2dK' <= 35/100 (../../crosscheck/covariance.py) and
Lemma 3.6' otherwise (covariance_36prime.py), and asserts in exact rational arithmetic that code base A satisfies
every worse-of-two display of ../params.json (rho_-, eta', Score, Score_c, Score'_c, theta_*, theta') and the
code-base-A columns of the research notes (noisy-cavity write-up Sec. 7.2).

    python crosscheck_global_A.py [LABEL ...]      # all rows, or only these labels (e.g. t=3/20)
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
sys.path.insert(0, str(HERE.parent.parent / "crosscheck"))

from covariance import cov_lower_sharp  # noqa: E402  Lemma 3.6 (asserts 2dK' <= 35/100)
from criterion import evaluate, iv, ivq, q_lo  # noqa: E402
from evaluate_part1 import enclose  # noqa: E402
import covariance_36prime as COV36  # noqa: E402  Lemma 3.6'

PARAMS = json.loads((HERE.parent / "params.json").read_text())
WORSE_OF_TWO = {"rho_minus", "eta_prime", "score", "score_c", "score_c_prime", "theta_star", "theta_prime"}
# Code-base-A columns of the noisy-cavity write-up, Sec. 7.2: c'_cov >= (A), Score'_c <= (A), theta' >= (A).
A_COLUMNS = {"t=7/50": ("0.0162550", "0.99754560", "0.000950073"),
             "t=29/200": ("0.0171782", "0.98757016", "0.004544209"),
             "t=3/20": ("0.0180554", "0.98241246", "0.006000179"),
             "t=31/200": ("0.0188910", "0.98119283", "0.005891535"),
             "t=4/25": ("0.0196274", "0.98511541", "0.004182459")}
# Components at t = 3/20 that the notes attribute to code base A (the larger of the two code bases).
A_COMPONENTS = {"U1": ("upper", "0.0574151"), "SU": ("upper", "11.06112"),
                "rho_c_prime": ("lower", "0.4411173"), "nine_rho_c_prime_p": ("lower", "1.0675481")}


def dn(x, k):
    return f"{float(Fr((x.numerator * 10 ** k) // x.denominator, 10 ** k)):.{k}f}"


def up(x, k):
    return f"{float(Fr(-((-x.numerator * 10 ** k) // x.denominator), 10 ** k)):.{k}f}"


def row_values(row):
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
    if 2 * d * Khi <= Fr(35, 100):
        lemma, (c2, _) = "3.6", cov_lower_sharp(d, Klo, Khi, Hhi, o["mbar"], Khi)
        Lstar = 40
    else:
        lemma, (c2, _) = "3.6'", COV36.cov_lower_sharp(d, Klo, Khi, Hhi, o["mbar"], Khi)
        Lstar = COV36.cov_lower_sharp.Lstar
        assert Lstar <= 41
    assert rho ** 2 >= c2 / 4
    rho_cp = rho + c2 / (4 * rho)
    score_cp = 1 / (d * rho_cp * p) + o["term2"] + o["term3"]
    theta = lambda s: Fr(99, 100) * rho * (1 - eta) * (1 - s) / (1 - Fr(1, d))
    used = score_c if row["score_used"] == "score_c" else score_cp
    assert eta < 1 and used < 1 and score_cp < 1
    return dict(rho_minus=rho, eta_prime=eta, score=o["score"], score_c=score_c, score_c_prime=score_cp,
                theta_star=theta(used), theta_prime=theta(score_cp), c_cov_prime=c2, rho_c_prime=rho_cp,
                nine_rho_c_prime_p=d * rho_cp * p, U1=o["U1"], SU=o["SU"], lemma=lemma, Lstar=Lstar)


def main():
    labels = sys.argv[1:]
    unknown = set(labels) - {row["label"] for row in PARAMS["rows"]}
    assert not unknown, f"unknown row labels: {sorted(unknown)}"
    checked = 0
    for row in PARAMS["rows"]:
        if labels and row["label"] not in labels:
            continue
        start = time.monotonic()
        v = row_values(row)
        for item in PARAMS["displays"]:
            if item["scope"] != row["label"] or item["quantity"] not in WORSE_OF_TWO:
                continue
            if row["score_used"] == "score_c_prime" and item["quantity"] == "theta_star":
                continue                   # at t = 7/50 theta_* is theta' (Score_c > 1 there)
            x, shown = v[item["quantity"]], Fr(item["value"])
            if item["direction"] == "upper":
                assert x <= shown, item
            elif item["direction"] == "lower":
                assert x >= shown, item
            checked += 1
        cA, sA, tA = (Fr(s) for s in A_COLUMNS[row["label"]])
        assert v["c_cov_prime"] >= cA and v["score_c_prime"] <= sA and v["theta_prime"] >= tA
        checked += 3
        if row["label"] == "t=3/20":
            for key, (direction, shown) in A_COMPONENTS.items():
                assert (v[key] <= Fr(shown)) if direction == "upper" else (v[key] >= Fr(shown)), (key, shown)
                checked += 1
        print(f"PASS A-global {row['label']}: rho_- >= {dn(v['rho_minus'], 7)}, eta' <= {up(v['eta_prime'], 7)}, "
              f"Score <= {up(v['score'], 7)}, Score_c <= {up(v['score_c'], 7)}, Score'_c <= {up(v['score_c_prime'], 8)} "
              f"(Lemma {v['lemma']}, L_* = {v['Lstar']}), theta_* >= {dn(v['theta_star'], 6)}, "
              f"theta' >= {dn(v['theta_prime'], 9)} ({time.monotonic() - start:.1f}s)", flush=True)
    print(f"PASS code base A satisfies {checked} displayed bounds")


if __name__ == "__main__":
    main()

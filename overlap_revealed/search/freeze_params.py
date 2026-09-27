"""Re-derive the frozen search outputs in ../params.json (optional; needs SciPy).

Nothing here is needed for verification: ../certify.py reads only the frozen
exact rationals and never calls this program. Floating point only chooses
admissible parameters, following the second of the two implementations
described in Section 4.7 of the manuscript:

* p_B tangent points. For the eight rows marked "direct", Nelder-Mead chooses
  (kappa, lambda) in the family
  c_U = kappa * omega(U)^lambda * (2 cosh^m beta)^(1 - lambda), and every c_U
  is rounded to 10 decimals. The row (12, 3/25) keeps the first
  implementation's values c_U = rho_U * omega(U) (Proposition 4.7).
* Holley tangent weights w_k: the best point of the grid
  {0} U {i/200 : 1 <= i <= 180} for the floating numerator bound. The
  row (12, 3/25) keeps the first implementation's weights
  (crosscheck/certify_local.py), for which the manuscript displays the slack.
* mbar: the smallest multiple of 1e-7 above the floating mean-field root
  that passes the exact root test.

Any positive tangent point and any weight in [0, 1] are admissible, and every
mbar that passes the exact test is valid. Small platform-dependent
differences in this search therefore cannot invalidate anything. The
row-level inputs (d, t, p, g_K, g_h) come from choose_rows.py and are copied
below, as are the displayed numbers of the manuscript.

    python freeze_params.py            # recompute; compare with ../params.json
    python freeze_params.py --write P  # write a fresh parameter file to P
"""

import argparse
import json
import sys
from fractions import Fraction as Fr
from math import atanh, comb, cosh, log, tanh
from pathlib import Path

from scipy.optimize import minimize

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import certify as exact  # noqa: E402

PARAMS = HERE.parent / "params.json"

# Tangent multipliers rho_U of the first implementation at (12, 3/25) (Proposition 4.7): c_U = rho_U * omega(U).
HEADLINE_RHO_U = {1: "1.080881", 3: "1.042249", 5: "0.984864", 7: "0.916142", 9: "0.842612",
                  11: "0.768993", 13: "0.698220", 15: "0.631872", 17: "0.570633", 19: "0.514653",
                  21: "0.463784", 23: "0.417731"}

MBAR_ORIGIN = "smallest multiple of 1e-7 above the floating mean-field root passing the exact test"
PART1 = "choose_rows.py (procedure of the first implementation)"
# label, d, t, p, g_K, g_h, criterion, published Table 2 bounds, origin of (p, g_K, g_h).
ROWS = [
    ("d12-t3/25", 12, "3/25", "2161/10000", "102634/100000", "83527/100000", "basic",
     ("0.4351284", "0.168984", "0.907068", "0.897470", "0.036292"),
     "p from the first implementation's procedure of choose_rows.py; Holley line = exploratory line K'=0.013, "
     "h=-0.09 with e^(2K'), e^(2h) rounded down to five decimals (Propositions 4.7 and 4.12)"),
    ("d12-t11/100", 12, "11/100", "2023/10000", "255577/250000", "873143/1000000", "basic",
     ("0.4540507", "0.127091", "0.926888", "0.918156", "0.031296"), PART1),
    ("d12-t23/200", 12, "23/200", "2094/10000", "1024283/1000000", "855173/1000000", "basic",
     ("0.4453800", "0.146071", "0.913716", "0.904575", "0.035441"), PART1),
    ("d12-t1/8", 12, "1/8", "2224/10000", "1028469/1000000", "407843/500000", "basic",
     ("0.4240888", "0.196297", "0.905149", "0.895074", "0.034915"), PART1),
    ("d12-t13/100", 12, "13/100", "2283/10000", "515339/500000", "794291/1000000", "basic",
     ("0.4111645", "0.230516", "0.910357", "0.899716", "0.030630"), PART1),
    ("d11-t13/100", 11, "13/100", "2310/10000", "6443/6250", "25553/31250", "basic",
     ("0.4252374", "0.215864", "0.952104", "0.940625", "0.017392"), PART1),
    ("d11-t3/25", 11, "3/25", "2180/10000", "513251/500000", "427717/500000", "basic",
     ("0.4454995", "0.162436", "0.960915", "0.950461", "0.015882"), PART1),
    ("d10-t27/200", 10, "27/200", "2399/10000", "516713/500000", "102809/125000", "refined_c",
     ("0.4279573", "0.236224", "1.007739", "0.994576", "0.001950"),
     "choose_rows.py (procedure of the second implementation, d = 10, t = 27/200)"),
    ("d10-t13/100", 10, "13/100", "2336/10000", "1031101/1000000", "840843/1000000", "refined_c",
     ("0.4379919", "0.205637", "1.009832", "0.997261", "0.001048"), PART1),
]
PUBLISHED_KEYS = ("rho_minus_lower", "eta_prime_upper", "score_upper", "score_c_upper",
                  "theta_star_lower")
# The inputs p, g_K, g_h as the terminating decimals of Table 2 of the manuscript.
TABLE_INPUTS = {
    "d12-t3/25": ("0.2161", "1.02634", "0.83527"),
    "d12-t11/100": ("0.2023", "1.022308", "0.873143"),
    "d12-t23/200": ("0.2094", "1.024283", "0.855173"),
    "d12-t1/8": ("0.2224", "1.028469", "0.815686"),
    "d12-t13/100": ("0.2283", "1.030678", "0.794291"),
    "d11-t13/100": ("0.2310", "1.03088", "0.817696"),
    "d11-t3/25": ("0.2180", "1.026502", "0.855434"),
    "d10-t27/200": ("0.2399", "1.033426", "0.822472"),
    "d10-t13/100": ("0.2336", "1.031101", "0.840843"),
}
TABLE_INPUT_KEYS = ("p_exact", "g_K_exact", "g_h_exact")
# Further displays of the manuscript's Section 4 (Sections 4.2 and 4.7, Remarks 4.6 and 4.40). 'margin' is
# 1 - Score (dimensions 11, 12) or 1 - Score_c (dimension 10), displayed as "about" a percentage: the
# certified margins range from about 7.3% to 9.5% in dimension 12 (about 9% at t = 3/25), are about 4.8%
# and 3.9% in dimension 11, and 0.54% and 0.27% in dimension 10; 'rho_c_gain' = rho_c - rho_- is about
# 0.006 in dimension 10. 'sg_kappa' is the parameter kappa_24 ~ 5.5 of the centered susceptibility bound.
# Every row: p_B is "about 0.9 tanh 2beta", checked as p / tanh 2beta >= 0.85 and witness / tanh 2beta <= 0.95.
COMMON_PUBLISHED = {"p_over_tanh_2beta_lower": "0.85", "pB_witness_over_tanh_2beta_upper": "0.95"}
EXTRA_PUBLISHED = {
    "d12-t3/25": {
        "sg_kappa_lower": "5.45", "sg_kappa_upper": "5.55",
        "kappa_upper": "1.29818", "alpha_upper": "0.311972",
        "term1_upper": "0.886231", "term2_upper": "0.019603", "term3_upper": "0.001235",
        "pB_witness_upper": "0.2210000",
        "p_A_exact": "75/196", "tanh_2beta_exact": "75/317", "closed_floor_exact": "75/392",
        "score_closed_floor_lower": "1.015", "score_closed_floor_upper": "1.025",
        "argmax_delta_exact": "20", "max_delta_lower": "-0.0001955", "max_delta_upper": "-0.0001945",
        "frozen_success_1_upper": "0.2186899", "frozen_success_2_upper": "0.2163712",
        "frozen_success_3_upper": "0.2140597", "frozen_success_4_upper": "0.2117739",
        "holley_ratio_S_max_loglower": "0.001827", "holley_ratio_S_min_loglower": "0.00192",
        "holley_ratio_interior_min_loglower": "0.0057",
        "K_prime_lo_lower": "0.0129995", "K_prime_hi_upper": "0.0130005",
        "abs_h_lo_lower": "0.0899995", "abs_h_hi_upper": "0.0900005",
        "gamma_min_lower": "0.305", "gamma_min_upper": "0.315",
        "G_d_lo_lower": "1.1011537095", "G_d_hi_upper": "1.1011552715",
        "F_d_minus_1_over_d_upper": "0.0085295258", "G2_d_upper": "1.2272959943",
        "margin_lower": "0.085", "margin_upper": "0.095",
    },
    "d12-t11/100": {"term3_upper": "0.0019", "margin_lower": "0.0725", "margin_upper": "0.0735"},
    "d12-t23/200": {"term3_upper": "0.0019", "margin_lower": "0.0725", "margin_upper": "0.0955"},
    "d12-t1/8": {"term3_upper": "0.0019", "margin_lower": "0.0945", "margin_upper": "0.0955"},
    "d12-t13/100": {"term3_upper": "0.0019", "margin_lower": "0.0725", "margin_upper": "0.0955"},
    "d11-t13/100": {"F_d_minus_1_over_d_upper": "0.0104053880", "margin_lower": "0.0475", "margin_upper": "0.0485"},
    "d11-t3/25": {"margin_lower": "0.0385", "margin_upper": "0.0395"},
    "d10-t27/200": {"F_d_minus_1_over_d_upper": "0.0130075189", "term3_upper": "0.0034", "score_lower": "1.0077",
                    "margin_lower": "0.00535", "margin_upper": "0.00545",
                    "rho_c_gain_lower": "0.0055", "rho_c_gain_upper": "0.0065"},
    "d10-t13/100": {"term3_upper": "0.0034", "score_lower": "1.0098",
                    "margin_lower": "0.00265", "margin_upper": "0.00275",
                    "rho_c_gain_lower": "0.0055", "rho_c_gain_upper": "0.0065"},
}
# Tangent weights W_K of the first implementation (crosscheck/certify_local.py); the manuscript's Holley
# slack at d = 12, t = 3/25 is displayed for these weights.
HEADLINE_WEIGHTS = ["0"] * 16 + ["0.0442", "0.0964", "0.1519", "0.2104", "0.2718", "0.3359",
                                 "0.4027", "0.4719", "0.5435"]
# Rows that the manuscript reports as not certified: label, d, t, p, g_K, g_h, displays.
UNCERTIFIED = [("d10-t3/25", 10, "3/25", "2199/10000", "1026679/1000000", "437233/500000",
                {"score_c_lower": "1.0192", "score_c_upper": "1.0193"})]


class Model(exact.Model):
    """The exact model of ../certify.py plus a floating beta for the searches."""

    def __init__(self, d, t):
        super().__init__(d, t)
        self.beta = atanh(float(self.t))


# ---------------------------------------------------------------- p_B tangent family
def float_max_delta(M, p, kappa, lx):
    n = M.m - 1; b = M.beta; a = float(M.a); w = float(M.w)
    r = float(p / M.pA); lam = w * w * (1 - r) / r
    K0 = 2 * cosh(b) ** M.m
    best = -1e9
    for j in range(M.m):
        tot = 0.0
        for al in range(j + 1):
            for bb in range(n - j + 1):
                A = 2 * al - j; B = 2 * bb - (n - j); U = A + B; z = A - B
                P = comb(j, al) * comb(n - j, bb) * a ** (al + bb) * (1 - a) ** (n - al - bb)
                kap = w ** (-2 * U) - lam
                om = 2 * cosh(b * (1 + z))
                if kap >= 0:
                    phi = om ** -2
                else:
                    c = kappa * (2 * cosh(b * (1 + U))) ** lx * K0 ** (1 - lx)
                    phi = 3 / c ** 2 - 2 * om / c ** 3
                tot += P * kap * phi
        best = max(best, tot)
    return best


def search_family(M, p):
    f = lambda x: float_max_delta(M, p, x[0], x[1])
    res = min((minimize(f, [k0, l0], method='Nelder-Mead',
                        options={'xatol': 1e-9, 'fatol': 1e-15, 'maxiter': 3000})
               for k0 in (0.98, 1.02) for l0 in (0.5, 0.6)), key=lambda r: r.fun)
    return float(res.x[0]), float(res.x[1])


def tangent_family(M, p, kappa, lx, digits=10):
    r = p / M.pA; lam = M.w ** 2 * (1 - r) / r
    K0 = 2 * cosh(M.beta) ** M.m
    cU = {}
    for U in range(-(M.m - 1), M.m, 2):
        if M.w ** (-2 * U) < lam:
            v = kappa * (2 * cosh(M.beta * (1 + U))) ** lx * K0 ** (1 - lx)
            cU[U] = Fr(round(v * 10 ** digits), 10 ** digits)
    return cU


# ---------------------------------------------------------------- Holley tangent weights
def float_numerator(M, k, wk):
    b = M.beta; a = float(M.a); km = M.m - k
    Gf = lambda z: 2 * cosh(b) ** km * cosh(b * z)
    if k == 0:
        return Gf(0) ** -2
    c = [(1 - wk) * Gf(2 * i - k) + wk * (i / k * Gf(2 * i - k - 2) + (k - i) / k * Gf(2 * i - k + 2))
         for i in range(k + 1)]
    vals = []
    for j in range(k + 1):
        tot = 0.0
        for al in range(j + 1):
            for bb in range(k - j + 1):
                A = 2 * al - j; B = 2 * bb - (k - j); ci = c[al + bb]
                tot += (comb(j, al) * comb(k - j, bb) * a ** (al + bb) * (1 - a) ** (k - al - bb)
                        * (3 / ci ** 2 - 2 * Gf(A - B) / ci ** 3))
        vals.append(tot)
    return min(vals)


def choose_wk(M, k):
    if k < 2:
        return Fr(0)
    best = (float_numerator(M, k, 0.0), 0.0)
    for i in range(1, 181):
        x = i / 200
        v = float_numerator(M, k, x)
        if v > best[0]:
            best = (v, x)
    return Fr(round(best[1] * 10 ** 4), 10 ** 4)


# ---------------------------------------------------------------- mean-field mbar
def mf_ok(d, gK, gh, xbar):
    """The exact root test of ../certify.py, as a boolean."""
    H_hi = exact.atanh_bounds((1 - gh) / (1 + gh))[1]
    K_hi = exact.atanh_bounds((gK - 1) / (gK + 1))[1]
    return H_hi + 2 * d * K_hi * xbar < exact.atanh_bounds(xbar)[0]


def find_xbar(d, gK, gh):
    Kp = log(float(gK)) / 2; H = -log(float(gh)) / 2; x = 0.0
    for _ in range(5000):
        x = tanh(H + 2 * d * Kp * x)
    xb = Fr(int(x * 10 ** 7) + 1, 10 ** 7)
    while not mf_ok(d, gK, gh, xb):
        xb += Fr(1, 10 ** 7)
    return xb


def decimal(x, digits):
    """Exact shortest decimal string of a rational with at most `digits` decimal places."""
    scaled = Fr(x) * 10 ** digits
    assert scaled.denominator == 1, x
    q = scaled.numerator
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}".rstrip("0").rstrip(".")


def published_block(label, table_bounds):
    """The displays of the manuscript that ../certify.py asserts for one row."""
    return {**dict(zip(TABLE_INPUT_KEYS, TABLE_INPUTS[label])), **dict(zip(PUBLISHED_KEYS, table_bounds)),
            **COMMON_PUBLISHED, **EXTRA_PUBLISHED.get(label, {})}


def build():
    rows = []
    for label, d, t, p, gK, gh, criterion, published, origin in ROWS:
        M = Model(d, Fr(t))
        entry = {"label": label, "d": d, "t": t, "p": p, "g_K": gK, "g_h": gh}
        if label == "d12-t3/25":
            tangent = {"form": "rho_times_omega", "values": {str(U): v for U, v in HEADLINE_RHO_U.items()}}
            tangent_origin = ("first-implementation values rho_U from the family kappa=1.01439, "
                              "lambda=0.56291, rounded to 6 decimals (Proposition 4.7)")
        else:
            kappa, lam = search_family(M, Fr(p))
            cU = tangent_family(M, Fr(p), kappa, lam)
            tangent = {"form": "direct", "values": {str(U): decimal(c, 10) for U, c in sorted(cU.items())}}
            tangent_origin = f"Nelder-Mead family kappa={kappa!r}, lambda={lam!r}; c_U rounded to 10 decimals"
            print(f"{label}: tangent family kappa={kappa:.6f} lambda={lam:.6f}", flush=True)
        if label == "d12-t3/25":
            weights = HEADLINE_WEIGHTS
            weights_origin = ("weights W_K of the first implementation, crosscheck/certify_local.py (bounded "
                              "floating maximisation of the numerator bound, rounded to 4 decimals; "
                              "Proposition 4.12); they give the Holley slack displayed in the manuscript")
        else:
            weights = [decimal(choose_wk(M, k), 4) for k in range(M.m + 1)]
            weights_origin = "best point of {0} U {i/200: 1<=i<=180} for the floating numerator bound"
        mbar = find_xbar(d, Fr(gK), Fr(gh))
        entry.update({
            "mbar": f"{mbar.numerator}/{mbar.denominator}",
            "criterion": criterion,
            "pB_tangent_points": tangent,
            "holley_weights": weights,
            "published": published_block(label, published),
            "provenance": {
                "row": origin,
                "pB_tangent_points": tangent_origin,
                "holley_weights": weights_origin,
                "mbar": MBAR_ORIGIN,
            },
        })
        print(f"{label}: mbar={entry['mbar']}; w_k={weights}", flush=True)
        rows.append(entry)
    uncertified = []
    for label, d, t, p, gK, gh, published in UNCERTIFIED:
        mbar = find_xbar(d, Fr(gK), Fr(gh))
        uncertified.append({
            "label": label, "d": d, "t": t, "p": p, "g_K": gK, "g_h": gh,
            "mbar": f"{mbar.numerator}/{mbar.denominator}", "published": published,
            "provenance": {"row": ("choose_rows.py (procedure of the first implementation; reported as not "
                                   "certified in Section 4.7)"),
                           "mbar": MBAR_ORIGIN},
        })
    return {
        "description": ("Frozen exact inputs of the overlap-revealed certificates. Numbers are exact decimal "
                        "or fraction strings, read by certify.py with fractions.Fraction. Each row's "
                        "'published' block lists numbers displayed in the manuscript that certify.py must "
                        "reproduce: the key is '<quantity>_<direction>', where the direction says whether the "
                        "display is an upper bound, a lower bound, a lower bound on the logarithm ('loglower') "
                        "or an exact value. 'uncertified_rows' lists rows at which only the global "
                        "criterion is evaluated, to confirm that it fails there. 'provenance' is documentation "
                        "only; search/freeze_params.py re-derives the searched values."),
        "tangent_forms": {
            "rho_times_omega": "c_U = value * omega(U), omega(U) = w^((1+U)/2) + w^(-(1+U)/2), w = (1+t)/(1-t)",
            "direct": "c_U = value",
        },
        "rows": rows,
        "uncertified_rows": uncertified,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", type=Path, help="write the recomputed parameters to this path")
    args = parser.parse_args()
    fresh = build()
    if args.write:
        args.write.write_text(json.dumps(fresh, indent=2) + "\n")
        print(f"wrote {args.write}")
        return
    frozen = json.loads(PARAMS.read_text())
    if fresh != frozen:
        keys = sorted({k for a in fresh["rows"] for k in a})
        diffs = [(a["label"], k) for a, b in zip(fresh["rows"], frozen["rows"]) for k in keys if a.get(k) != b.get(k)]
        print("Different from params.json (still admissible, but not the frozen values):", diffs or "top level")
        sys.exit(1)
    print("PASS: the search reproduces params.json exactly")


if __name__ == "__main__":
    main()

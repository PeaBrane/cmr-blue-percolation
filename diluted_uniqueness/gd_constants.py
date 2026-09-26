"""Exact constants of the uniqueness window for the diluted +-1 EA model.

With surgery radius rho = 4 (largest surgery box B_5, side 11), the
diminishment argument uses

    K_F = 12 d 11^(d-1),  K_Q = C(d,2) 14^2 13^(d-2),  N_1 = d (4 * 11^d + 9^d),
    g_d = (p_c/2)^K_F / (N_1 2^(K_F + K_Q)),

and uniqueness at every temperature holds for p < p_c + g_d/2. This program
checks, with exact integers and rationals:

1. For d = 2, 3, 4, K_F equals the number of edges meeting B_5 and K_Q the
   number of plaquettes with a corner in B_6, both counted by enumeration, as
   in the manuscript's lemma on the sizes of the surgery sets (N_1 is a union
   bound there, not a count). The displayed values hold:
   (K_F, K_Q, N_1, 3K_F+K_Q) = (264, 196, 1130, 988) for d = 2 and
   (4356, 7644, 18159, 20712) for d = 3.
2. Since p_c >= 1/(2d-1), g_d >= 1/den with den = (4d-2)^K_F 2^(K_F+K_Q) N_1,
   and 10^(k-1) <= den < 10^k with k = 347, 7973, 152290 for d = 2, 3, 4
   (for d = 3, den = 10^K_F 2^(K_F+K_Q) N_1 has 7973 digits).
3. For d = 2 with p_c(Z^2) = 1/2: g_2 = 4^-264/(1130 * 2^460) = 2^-988/1130,
   the window width is g_2/2 = 2^-989/1130, log10(1/g_2) = 300.4707..., and
   3.38e-301 < g_2 < 3.39e-301.
4. The classical window bound beta_*(d) = (1/2)[(3K_F+K_Q) ln 2 + ln N_1]
   satisfies beta_*(2) in [345.929693652, 345.929693653] and beta_*(3) in
   [7183.13566267, 7183.13566268], so 1/beta_*(2) < 2.8908e-3 (hence
   < 2.9e-3) and 1/beta_*(3) < 1.3922e-4. Logarithms are enclosed by the
   atanh series ln x = 2 atanh((x-1)/(x+1)) with a geometric tail bound.

Assertions decide every comparison; only the standard library is used.
"""

from fractions import Fraction as Fr
from itertools import combinations, product
from math import comb
import json
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

SERIES_TERMS = 60
DIGITS_OF_DEN = {2: 347, 3: 7973, 4: 152290}
DISPLAYED = {2: (264, 196, 1130, 988), 3: (4356, 7644, 18159, 20712)}
BETA_STAR = {2: ("345.929693652", "345.929693653"), 3: ("7183.13566267", "7183.13566268")}
INVERSE_BETA_STAR = {2: "2.8908e-3", 3: "1.3922e-4"}


def constants(d):
    KF = 12 * d * 11 ** (d - 1)
    KQ = comb(d, 2) * 14 ** 2 * 13 ** (d - 2)
    N1 = d * (4 * 11 ** d + 9 ** d)
    return KF, KQ, N1


def recount(d, r=5):
    """Enumerate the edges meeting the box B_r and the plaquettes with a corner in B_(r+1).

    These are the two counts of the lemma on the sizes of the surgery sets, at the largest
    surgery box B_5 (side 11)."""
    edges = set()
    for x in product(range(-r, r + 1), repeat=d):
        for i in range(d):
            for shift in (-1, 0):                    # the two edges at x in direction i
                low = list(x)
                low[i] += shift
                edges.add((tuple(low), i))
    R = r + 1
    plaquettes = 0
    for x in product(range(-R - 1, R + 1), repeat=d):  # lowest corner of the plaquette
        for i, j in combinations(range(d), 2):
            for a, b in ((0, 0), (1, 0), (0, 1), (1, 1)):
                corner = list(x)
                corner[i] += a
                corner[j] += b
                if max(map(abs, corner)) <= R:
                    plaquettes += 1
                    break
    return len(edges), plaquettes


def atanh_bounds(z):
    """lo <= atanh(z) <= hi for rational 0 <= z < 1."""
    z = Fr(z)
    assert 0 <= z < 1
    lo = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(SERIES_TERMS))
    return lo, lo + z ** (2 * SERIES_TERMS + 1) / ((2 * SERIES_TERMS + 1) * (1 - z * z))


def log_bounds(n):
    """Enclosure of ln n for an integer n >= 2: e ln 2 + ln(n / 2^e) with 2^e <= n < 2^(e+1)."""
    e = n.bit_length() - 1
    two_lo, two_hi = (2 * x for x in atanh_bounds(Fr(1, 3)))
    rest_lo, rest_hi = (2 * x for x in atanh_bounds(Fr(n - 2 ** e, n + 2 ** e)))
    return e * two_lo + rest_lo, e * two_hi + rest_hi


def dec(x, digits, up):
    """x rounded down (up=False) or up (up=True) to `digits` decimals, exactly."""
    q = x.numerator * 10 ** digits
    q = -((-q) // x.denominator) if up else q // x.denominator
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def main():
    report = {}
    for d, shown in DISPLAYED.items():
        KF, KQ, N1 = constants(d)
        assert (KF, KQ, N1, 3 * KF + KQ) == shown
    for d, k in DIGITS_OF_DEN.items():
        KF, KQ, N1 = constants(d)
        assert recount(d) == (KF, KQ)
        den = (4 * d - 2) ** KF * 2 ** (KF + KQ) * N1
        assert 10 ** (k - 1) <= den < 10 ** k
        report[f"d{d}"] = {"K_F": KF, "K_Q": KQ, "N_1": N1, "g_d_lower": f"10^-{k}"}
        print(f"PASS d={d}: K_F={KF}, K_Q={KQ} (both recounted), N_1={N1}; g_d >= 1/den > 10^-{k}", flush=True)

    KF, KQ, N1 = constants(2)
    inv_g2 = 2 ** (3 * KF + KQ) * N1
    g2 = Fr(1, 4) ** KF / (N1 * 2 ** (KF + KQ))              # (p_c/2)^K_F / (N_1 2^(K_F+K_Q)), p_c = 1/2
    assert g2 == Fr(1, 4 ** 264 * 1130 * 2 ** 460) == Fr(1, 2 ** 988 * 1130) == Fr(1, inv_g2)
    assert g2 / 2 == Fr(1, 2 ** 989 * 1130)
    assert 10 ** 300 < inv_g2 < 10 ** 301
    assert 338 * inv_g2 < 10 ** 303 < 339 * inv_g2          # 3.38e-301 < g_2 < 3.39e-301
    ln_inv_lo, ln_inv_hi = log_bounds(inv_g2)
    ln10_lo, ln10_hi = log_bounds(10)
    assert Fr("300.4707") <= ln_inv_lo / ln10_hi <= ln_inv_hi / ln10_lo < Fr("300.4708")
    report["d2"]["g_2_exact"] = "2^-988/1130"
    report["d2"]["g_2_bracket"] = ["3.38e-301", "3.39e-301"]
    print("PASS d=2 with p_c = 1/2: g_2 = 2^-988/1130, 3.38e-301 < g_2 < 3.39e-301", flush=True)

    ln2_lo, ln2_hi = (2 * x for x in atanh_bounds(Fr(1, 3)))
    for d, (shown_lo, shown_hi) in BETA_STAR.items():
        KF, KQ, N1 = constants(d)
        lnN_lo, lnN_hi = log_bounds(N1)
        lo = ((3 * KF + KQ) * ln2_lo + lnN_lo) / 2
        hi = ((3 * KF + KQ) * ln2_hi + lnN_hi) / 2
        assert Fr(shown_lo) <= lo <= hi <= Fr(shown_hi)
        report[f"d{d}"]["beta_star"] = [dec(lo, 15, False), dec(hi, 15, True)]
        report[f"d{d}"]["inverse_beta_star_upper"] = dec(1 / lo, 12, True)
        assert 1 / lo < Fr(INVERSE_BETA_STAR[d])
        print(f"PASS beta_*({d}) in [{shown_lo}, {shown_hi}], 1/beta_*({d}) < {INVERSE_BETA_STAR[d]}", flush=True)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

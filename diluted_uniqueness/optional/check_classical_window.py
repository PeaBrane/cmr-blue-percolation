"""Interval-arithmetic cross-check of the window constants (optional; needs mpmath).

The decisive version is ../gd_constants.py (standard library only). Manuscript locators: the constant
(1.4), the threshold (1.7) and Appendix A.1 of the diluted-model manuscript.

Closed-form quantities only (4 expressions per dimension); no case-by-case certificate.

(1) g_d = (p_c/2)^K_F / (N_1 2^(K_F+K_Q)), eq. (1.4), with
      K_F = 12 d 11^(d-1),  K_Q = C(d,2) 14^2 13^(d-2),  N_1 = d (4*11^d + 9^d).
    d = 2: p_c(Z^2) = 1/2 (Kesten 1980, Thm 1), so 1/g_2 = 2^(3 K_F + K_Q) N_1 is an integer.
    d = 3: 1/5 <= p_c(Z^3) <= 1/2 (Grimmett 1999 (1.13) with lambda(3) <= 5; (1.9) and Kesten),
           so g_3 >= 1/den3 with den3 = 10^K_F 2^(K_F+K_Q) N_1 (exact integer).
(2) Classical window (Newman 1997, Thms 3.10-3.11; ACCN 1987): for mu_p the EA model is unique at
    every beta < beta_l(p) = (1/2) ln(p/(p-p_c)).  For p in (p_c, p_c + g_d/2),
      beta_l(p) > (1/2) ln(2 p_c / g_d) = (1/2)[ln 2 + ln p_c + K_F ln(2/p_c) + (K_F+K_Q) ln 2 + ln N_1],
    which is decreasing in p_c (K_F > 1); at the admissible maximum p_c = 1/2 it equals
      beta_*(d) := (1/2)[(3 K_F + K_Q) ln 2 + ln N_1].
    So beta_l(p) > beta_*(d) for every p in the window, using only p_c <= 1/2.
All logarithms: mpmath interval arithmetic (outward rounding), 50 digits. The program asserts that
the enclosures lie inside the manuscript's displays: log10(1/g_2) = 300.4707..., g_2 in
(3.38e-301, 3.39e-301), 7973 digits for d = 3, the beta_*(d) intervals and the bounds on 1/beta_*(d).
"""
from fractions import Fraction
from math import ceil, comb, floor
import sys

import mpmath
from mpmath import iv

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(10 ** 6)
iv.dps = 50


def consts(d):
    KF = 12 * d * 11 ** (d - 1)
    KQ = comb(d, 2) * 14 ** 2 * 13 ** (d - 2)
    N1 = d * (4 * 11 ** d + 9 ** d)
    return KF, KQ, N1


def lo(x):
    return mpmath.mpf(x.a)


def hi(x):
    return mpmath.mpf(x.b)


def exact(x):
    """The binary mpf x as an exact Fraction."""
    man, exp = mpmath.mpf(x).man_exp
    return Fraction(man) * Fraction(2) ** exp


def directed(x, digits, up):
    """Decimal string of the mpf x rounded down (up=False) or up (up=True), exactly."""
    q = exact(x) * 10 ** digits
    q = ceil(q) if up else floor(q)
    assert q >= 0
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


# Displays of the manuscript (appendix on explicit constants), checked against the enclosures below.
BETA_STAR = {2: ("345.929693652", "345.929693653", "2.8908e-3"), 3: ("7183.13566267", "7183.13566268", "1.3922e-4")}


ln2 = iv.log(iv.mpf(2))
ln10 = iv.log(iv.mpf(10))

for d in (2, 3):
    KF, KQ, N1 = consts(d)
    # independent recount of K_F, K_Q from the box geometry (r = 5, n = 11)
    assert KF == d * 11 ** (d - 1) * 12
    assert KQ == comb(d, 2) * (2 * 5 + 4) ** 2 * (2 * 5 + 3) ** (d - 2)
    print(f'd={d}: K_F={KF}  K_Q={KQ}  N_1={N1}')
    if d == 2:
        inv_g = 2 ** (3 * KF + KQ) * N1                     # exact, p_c = 1/2
        l10 = (iv.mpf(3 * KF + KQ) * ln2 + iv.log(iv.mpf(N1))) / ln10
        print(f'   1/g_2 = 2^{3 * KF + KQ} * {N1};  log10(1/g_2) in [{directed(lo(l10), 9, False)}, {directed(hi(l10), 9, True)}]')
        assert inv_g < 10 ** 301 and inv_g > 10 ** 300
        assert Fraction("300.4707") <= exact(lo(l10)) and exact(hi(l10)) < Fraction("300.4708")
        # mantissa: g_2 = 10^-301 * m with m = 10^301 / inv_g
        m = iv.mpf(10) ** 301 / iv.mpf(inv_g)
        print(f'   g_2 = m * 1e-301 with m in [{directed(lo(m), 7, False)}, {directed(hi(m), 7, True)}]  (so g_2 > 3e-301)')
        assert lo(m) > 3 and Fraction("3.38") < exact(lo(m)) and exact(hi(m)) < Fraction("3.39")
    else:
        den = 10 ** KF * 2 ** (KF + KQ) * N1               # g_3 >= 1/den (p_c >= 1/5)
        k = len(str(den))
        assert k == 7973
        print(f'   g_3 >= 1/den, den has {k} digits, so g_3 > 10^-{k}')
    beta_star = (iv.mpf(3 * KF + KQ) * ln2 + iv.log(iv.mpf(N1))) / 2
    T_star = 1 / beta_star
    print(f'   beta_*({d}) in [{directed(lo(beta_star), 9, False)}, {directed(hi(beta_star), 9, True)}];'
          f'  1/beta_* < {directed(hi(T_star), 10, True)}')
    shown_lo, shown_hi, shown_inv = BETA_STAR[d]
    assert Fraction(shown_lo) <= exact(lo(beta_star)) and exact(hi(beta_star)) <= Fraction(shown_hi)
    assert exact(hi(T_star)) < Fraction(shown_inv)
    # monotonicity in p_c of ln(2 p_c/g_d): derivative (1 - K_F)/p_c < 0
    assert KF > 1

# beta_u - beta_l = (1/2) ln((p + p_c)/p) in (0, (1/2) ln 2) for p > p_c: closed form, no numerics needed.
print('beta_u(p) - beta_l(p) = (1/2) ln(1 + p_c/p) < (1/2) ln 2 for p > p_c')

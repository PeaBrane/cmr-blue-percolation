"""Lemma 4.28 of the manuscript (c'_cov, rho'_c) and the least torus size L_* of its variant, in exact rational arithmetic.

This is the second of the two implementations of Lemma 4.28 (Section 4.7), written from the statement of the lemma
without code from ../crosscheck/covariance.py. Standard library only.

    c'_cov = 2 sinh(2K') / ( e^{K'} (1 + E_+ c(U_+)) + e^{-K'} (1 + E_- c(U_-)) )^2,   c(U) = (cosh U - 1)/U^2,
    U_+ = 2(Hb + (2d-1)Kb) >= 2M,   U_- = 2(2d-1)Kb >= 2M',
    E_+ = 4 Hb^2 + 4 Hb Kb (4d-2) mbar + Kb^2 E_T,   E_- = 2 Kb^2 E_1,
    E_1 = (2d-1) + (2d-1)(2d-2)(mbar^2 + delta_2),
    E_T = 2 E_1 + 2[(2d-2)(mbar^2 + delta_1) + ((2d-1)^2 - (2d-2))(mbar^2 + delta_3)],
    delta_r = N^max_r tbar^r Gamma_r(tbar) + 10^-15,   N^max = (1, 2, 6),
    Gamma_r(t) = 2^-d sum_j C(d,j) (1 - 2t(d-2j))^-(r+1).

Evaluation choices, each in the safe direction (all quantities are monotone in the parameters in the direction
used): 2 sinh(2K') = g_K - 1/g_K and e^{+-K'} = g_K^{+-1/2} exactly, with rational bounds on the square roots;
Kb >= K' and Hb >= |h| from the atanh series with a geometric remainder; tbar = tanh K' = (g_K-1)/(g_K+1) exactly;
mbar = the frozen xbar, which passes the exact root test; c(U) bounded by its power series with a geometric
remainder.

Lemma 4.28 needs 2d Kb <= 35/100 (and L >= 40). Its variant replaces this by L >= L_*, where L_* >= 4 is an integer
with abar^(L_*-3)/(1-abar) <= 10^-15 for a rational abar >= 2d tanh K', abar < 1; least_Lstar returns the least one.
"""
from fractions import Fraction as Fr
from math import comb, factorial, isqrt

EPS_WRAP = Fr(1, 10 ** 15)


def atanh_hi(z, K=40):
    z = Fr(z)
    assert 0 <= z < 1
    s = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(K + 1))
    return s + z ** (2 * K + 3) / ((2 * K + 3) * (1 - z * z))


def sqrt_hi(x, digits=40):
    x = Fr(x)
    S = 10 ** digits
    r = Fr(isqrt(x.numerator * S * S // x.denominator) + 1, S)
    assert r * r >= x
    return r


def sqrt_lo(x, digits=40):
    x = Fr(x)
    S = 10 ** digits
    r = Fr(isqrt(x.numerator * S * S // x.denominator), S)
    assert r * r <= x
    return r


def c_hi(U, K=30):
    """Rational >= c(U) = (cosh U - 1)/U^2 = sum_{k>=1} U^{2k-2}/(2k)!."""
    U = Fr(U)
    U2 = U * U
    s = sum(U2 ** (k - 1) / factorial(2 * k) for k in range(1, K + 1))
    first = U2 ** K / factorial(2 * K + 2)
    ratio = U2 / ((2 * K + 3) * (2 * K + 4))
    assert ratio < 1
    return s + first / (1 - ratio)


def gamma_r(d, t, r):
    return sum(comb(d, j) / (1 - 2 * t * (d - 2 * j)) ** (r + 1) for j in range(d + 1)) / 2 ** d


def least_Lstar(abar):
    abar = Fr(abar)
    assert 0 < abar < 1
    L = 4
    while abar ** (L - 3) / (1 - abar) > EPS_WRAP:
        L += 1
    return L


def lemma36(d, gK, gh, mbar):
    """Lower bound on c'_cov (Lemma 4.28) and the auxiliary quantities."""
    gK, gh, mbar = Fr(gK), Fr(gh), Fr(mbar)
    n1 = 2 * d - 1
    Kb = atanh_hi((gK - 1) / (gK + 1))
    Hb = atanh_hi((1 - gh) / (1 + gh))
    tbar = (gK - 1) / (gK + 1)
    assert 2 * d * tbar < 1
    delta = {r: N * tbar ** r * gamma_r(d, tbar, r) + EPS_WRAP for r, N in ((1, 1), (2, 2), (3, 6))}
    m2 = mbar * mbar
    E1 = n1 + n1 * (n1 - 1) * (m2 + delta[2])
    ET = 2 * E1 + 2 * ((2 * d - 2) * (m2 + delta[1]) + (n1 * n1 - (2 * d - 2)) * (m2 + delta[3]))
    Eplus = 4 * Hb ** 2 + 4 * Hb * Kb * (4 * d - 2) * mbar + Kb ** 2 * ET
    Eminus = 2 * Kb ** 2 * E1
    Up, Um = 2 * (Hb + n1 * Kb), 2 * n1 * Kb
    den = sqrt_hi(gK) * (1 + Eplus * c_hi(Up)) + (1 + Eminus * c_hi(Um)) / sqrt_lo(gK)
    ccov = (gK - 1 / gK) / den ** 2
    return ccov, dict(Kb=Kb, Hb=Hb, tbar=tbar, Eplus=Eplus, Eminus=Eminus, abar=2 * d * Kb, alpha=2 * d * tbar)

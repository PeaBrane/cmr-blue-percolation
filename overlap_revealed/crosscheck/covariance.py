"""Sharper lower bound on Cov(sigma_x, sigma_y), x~y, for the torus Ising law (global.md Lemma 3.6).

Cov >= 2 sinh(2K') / (E[P+Q])^2  (Jensen), P+Q = e^{K'} cosh(a+b) + e^{-K'} cosh(a-b),
E cosh(u) <= 1 + E[u^2] (cosh U - 1)/U^2 for |u| <= U,
E[(a+b)^2] <= 4H^2 + 4 H K' (4d-2) mbar + K'^2 E[T^2],   E[(a-b)^2] <= 2 K'^2 E[T_x^2],
E[sigma_i sigma_j] <= mbar^2 + D_ij (Lemma 3.6),  D_ij <= G(i-j) + 1e-15 for L >= 40.
Inputs are rationals: Klo <= K' <= Khi, H <= Hhi (H = |h|), mbar >= m_L, tbar >= tanh K'.
"""
from fractions import Fraction as Fr
from criterion import Gamma, iv, ivq, q_lo, q_up
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

EPS_WRAP = Fr(1, 10 ** 15)


def cov_lower_sharp(d, Klo, Khi, Hhi, mbar, tbar):
    n1 = 2 * d - 1                       # |N(x)\{y}|
    assert 2 * d * Khi <= Fr(35, 100)       # alpha <= 0.35: wrap term 0.35^37/0.65 < 1e-15 for L >= 40
    D1 = tbar * Gamma(d, tbar, 1) + EPS_WRAP          # |i-j|_1 = 1
    D2 = 2 * tbar ** 2 * Gamma(d, tbar, 2) + EPS_WRAP  # |i-j|_1 = 2 (N_2 <= 2)
    D3 = 6 * tbar ** 3 * Gamma(d, tbar, 3) + EPS_WRAP  # |i-j|_1 = 3 (N_3 <= 6)
    m2 = mbar ** 2
    # E[T_x^2] over one group (pairs within N(x)\{y} are at distance 2)
    ETx2 = n1 + n1 * (n1 - 1) * (m2 + D2)
    # cross pairs: 2d-2 at distance 1, rest at distance 3 (ordered pairs counted twice below)
    cross = (2 * d - 2) * (m2 + D1) + (n1 * n1 - (2 * d - 2)) * (m2 + D3)
    ET2 = 2 * ETx2 + 2 * cross
    Eplus = 4 * Hhi ** 2 + 4 * Hhi * Khi * (4 * d - 2) * mbar + Khi ** 2 * ET2
    Eminus = 2 * Khi ** 2 * ETx2
    M = Hhi + n1 * Khi
    M2 = n1 * Khi
    ch = lambda x: (iv.exp(x) + iv.exp(-x)) / 2
    sh = lambda x: (iv.exp(x) - iv.exp(-x)) / 2
    cM = (ch(2 * ivq(M)) - 1) / (2 * ivq(M)) ** 2
    cM2 = (ch(2 * ivq(M2)) - 1) / (2 * ivq(M2)) ** 2
    EPQ = iv.exp(ivq(Khi)) * (1 + ivq(Eplus) * cM) + iv.exp(-ivq(Klo)) * (1 + ivq(Eminus) * cM2)
    c = q_lo(2 * sh(2 * ivq(Klo)) / EPQ ** 2)
    return c, dict(Eplus=Eplus, Eminus=Eminus, ET2=ET2)

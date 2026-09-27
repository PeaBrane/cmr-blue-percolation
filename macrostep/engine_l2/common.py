"""Rigorous helpers for the macrostep second-moment engine (the directed-rounding program).

Rounding conventions.  Exact inputs are fractions.Fraction.  Heavy sums are done in IEEE-754 binary64 with
UPWARD rounding emulated by math.nextafter after every round-to-nearest operation (x (+) y and x (*) y are within
half an ulp of the exact result, so nextafter(fl(.), +inf) is an upper bound).  Where a long nonnegative dot
product is evaluated in round-to-nearest (numba kernels), the a priori bound |fl(sum) - sum| <= gamma_k * sum
(Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., eq. (3.5)) is applied explicitly, with a safety
factor, by an inflation (1 + INFL) documented at each use.
"""
import sys
import math
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
import os
IG_DIR = os.environ.get("IG_DIR", os.path.dirname(os.path.abspath(__file__)))   # indep_global.py lives next to this file
if IG_DIR not in sys.path:
    sys.path.insert(0, IG_DIR)
import indep_global as IG  # noqa: E402  (local constants; used read-only)

INF = float("inf")


def fup(x):
    """float >= x (x a Fraction or int)."""
    x = Fr(x)
    f = float(x)
    if Fr(f) < x:
        f = math.nextafter(f, INF)
    return f


def fdn(x):
    x = Fr(x)
    f = float(x)
    if Fr(f) > x:
        f = math.nextafter(f, -INF)
    return f


def add_up(a, b):
    return math.nextafter(a + b, INF)


def mul_up(a, b):
    return math.nextafter(a * b, INF)


def sum_up(xs):
    s = 0.0
    for x in xs:
        s = math.nextafter(s + x, INF)
    return s


def div_up(a, b):
    return math.nextafter(a / b, INF)


def local_constants(d, p, gK, gh, xbar=None):
    """Exact rational local constants from indep_global.evaluate:
    rho_- (mean-field lower bound, Lemma 4.22), rho_c (Lemma 4.27), t = tanh K' (exact), kappa = (1-rho)/rho,
    alpha = 2 d t, q_* = d t/(1 - 2 d t).  Asserts the hypotheses (33) and the mean-field root test."""
    p, gK, gh = Fr(p), Fr(gK), Fr(gh)
    if xbar is None:
        xbar = IG.find_xbar(d, gK, gh)
    o = IG.evaluate(d, p, gK, gh, Fr(xbar))
    return dict(d=d, p=p, gK=gK, gh=gh, xbar=Fr(xbar), rho=o["rho"], rho_c=o["rho_c"], t=o["t"],
                kappa=o["kappa"], alpha=o["alpha"], qs=o["qs"], score_oriented=o["score"],
                score_c_oriented=o["score_c"])


def b_exact(d, t, m):
    """Lemma 4.26(iv): b(m) = t^m M_m Gamma_m(t) >= sup_{|u|_1 = m} G_t(u)."""
    return t ** m * IG.Mmax_multinomial(d, m) * IG.Gam(d, t, m)


# local points (exact inputs); the points used by the stored certificates are decided exactly by
# ../certify_macrostep.py
POINTS = {
    (9, "7/50"): dict(p=Fr(623, 2500), gK=Fr(1036137, 10 ** 6), gh=Fr(83033, 10 ** 5), xbar=Fr(676921, 5000000)),
    (9, "29/200"): dict(p=Fr(2553, 10000), gK=Fr(1038647, 10 ** 6), gh=Fr(812491, 10 ** 6), xbar=Fr(1556769, 10 ** 7)),
    (9, "3/20"): dict(p=Fr(2611, 10000), gK=Fr(1041247, 10 ** 6), gh=Fr(12403, 15625), xbar=Fr(178451, 10 ** 6)),
    (9, "27/200"): dict(p=Fr(2427, 10000), gK=Fr(516847, 500000), gh=Fr(847407, 10 ** 6), xbar=Fr(293003, 2500000)),
    (9, "13/100"): dict(p=Fr(59, 250), gK=Fr(51567, 50000), gh=Fr(172699, 200000), xbar=Fr(1011213, 10 ** 7)),
    (8, "3/20"): dict(p=Fr(2651, 10000), gK=Fr(520811, 500000), gh=Fr(206239, 250000), xbar=Fr(1413829, 10 ** 7)),
    (8, "31/200"): dict(p=Fr(271, 1000), gK=Fr(208871, 200000), gh=Fr(25247, 31250), xbar=Fr(806041, 5000000)),
    (8, "4/25"): dict(p=Fr(2767, 10000), gK=Fr(1047157, 10 ** 6), gh=Fr(395043, 500000), xbar=Fr(1832723, 10 ** 7)),
    (8, "7/50"): dict(p=Fr(1261, 5000), gK=Fr(32389, 31250), gh=Fr(856969, 10 ** 6), xbar=Fr(1075653, 10 ** 7)),
    # [macrostep route] noisy-cavity inputs of Section 4.4 (Lemma 4.15, pair certificate Lemma 4.19, symmetrised
    # Holley line Lemma 4.20); their local certificates are ../local_sharpened.json, decided by
    # ../certify_macrostep.py.
    (8, "nc-3/20"): dict(p=Fr(2707, 10000), gK=Fr(1039407, 10 ** 6), gh=Fr(860671, 10 ** 6), xbar=Fr(539947, 5000000)),
    (8, "nc-31/200"): dict(p=Fr(1389, 5000), gK=Fr(65119, 62500), gh=Fr(845583, 10 ** 6), xbar=Fr(154899, 1250000)),
    (7, "nc-3/20"): dict(p=Fr(273, 1000), gK=Fr(519989, 500000), gh=Fr(44449, 50000), xbar=Fr(808479, 10 ** 7)),
    (7, "nc-31/200"): dict(p=Fr(1401, 5000), gK=Fr(104253, 100000), gh=Fr(175269, 200000), xbar=Fr(463903, 5000000)),
}

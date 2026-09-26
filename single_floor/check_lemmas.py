"""Supplementary exact checks behind verify_single_floor.py (standard library only).

1. Lemma J at all 256 (d, t, k) triples of the seven corollary points: the
   cavity maximum T_s over every orientation index equals its balanced value.
   With it, the complete h = 0 row p(1), ..., p(2d) is evaluated directly and
   its minimum is p(1) = p_*, independently of the block cover of
   Proposition K.
2. Identity (K1) for every k, and the monotonicity facts of Lemma M (i)-(v),
   at (d, t) = (16, 13/125) and (17, 1/10).
3. The count enumeration of the signed residual (paper eq. (50)) agrees with
   brute-force enumeration over all eps, x in {-1, 1}^Delta for
   Delta = 4, 6, k = 1, 2, 3 and three (t, lambda, c), and with a separate
   k = 1 formula at d = 16 and d = 22.
4. At d = 22 the signed residuals M_(0,1), M_(0,2) equal, as exact
   rationals, the maxima computed by the dimension-22 certificate
   ../verify_signed_star.py, an independent implementation.
"""

from fractions import Fraction as Fr
from itertools import product
from math import comb
from pathlib import Path
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from verify_single_floor import (POINTS, Reduction, Star, cosh_2beta, h_ell,  # noqa: E402
                                 nonneg_floor, signed_residual)
import verify_signed_star  # noqa: E402


def nonneg_floor_all_j(S, k):
    """p(k) with the cavity maximum over every orientation index; asserts Lemma J."""
    D, t, C, v, u, eta = S.D, S.t, S.C, S.v, S.u, S.eta
    s = D - k
    N = t * (1 - t * t) ** (D - 1) * C ** (k - 1) / h_ell(S, k - 1) ** 2
    cs = Fr(1) if s % 2 == 0 else C
    lam = (1 + eta ** (2 * k) + 2 * cs * eta ** k) / S.alpha
    T_all = [S.e_psi(D, s, j) for j in range(s + 1)]
    assert max(T_all) == T_all[s // 2], ("Lemma J", D, S.t, k)
    B = max((u - lam * v) * S.e_psi(D - 1, k - 1, j, 1) + (1 - u) * S.e_psi(D - 1, k - 1, j, -1)
            for j in range(k))
    M = C ** s * T_all[s // 2] / 2 + C ** k / 2 * B
    assert M >= 0
    return N / (lam * N + M)


def check_identity_k1(d, t):
    R = Reduction(d, t)
    S, D, C, eta = R.S, R.D, R.C, R.eta
    E = {}
    for k in range(1, D + 1):
        s = D - k
        E[k] = S.e_psi(D, s, s // 2)
        kap = R.kap(k)
        B_hat = max(eta * S.e_psi(D - 1, k - 1, j, -1) - kap * S.e_psi(D - 1, k - 1, j, 1) for j in range(k))
        gamma_k = R.lt(k) + R.c2D * R.h(k) * (C ** (D - 2 * k + 1) * E[k] + B_hat / 2)
        assert gamma_k == R.sh2 / nonneg_floor(S, k), ("identity K1", d, t, k)
        assert B_hat <= eta * R.phistar(k) - kap * R.phik(k)                   # Lemma M(v)
    assert E[1] == E[2] and all(E[k + 1] <= E[k] for k in range(1, D))          # Lemma M(i)
    assert all(E[2 * m - 1] == E[2 * m] for m in range(1, D // 2 + 1))
    assert all(R.r(k + 1) < R.r(k) and R.h(k + 1) > R.h(k) for k in range(1, D))        # (ii)
    assert all(R.lt(k + 1) < R.lt(k) and R.kap(k + 1) < R.kap(k) for k in range(1, D))  # (iii)
    assert all(R.phistar(k + 1) >= R.phistar(k) and R.phik(k + 1) <= R.phik(k) for k in range(1, D))  # (iv)


def signed_brute(t, D, k, lam, c):
    """Definition (50) at h = 0: max over x of E_eps[H_+ psi(eps.x) + H_- (2cosh(beta eps.x)/c^3 - 3/c^2)]."""
    S = Star(t, D)
    signs = list(product((1, -1), repeat=D))
    best = None
    for x in signs:
        total = Fr(0)
        for eps in signs:
            XU, XV = sum(eps[:k]), sum(eps[k:])
            H = (cosh_2beta(S, XU) + cosh_2beta(S, XV)) / 2 - lam * S.alpha / 4 * S.eta ** (-eps[0] * XU)
            L = sum(a * b for a, b in zip(eps, x))
            if H >= 0:
                total += H * S.psi(L)
            else:
                total += (-H) * (2 * S.cosh_even(L) / c ** 3 - Fr(3) / c ** 2)
        value = total / len(signs)
        best = value if best is None or value > best else best
    return best


def signed_k1(S, lam, c):
    """k = 1: H = (C + cosh(2 beta X_V))/2 - lam t cosh^2 beta depends on X_V only."""
    s = S.D - 1
    H = {m: (S.C + S.cosh_even(2 * (2 * m - s))) / 2 - lam * S.tc02 for m in range(s + 1)}
    values = []
    for j in range(s + 1):
        total = Fr(0)
        for n_plus in range(j + 1):
            for n_minus in range(s - j + 1):
                mult = comb(j, n_plus) * comb(s - j, n_minus)
                h = H[n_plus + n_minus]
                for e0 in (1, -1):
                    L = e0 + 2 * n_plus - 2 * n_minus - 2 * j + s
                    if h >= 0:
                        total += mult * h * S.psi(L)
                    else:
                        total += mult * (-h) * (2 * S.cosh_even(L) / c ** 3 - Fr(3) / c ** 2)
        values.append(total / 2 ** (s + 1))
    return max(values)


def main():
    start = time.monotonic()
    triples = 0
    for d, t in POINTS:
        S = Star(t, 2 * d)
        row = [nonneg_floor_all_j(S, k) for k in range(1, 2 * d + 1)]
        triples += len(row)
        assert min(row) == row[0] == Reduction(d, t).p_star
        print(f"PASS full h=0 row d={d} t={t}: Lemma J at all {len(row)} k, minimum at k=1 "
              f"({time.monotonic() - start:.1f}s)", flush=True)
    assert triples == 256
    for d, t in [(16, Fr(13, 125)), (17, Fr(1, 10))]:
        check_identity_k1(d, t)
        print(f"PASS identity (K1) and Lemma M at d={d} t={t} ({time.monotonic() - start:.1f}s)", flush=True)
    cases = 0
    for D in (4, 6):
        for k in (1, 2, 3):
            for t, lam, c in [(Fr(1, 10), Fr(10), Fr(6, 5)), (Fr(1, 4), Fr(3), Fr(11, 10)),
                              (Fr(11, 125), Fr(100000, 5860), Fr(119, 100))]:
                assert signed_residual(Star(t, D), k, lam, c)[0] == signed_brute(t, D, k, lam, c)
                cases += 1
    for d, t, floor, c in [(16, Fr(13, 125), Fr(689, 10000), Fr(6, 5)), (22, Fr(11, 125), Fr(586, 10000), Fr(119, 100))]:
        S = Star(t, 2 * d)
        assert signed_k1(S, 1 / floor, c) == signed_residual(S, 1, 1 / floor, c)[0]
    print(f"PASS signed count enumeration == brute force in {cases} cases; k=1 formula agrees "
          f"({time.monotonic() - start:.1f}s)", flush=True)
    S22 = Star(Fr(11, 125), 44)
    for k in (1, 2):
        ours = signed_residual(S22, k, 1 / Fr(586, 10000), Fr(119, 100))[0]
        assert ours == max(verify_signed_star.signed_residual(0, k))
    print(f"PASS M_(0,1), M_(0,2) at d=22 equal the dimension-22 certificate's values exactly "
          f"({time.monotonic() - start:.1f}s)")


if __name__ == "__main__":
    main()

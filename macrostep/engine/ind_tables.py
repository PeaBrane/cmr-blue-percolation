"""Macrostep route, independent implementation (code base R0), part 1: exact boost tables.

Written from the statements of the macrostep write-up (Lemma E) and global.md Lemma 3.4; it does not import the
L2 rig code.  Every table entry is an exact rational upper bound (a dyadic rational with denominator 2^SC,
obtained by rounding each exact b(r) UP to that grid; sums of such numbers are exact), and is converted to the
smallest binary64 number >= that rational.

Notation (macrostep write-up Sec. 3): t = tanh K' (exact rational (g_K-1)/(g_K+1)), q = d t/(1-2 d t) < 1,
  b(r)       = t^r M_r Gamma_r(t)   (global.md Lemma 3.4(b),(d)),  r >= 1,
  T1(F, r)   = sum_{k=0}^{c} b(F + |r - k|),                        F >= 1,
  Phi2(r, F) = max_{rho >= r} sum_{k=0}^{c} b(F + |rho - k|),       F >= 1,
  psi3(r, D) = sum_{s >= 3} Phi2(max(0, r - c(s-1)), max(s, D - s)).
Monotonicity used for clamping: b, Phi2, psi3 are nonincreasing in every argument; T1 is nonincreasing in F,
and in r for r >= c.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import numpy as np

SC = 300                     # dyadic scale: all exact table values are integers / 2^SC


def fup(x):
    """smallest float >= the rational x (exact comparison)."""
    x = Fr(x)
    f = float(x)
    while Fr(f) < x:
        f = float(np.nextafter(f, np.inf))
    return f


def ceil_scaled(x):
    """smallest integer N with N / 2^SC >= x."""
    x = Fr(x)
    return -((-x.numerator << SC) // x.denominator)


def int_to_fup(N):
    return fup(Fr(N, 1 << SC))


def largest_multinomial(d, r):
    """max over compositions (k_1..k_d) of r of r!/prod k_i!  (attained at the balanced composition)."""
    q, j = divmod(r, d)
    return Fr(factorial(r), factorial(q + 1) ** j * factorial(q) ** (d - j))


def Gamma(d, t, r):
    """Gamma_r(t) = 2^-d sum_j C(d,j) (1 - 2t(d-2j))^{-(r+1)}  (global.md Lemma 3.4(b))."""
    return sum(Fr(comb(d, j)) / (1 - 2 * t * (d - 2 * j)) ** (r + 1) for j in range(d + 1)) / 2 ** d


class BoostTables:
    def __init__(self, d, t, c, Rb=150, Fmax=300, Rmax=100, Dmax=60, S=260):
        t = Fr(t)
        self.d, self.t, self.c = d, t, c
        assert 0 < 2 * d * t < 1
        self.q = q = Fr(d) * t / (1 - 2 * d * t)
        assert q < 1
        be = [None] + [t ** r * largest_multinomial(d, r) * Gamma(d, t, r) for r in range(1, Rb + 1)]
        for r in range(1, Rb):
            assert be[r + 1] <= q * be[r], r               # Lemma 3.4(d) on the computed range
        top = max(Fmax + Rmax + 2 * c + 10, S + 2 * c + 10)
        # integer upper bounds B[r] >= 2^SC b(r); beyond Rb use b(r) <= b(Rb) q^{r-Rb}
        B = [0] * (top + 1)
        for r in range(1, top + 1):
            B[r] = ceil_scaled(be[r] if r <= Rb else be[Rb] * q ** (r - Rb))
        self.B = B
        self.bt = np.array([np.inf] + [int_to_fup(B[r]) for r in range(1, top + 1)])

        def row_sum(F, rho):
            return sum(B[F + abs(rho - k)] for k in range(c + 1))

        T1 = np.full((Fmax + 1, Rmax + 1), np.inf)
        for F in range(1, Fmax + 1):
            for r in range(Rmax + 1):
                T1[F, r] = int_to_fup(row_sum(F, r))
        self.T1 = T1

        def phi2_int(r, F):
            return max(row_sum(F, rho) for rho in range(r, max(r, c) + 1))
        self.phi2_int = phi2_int
        P2i = {}
        P2 = np.full((Rmax + 1, Fmax + 1), np.inf)
        for F in range(1, Fmax + 1):
            for r in range(Rmax + 1):
                P2i[(r, F)] = phi2_int(r, F)
                P2[r, F] = int_to_fup(P2i[(r, F)])
        self.PHI2 = P2
        # psi3: exact sum over 3 <= s <= S plus tail for s > S (valid when S >= D and c(S-1) >= r):
        #   term_s = Phi2(0, s) <= b(s) (1 + 2q/(1-q)),   sum_{s > S} b(s) <= b(S+1)/(1-q).
        g = 1 + 2 * q / (1 - q)
        tailI = ceil_scaled(g * Fr(B[S + 1], 1 << SC) / (1 - q))
        PS = np.full((Rmax + 1, Dmax + 1), np.inf)
        for r in range(Rmax + 1):
            for D in range(Dmax + 1):
                assert S >= D and c * (S - 1) >= r
                tot = 0
                for s in range(3, S + 1):
                    rr, FF = max(0, r - c * (s - 1)), max(s, D - s)
                    tot += P2i[(rr, FF)] if (rr <= Rmax and FF <= Fmax) else phi2_int(rr, FF)
                PS[r, D] = int_to_fup(tot + tailI)
        self.PSI3 = PS
        self.Fmax, self.Rmax, self.Dmax = Fmax, Rmax, Dmax
        assert Rmax >= c

    def psi3_00(self):
        return self.PSI3[0, 0]


def check_monotone(T):
    """float sanity of the monotonicity used for clamping."""
    assert np.all(np.diff(T.bt[1:]) <= 0)
    assert np.all(np.diff(T.T1[1:, :], axis=0) <= 0)
    assert np.all(np.diff(T.T1[1:, T.c:], axis=1) <= 0)
    assert np.all(np.diff(T.PHI2[:, 1:], axis=0) <= 0) and np.all(np.diff(T.PHI2[:, 1:], axis=1) <= 0)
    assert np.all(np.diff(T.PSI3, axis=0) <= 0) and np.all(np.diff(T.PSI3, axis=1) <= 0)
    return True

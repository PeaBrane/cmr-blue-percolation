"""Exact certificates for the single-floor oriented comparison (manuscript appendix).

Theorem (single-floor oriented comparison): if a uniform fresh-star floor p
exceeds F_d = 1 - 1/G_d, where G_d is the collision Green function of two
uniform oriented walks, then both overlap signs percolate and
nu(o <-> infinity blue, q_o = s) >= (1/2)(p - F_d)/(p (1 - F_d)). The row
h = 0 of the fresh-star floor, p_0(2d, beta), is such a floor.

This program checks, in exact rational arithmetic (research-note names in
brackets):

(G)   collision Green function bounds [Lemma G]: u_0..u_3 in closed form,
      exact head sums over n < 4d plus the maximal-atom tail (part (e)), the
      eight-term closed form hat G_d (part (f)) and the four-term bound on
      G_22 (part (g));
(K)   the row h = 0 as a single ratio [identity (K1)]: the closed form p_* at
      k = 1 equals the nonnegative-residual value p(1), and the block cover of
      the reduction to k = 1 by blocks [Proposition K] gives p(k) >= p(1) for
      every overlap count k, so p_0(2d, beta) = p_*;
(S)   the corollary for dimensions 16 to 20 and 22 [Corollaries S16-S20, S22
      route (ii)]: p_* > F_d at seven points, with the table displays;
(S22) the signed refinement at d = 22 [route (i)]: the signed residuals
      M_(0,1), M_(0,2) < 0 at lambda = 1/0.0586, c = 119/100, the nonnegative
      residual min_(k >= 3) p(k) > 0.0586, and G_22 < 1/(1 - 0.0586);
(S16) the signed floor 0.0689 at d = 16, t = 13/125;
(P)   the closed-form star floor and the certificate in dimension 25
      [Proposition P]: min(p_pair(1), p_pair(2d)) exceeds 1 - 1/hat G_d at
      d = 25 and 26 (and at d = 28, a research-note value), with no finite
      sums.

The cavity maximum of the nonnegative residual is evaluated at the balanced
index, which the lemma "the balanced tilt maximizes" [Lemma J] justifies;
check_lemmas.py in this directory evaluates every index. All decisive
comparisons are assertions on exact rationals; decimals in the output are
directed summaries. Only the standard library is used.
"""

from fractions import Fraction as Fr
from itertools import product
from math import comb, factorial
import hashlib
import json
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

POINTS = [(16, Fr(13, 125)), (16, Fr(51, 500)), (17, Fr(1, 10)), (18, Fr(12, 125)),
          (19, Fr(7, 75)), (20, Fr(1, 11)), (22, Fr(11, 125))]
# Proposition K block covers of {2, ..., 2d} (found by an exact greedy search).
BLOCKS = {
    (16, Fr(13, 125)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 7), (8, 13), (14, 32)],
    (16, Fr(51, 500)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 7), (8, 13), (14, 32)],
    (17, Fr(1, 10)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 7), (8, 12), (13, 32), (33, 34)],
    (18, Fr(12, 125)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 7), (8, 11), (12, 30), (31, 36)],
    (19, Fr(7, 75)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 7), (8, 11), (12, 30), (31, 38)],
    (20, Fr(1, 11)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 7), (8, 10), (11, 25), (26, 40)],
    (22, Fr(11, 125)): [(2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 8), (9, 13), (14, 37), (38, 44)],
}
# Displayed bounds (short-route.md Sec. 10.2; appendix tables of the revised manuscript):
# p_* >=, p(2) >=, p_*/F_d >= (head, hat G), theta_* = (1/2)(p_* - F)/(p_* (1 - F)) >= (head, hat G),
# smallest block slack Gamma_1 - B >=.
SHOWN = {
    (16, Fr(13, 125)): ("0.0679971076", "0.0690100", "1.0146", "1.0080", "0.0077", "0.0042", "0.0267"),
    (16, Fr(51, 500)): ("0.0680421165", "0.0690180", "1.0153", "1.0087", "0.0080", "0.0046", "0.0065"),
    (17, Fr(1, 10)): ("0.0659393703", "0.0668298", "1.0502", "1.0447", "0.0255", "0.0228", "0.0073"),
    (18, Fr(12, 125)): ("0.0640504639", "0.0648326", "1.0845", "1.0798", "0.0414", "0.0393", "0.0052"),
    (19, Fr(7, 75)): ("0.0623013944", "0.0630062", "1.1174", "1.1135", "0.0556", "0.0539", "0.0064"),
    (20, Fr(1, 11)): ("0.0606881634", "0.0613263", "1.1494", "1.1460", "0.0686", "0.0672", "0.0044"),
    (22, Fr(11, 125)): ("0.0577808419", "0.0583248", "1.2103", "1.2077", "0.0912", "0.0903", "0.0099"),
}
GAMMA_1_RANGE = (Fr("3.020"), Fr("3.093"))
P2_OVER_P1_RANGE = (Fr("1.0094"), Fr("1.0149"))
# Binomial sums used by the block check: E_1, E_2d and two per block (appendix block table).
SUMS = {(16, Fr(13, 125)): 16, (16, Fr(51, 500)): 16, (17, Fr(1, 10)): 18, (18, Fr(12, 125)): 18,
        (19, Fr(7, 75)): 18, (20, Fr(1, 11)): 18, (22, Fr(11, 125)): 20}
# Lemma G displays (short-route.md Sec. 8): head G_d <=, F_d <= (head), hat G_d <=, F_d <= (hat G), F_d >=.
SHOWN_G = {
    16: ("1.0718281723", "0.06701464", "1.07233323", "0.067455", "0.06701463"),
    17: ("1.0669890897", "0.06278330", "1.06736610", "0.063115", "0.06278329"),
    18: ("1.0627639544", "0.05905729", "1.06305199", "0.059313", "0.05905728"),
    19: ("1.0590423693", "0.05575072", "1.05926631", "0.055951", "0.05575071"),
    20: ("1.0557390227", "0.05279622", "1.05591558", "0.052955", "0.05279621"),
    22: ("1.0501326858", "0.04773939", "1.05024614", "0.047843", "0.04773938"),
}
# Proposition P (short-route.md Sec. 9.4): p_pair(1) >=, p_pair(2d) >=, F_d <= (hat G), ratio >=,
# (1/2)(p - F)/(p (1 - F)) >=.
SHOWN_P = {
    (25, Fr(27, 400)): ("0.04196433", "0.05133661", "0.04180338", "1.0038", "0.0020"),
    (26, Fr(33, 500)): ("0.04113940", "0.05037252", "0.04011778", "1.0254", "0.0129"),
    (28, Fr(8, 125)): ("0.03962652", "0.04860885", "0.03712607", "1.0673", "0.0327"),
}
FLOOR_22, TANGENT_22 = Fr(586, 10000), Fr(119, 100)
FLOOR_16, TANGENT_16 = Fr(689, 10000), Fr(6, 5)


def dec_up(x, digits):
    x = Fr(x)
    return _decimal(-((-x.numerator * 10 ** digits) // x.denominator), digits)


def dec_dn(x, digits):
    x = Fr(x)
    return _decimal((x.numerator * 10 ** digits) // x.denominator, digits)


def _decimal(q, digits):
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


# ------------------------------------------------------------------ (G) collision Green function
def u_exact(d, N):
    """u_n(d) = sum over k_1+..+k_d = n of (n!/(k_1!..k_d!) d^-n)^2, n = 0..N."""
    a = [Fr(1, factorial(k) ** 2) for k in range(N + 1)]
    p = [Fr(1)] + [Fr(0)] * N
    for _ in range(d):
        q = [Fr(0)] * (N + 1)
        for i, pi in enumerate(p):
            if pi:
                for j in range(N + 1 - i):
                    q[i + j] += pi * a[j]
        p = q
    return [p[n] * factorial(n) ** 2 / Fr(d) ** (2 * n) for n in range(N + 1)]


def m_atom(d, n):
    """Largest atom of Multinomial(n; 1/d, ..., 1/d)."""
    a, j = divmod(n, d)
    return Fr(factorial(n), d ** n * factorial(a + 1) ** j * factorial(a) ** (d - j))


def tail_bound(d, A):
    """sum_(n >= dA) u_n <= d M_(dA) (1 + 2(A+1)/(d-3)) (Lemma G(e); d >= 4, A >= 1)."""
    assert d >= 4 and A >= 1
    return d * m_atom(d, d * A) * (1 + Fr(2 * (A + 1), d - 3))


def g_closed(d):
    """Eight-term closed-form upper bound hat G_d (Lemma G(f); d >= 6)."""
    assert d >= 6
    return (1 + Fr(1, d) + Fr(2 * d - 1, d ** 3) + Fr(6 * d * d - 9 * d + 4, d ** 5)
            + Fr(24, d ** 4) + Fr(120, d ** 5) + Fr((d - 6) * 720, d ** 6)
            + d * Fr(factorial(d), d ** d) * (1 + Fr(4, d - 3)))


def g_head(d, A=4):
    N = d * A - 1
    u = u_exact(d, N)
    assert u[0] == 1 and u[1] == Fr(1, d) and u[2] == Fr(2 * d - 1, d ** 3)
    assert u[3] == Fr(6 * d * d - 9 * d + 4, d ** 5)
    for n in range(N + 1):
        assert u[n] <= m_atom(d, n)
        assert n == 0 or m_atom(d, n) <= m_atom(d, n - 1)
    head = sum(u)
    return head, head + tail_bound(d, A)


def g22_four_term():
    """Lemma G(g): u_2 <= 2/d^2, u_n <= 6/d^3 for 3 <= n <= 21, tail with A = 1."""
    d = 22
    return 1 + Fr(1, d) + Fr(2, d * d) + (d - 3) * Fr(6, d ** 3) + tail_bound(d, 1)


# ------------------------------------------------------------------ (L) the h = 0 fresh-star floor
class Star:
    def __init__(self, t, degree):
        self.t = t = Fr(t)
        self.D = degree
        self.eta = (1 - t) / (1 + t)           # e^{-2 beta}
        self.C = (1 + t * t) / (1 - t * t)     # cosh 2 beta
        self.v = 2 * t / (1 + t * t)           # tanh 2 beta
        self.u = (1 + self.v) / 2
        self.alpha = 1 - self.eta ** 2
        self.tc02 = t / (1 - t * t)            # t cosh^2 beta = sinh(2 beta)/2
        self._psi, self._dist = {}, {}

    def psi(self, x):
        """sech^2(beta x) = 4 eta^|x| / (1 + eta^|x|)^2."""
        x = abs(x)
        if x not in self._psi:
            e = self.eta ** x
            self._psi[x] = 4 * e / (1 + e) ** 2
        return self._psi[x]

    def cosh_even(self, L):
        """cosh(beta L) for even L."""
        assert L % 2 == 0
        e = self.eta ** (abs(L) // 2)
        return (1 / e + e) / 2

    def dist(self, n, l, j):
        """Law of the number of +1 among n signs: j with P(+) = u, l - j with P(+) = 1 - u, n - l fair."""
        key = (n, l, j)
        if key not in self._dist:
            u, P = self.u, [Fr(1)]
            for a, b in [(1 - u, u)] * j + [(u, 1 - u)] * (l - j) + [(Fr(1, 2), Fr(1, 2))] * (n - l):
                Q = [Fr(0)] * (len(P) + 1)
                for i, c in enumerate(P):
                    Q[i] += c * a
                    Q[i + 1] += c * b
                P = Q
            self._dist[key] = P
        return self._dist[key]

    def e_psi(self, n, l, j, shift=0):
        return sum(c * self.psi(2 * m - n + shift) for m, c in enumerate(self.dist(n, l, j)))


def h_ell(S, l):
    t, v = S.t, S.v
    return ((1 + t) * (1 + t * v) ** l + (1 - t) * (1 - t * v) ** l) / 2


def nonneg_floor(S, k):
    """Nonnegative-residual floor p(k) at h = 0 (paper App. A, eq. (49)); T_s at the balanced index."""
    D, t, C, v, u, eta = S.D, S.t, S.C, S.v, S.u, S.eta
    s = D - k
    N = t * (1 - t * t) ** (D - 1) * C ** (k - 1) / h_ell(S, k - 1) ** 2
    cs = Fr(1) if s % 2 == 0 else C
    lam = (1 + eta ** (2 * k) + 2 * cs * eta ** k) / S.alpha
    Ts = C ** s * S.e_psi(D, s, s // 2)
    B = max((u - lam * v) * S.e_psi(D - 1, k - 1, j, 1) + (1 - u) * S.e_psi(D - 1, k - 1, j, -1)
            for j in range(k))
    M = Ts / 2 + C ** k / 2 * B
    assert M >= 0
    return N / (lam * N + M)


def cosh_2beta(S, X):
    """cosh(2 beta X) = (eta^-X + eta^X)/2."""
    return (S.eta ** (-X) + S.eta ** X) / 2


def signed_residual(S, k, lam, c):
    """max over cavity orientations of paper eq. (50) at h = 0, overlap count k (count enumeration)."""
    D = S.D
    s = D - k
    K3, K2 = Fr(3) / c ** 2, 2 / c ** 3
    H_cache, neg_kernel = {}, {}
    best, orientations = None, 0
    for xi in product((1, -1), repeat=k - 1):
        for j in range(s + 1):
            orientations += 1
            total = Fr(0)
            for e in product((1, -1), repeat=k - 1):
                XU = 1 + sum(e)
                LU = 1 + sum(a * b for a, b in zip(xi, e))
                a_val = S.alpha / 4 * S.eta ** (-XU)
                for n_plus in range(j + 1):
                    for n_minus in range(s - j + 1):
                        XV = 2 * (n_plus + n_minus) - s
                        L = LU + 2 * n_plus - 2 * n_minus + s - 2 * j
                        if (XU, XV) not in H_cache:
                            H_cache[XU, XV] = (cosh_2beta(S, XU) + cosh_2beta(S, XV)) / 2 - lam * a_val
                        H = H_cache[XU, XV]
                        weight = comb(j, n_plus) * comb(s - j, n_minus)
                        if H >= 0:
                            total += weight * H * S.psi(L)
                        else:
                            if L not in neg_kernel:
                                neg_kernel[L] = K2 * S.cosh_even(L) - K3
                            total += weight * (-H) * neg_kernel[L]
            value = total / 2 ** (D - 1)
            if best is None or value > best:
                best = value
    return best, orientations


class Reduction:
    """Identity (K1), the closed form p_* and the block bound of Proposition K."""

    def __init__(self, d, t):
        self.D = D = 2 * d
        self.S = S = Star(t, D)
        self.C, self.eta = S.C, S.eta
        self.c2D = (1 / (1 - S.t ** 2)) ** D
        self.sh2 = 2 * S.t / (1 - S.t ** 2)
        self.E1 = self.pair_sum_e1(d)
        assert self.E1 == S.e_psi(D, D - 1, (D - 1) // 2)
        self.Gf = self.fair_sum(D, 0)
        self.N1 = 2 * self.C + self.c2D * (self.C ** (D - 1) * self.E1 - self.C * self.Gf)
        self.p_star = self.sh2 / self.N1

    def fair_sum(self, m, shift):
        return sum(Fr(comb(m, i), 2 ** m) * self.S.psi(2 * i - m + shift) for i in range(m + 1))

    def pair_sum_e1(self, d):
        """E psi(P_(d-1) + F_2), P a sum of d-1 balanced pairs in {-2, 0, 2}."""
        w = self.S.u * (1 - self.S.u)
        dist = {0: Fr(1)}
        for _ in range(d - 1):
            new = {}
            for x, c in dist.items():
                for dx, pr in ((-1, w), (0, 1 - 2 * w), (1, w)):
                    new[x + dx] = new.get(x + dx, 0) + c * pr
            dist = new
        F2 = {-2: Fr(1, 4), 0: Fr(1, 2), 2: Fr(1, 4)}
        return sum(c * pr * self.S.psi(2 * x + y) for x, c in dist.items() for y, pr in F2.items())

    def cs(self, k):
        return Fr(1) if (self.D - k) % 2 == 0 else self.C

    def lt(self, k):
        return (1 / self.eta + self.eta ** (2 * k - 1)) / 2 + self.cs(k) * self.eta ** (k - 1)

    def kap(self, k):
        return self.eta ** (2 * k - 1) + 2 * self.cs(k) * self.eta ** (k - 1)

    def h(self, k):
        return h_ell(self.S, k - 1) ** 2

    def r(self, k):
        return self.h(k) / self.C ** (2 * (k - 1))

    def phistar(self, k):
        """max of phi_k(y) = E psi(y + F_(2d-k)), attained at y = k mod 2."""
        return self.fair_sum(self.D - k, k % 2)

    def phik(self, k):
        """phi_k(k)."""
        return self.fair_sum(self.D - k, k)

    def block(self, ka, kb):
        """B(ka, kb) >= Gamma_k = sinh(2 beta)/p(k) for every k in [ka, kb]."""
        return self.lt(ka) + self.c2D * (self.r(ka) * self.C ** (self.D - 1) * self.E1
                                         + (self.h(kb) * self.eta * self.phistar(kb)
                                            - self.h(ka) * self.kap(kb) * self.phik(kb)) / 2)


def root_constant(p, F):
    return (p - F) / (p * (1 - F)) / 2


# ------------------------------------------------------------------ (P) the closed-form floor
def p_pair_ends(d, t):
    D = 2 * d
    C, v = (1 + t * t) / (1 - t * t), 2 * t / (1 + t * t)
    H = ((1 + t) * (1 + t * v) ** (D - 1) + (1 - t) * (1 - t * v) ** (D - 1)) / 2
    p1 = 2 * t * (1 - t * t) ** (D - 1) / (C + C ** (D - 1))
    pD = 2 * t * (1 - t * t) ** (D - 1) * C ** (D - 1) / (H * H * (C ** D + 1))
    return p1, pD


def main():
    start = time.monotonic()
    report, vector = {}, []
    green = {}
    for d in sorted({d for d, _ in POINTS}):
        head, up = g_head(d)
        hat = g_closed(d)
        shown = SHOWN_G[d]
        assert up <= Fr(shown[0]) and 1 - 1 / up <= Fr(shown[1])
        assert hat <= Fr(shown[2]) and 1 - 1 / hat <= Fr(shown[3]) and 1 - 1 / head >= Fr(shown[4])
        assert (1 - 1 / up) - (1 - 1 / head) <= Fr(1, 10 ** 8)            # F_d determined to within 1e-8
        green[d] = (head, up, hat)
        vector += [head, up, hat]
        report[f"G_{d}"] = {"head_upper": dec_up(up, 10), "F_d_upper_head": dec_up(1 - 1 / up, 8),
                            "hat_G_upper": dec_up(hat, 8), "F_d_upper_hat": dec_up(1 - 1 / hat, 6),
                            "F_d_lower": dec_dn(1 - 1 / head, 8)}
    print(f"PASS collision Green function bounds for d = 16-20, 22 ({time.monotonic() - start:.1f}s)",
          flush=True)

    corollaries, p2_over_p1 = [], []
    for d, t in POINTS:
        R = Reduction(d, t)
        p1 = nonneg_floor(R.S, 1)
        assert R.p_star == p1                                  # closed form (K1) at k = 1
        blocks = BLOCKS[d, t]
        assert blocks[0][0] == 2 and blocks[-1][1] == 2 * d
        assert all(b[1] + 1 == c[0] and b[0] <= b[1] for b, c in zip(blocks, blocks[1:] + [(2 * d + 1, 0)]))
        assert 2 + 2 * len(blocks) == SUMS[d, t]
        slacks = [R.N1 - R.block(ka, kb) for ka, kb in blocks]
        assert min(slacks) >= 0                                # Proposition K: p(k) >= p(1) for all k
        p2 = nonneg_floor(R.S, 2)
        p2_over_p1.append(p2 / R.p_star)
        head, up, hat = green[d]
        F_head, F_hat = 1 - 1 / up, 1 - 1 / hat
        assert R.p_star > F_head and R.p_star > F_hat          # Theorem S hypothesis p_0 > F_d
        shown = [Fr(x) for x in SHOWN[d, t]]
        assert R.p_star >= shown[0] and p2 >= shown[1]
        assert R.p_star / F_head >= shown[2] and R.p_star / F_hat >= shown[3]
        assert root_constant(R.p_star, F_head) >= shown[4] and root_constant(R.p_star, F_hat) >= shown[5]
        assert min(slacks) >= shown[6] and GAMMA_1_RANGE[0] <= R.N1 <= GAMMA_1_RANGE[1]
        vector += [R.p_star, p2] + slacks
        corollaries.append({"d": d, "t": str(t), "p_star_lower": dec_dn(R.p_star, 10), "p2_lower": dec_dn(p2, 7),
                            "blocks": len(blocks), "min_block_slack_lower": dec_dn(min(slacks), 4),
                            "ratio_head_lower": dec_dn(R.p_star / F_head, 4),
                            "ratio_hat_lower": dec_dn(R.p_star / F_hat, 4),
                            "root_constant_head_lower": dec_dn(root_constant(R.p_star, F_head), 4),
                            "root_constant_hat_lower": dec_dn(root_constant(R.p_star, F_hat), 4)})
        print(f"PASS single-floor corollary d={d} t={t}: p_0 = p_* >= {dec_dn(R.p_star, 10)} > F_{d}, "
              f"{len(blocks)} blocks ({time.monotonic() - start:.1f}s)", flush=True)
    assert P2_OVER_P1_RANGE[0] <= min(p2_over_p1) and max(p2_over_p1) <= P2_OVER_P1_RANGE[1]
    report["corollaries"] = corollaries
    report["p2_over_p1"] = [dec_dn(min(p2_over_p1), 4), dec_up(max(p2_over_p1), 4)]

    # Corollary S22, route (i): the floor 0.0586 from signed k = 1, 2 and nonnegative k >= 3.
    S22 = Star(Fr(11, 125), 44)
    assert S22.eta == Fr(57, 68)
    signed = [signed_residual(S22, k, 1 / FLOOR_22, TANGENT_22) for k in (1, 2)]
    assert all(value < 0 for value, _ in signed) and [n for _, n in signed] == [44, 86]
    assert signed[0][0] <= Fr("-0.000374848573757") and signed[1][0] <= Fr("-0.019889903399663")
    row22 = [nonneg_floor(S22, k) for k in range(1, 45)]
    rest = min(row22[2:])
    assert rest > FLOOR_22 and rest >= Fr("0.059010354650604") and rest == row22[2]
    assert min(row22) == row22[0] == Reduction(22, Fr(11, 125)).p_star
    assert Fr("0.0583248") <= row22[1] < Fr("0.0583249") < FLOOR_22    # p(2) alone is below the floor 0.0586
    four = g22_four_term()
    assert tail_bound(22, 1) < Fr("8.8e-8")
    F4 = 1 - 1 / four
    assert four < 1 / (1 - FLOOR_22) and four < Fr("1.0602932") and F4 < Fr("0.0568646")
    assert root_constant(FLOOR_22, F4) >= Fr("0.0157")
    _, _, hat22 = green[22]
    assert 2 * root_constant(Reduction(22, Fr(11, 125)).p_star, 1 - 1 / hat22) >= Fr("0.1806")
    vector += [value for value, _ in signed] + row22 + [four]
    report["S22_route_i"] = {"M01_upper": dec_up(signed[0][0], 15), "orientations_01": signed[0][1],
                             "M02_upper": dec_up(signed[1][0], 15), "orientations_02": signed[1][1],
                             "min_k_ge_3_lower": dec_dn(rest, 15), "four_term_G22_upper": dec_up(four, 10),
                             "F22_upper": dec_up(F4, 7), "root_constant_lower": dec_dn(root_constant(FLOOR_22, F4), 4)}
    print(f"PASS signed refinement d=22: floor 0.0586 > F_22, G_22 <= {dec_up(four, 10)} "
          f"({time.monotonic() - start:.1f}s)", flush=True)

    # Optional signed k = 1 floor at d = 16: 0.0689.
    S16 = Star(Fr(13, 125), 32)
    m16, n16 = signed_residual(S16, 1, 1 / FLOOR_16, TANGENT_16)
    rest16 = min(nonneg_floor(S16, k) for k in range(2, 33))
    _, up16, hat16 = green[16]
    assert m16 < 0 and m16 <= Fr("-0.00136416316276") and n16 == 32 and rest16 >= Fr("0.0690100") > FLOOR_16
    assert FLOOR_16 / (1 - 1 / up16) >= Fr("1.0281") and FLOOR_16 / (1 - 1 / hat16) >= Fr("1.0214")
    vector += [m16, rest16]
    report["S16_signed"] = {"M01_upper": dec_up(m16, 15), "orientations": n16,
                            "min_k_ge_2_lower": dec_dn(rest16, 7),
                            "ratio_head_lower": dec_dn(FLOOR_16 / (1 - 1 / up16), 4),
                            "ratio_hat_lower": dec_dn(FLOOR_16 / (1 - 1 / hat16), 4)}
    print(f"PASS signed refinement d=16: floor 0.0689 > F_16 ({time.monotonic() - start:.1f}s)", flush=True)

    # Proposition P: certificate-free dimensions.
    closed = []
    for (d, t), shown in SHOWN_P.items():
        p1, pD = p_pair_ends(d, t)
        hat = g_closed(d)
        F = 1 - 1 / hat
        assert min(p1, pD) == p1 > F
        assert p1 >= Fr(shown[0]) and pD >= Fr(shown[1]) and F <= Fr(shown[2]) and p1 / F >= Fr(shown[3])
        assert root_constant(p1, F) >= Fr(shown[4])
        vector += [p1, pD, hat]
        closed.append({"d": d, "t": str(t), "p_pair_1_lower": dec_dn(p1, 8), "p_pair_2d_lower": dec_dn(pD, 8),
                       "hat_G_upper": dec_up(hat, 8), "F_d_upper": dec_up(F, 8), "ratio_lower": dec_dn(p1 / F, 4),
                       "root_constant_lower": dec_dn(root_constant(p1, F), 4)})
    report["proposition_P"] = closed
    print(f"PASS closed-form star floor above F_d at d = 25, 26, 28 ({time.monotonic() - start:.1f}s)", flush=True)
    report["value_vector_sha256"] = hashlib.sha256("\n".join(str(x) for x in vector).encode()).hexdigest()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

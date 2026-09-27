"""Macrostep route, code base R0t: tightened rigorous forward return probabilities u_n(0) and tails T_u(N).

u_n(0) = P(Y_n = 0) for the difference Y_n of two independent uniform oriented walks in Z^f, i.e.
    u_n(0) = V_f(n) / f^(2n),   V_f(n) = sum_{k in N^f, |k| = n} multinom(n; k)^2   (exact integers).
Block recursion (exact): V_1(n) = 1, V_2(n) = C(2n, n), V_{a+b}(n) = sum_m C(n, m)^2 V_a(m) V_b(n - m).
Each u_n(0), n <= NEX, is stored as the integer U_n = ceil(2^SC V_f(n) / f^(2n)) (so U_n / 2^SC >= u_n(0)).

Tails beyond NEX (macrostep write-up, Lemma R.1):
  f = 4 (2+2 split):  u_n(0) <= C4 n^(-3/2) + 2 exp(-2 eps^2 n),  C4 = 2 pi^(-3/2) (1 - 4 eps^2)^(-1/2),
                      valid for n >= 1/(1/2 - eps);  sum_{n >= N1} n^(-3/2) <= N1^(-3/2) + 2 N1^(-1/2).
  f >= 5 (atom bound, global.md Lemma 4.2(ii) as in R0):
                      sum_{n >= f A1} u_n(0) <= f sqrt(f) (2 pi)^(-s) (A1^(-s) + A1^(1-s)/(s-1)),  s = (f-1)/2.
Every constant is an exact rational upper bound; no floating point enters.

Usage:  python tail_tight.py f NEX out.pkl       (builds the cache: U_n for n <= NEX, and the tail bound)
"""
import sys, time, pickle
from fractions import Fraction as Fr
from math import factorial, isqrt

SC = 256
PI_LO = Fr(314159265358979, 10 ** 14)         # < pi


def sqrt_lo(x, digits=30):
    x = Fr(x)
    S = 10 ** digits
    r = Fr(isqrt(x.numerator * S * S // x.denominator), S)
    assert r * r <= x
    return r


def sqrt_hi(x, digits=30):
    r = sqrt_lo(x, digits) + Fr(1, 10 ** digits)
    assert r * r >= x
    return r


def exp_neg_hi(x, K=400):
    """rational upper bound for e^(-x), x >= 0:  e^(-x) <= 1 / sum_{k <= K} x^k / k!."""
    x = Fr(x)
    assert x >= 0
    s, term = Fr(0), Fr(1)
    for k in range(K + 1):
        s += term
        term = term * x / (k + 1)
    return 1 / s


# ------------------------------------------------------------------------------------ exact V_f(n)
def V_table(f, NEX, log=None):
    """[V_f(0), ..., V_f(NEX)] as exact Python integers (block recursion)."""
    t0 = time.time()
    # binomial squares are generated row by row; V_2 is central binomials
    V1 = [1] * (NEX + 1)
    V2 = [1] * (NEX + 1)
    for n in range(1, NEX + 1):
        V2[n] = V2[n - 1] * 2 * (2 * n - 1) // n          # C(2n, n)
    need3 = f in (5, 6)
    if need3:
        V3 = [0] * (NEX + 1)
    Vf = [0] * (NEX + 1)
    row = [1]
    for n in range(NEX + 1):
        if n > 0:
            row = [1] + [row[k - 1] + row[k] for k in range(1, n)] + [1]
        sq = [r * r for r in row]
        if need3:
            V3[n] = sum(sq[m] * V2[n - m] for m in range(n + 1))       # V_3 = V_1 (x) V_2
        if f == 4:
            Vf[n] = sum(sq[m] * V2[m] * V2[n - m] for m in range(n + 1))
        elif f == 5:
            Vf[n] = sum(sq[m] * V2[m] * V3[n - m] for m in range(n + 1))
        elif f == 6:
            Vf[n] = sum(sq[m] * V3[m] * V3[n - m] for m in range(n + 1))
        else:
            raise ValueError(f)
        if log and n % 500 == 0:
            print(f"  f={f} n={n} [{time.time()-t0:.0f}s]", file=log, flush=True)
    return Vf


def V_brute(f, n):
    """sum over compositions of multinom^2 (validation, small n)."""
    import itertools
    tot = 0
    for k in itertools.product(range(n + 1), repeat=f - 1):
        s = sum(k)
        if s <= n:
            ks = list(k) + [n - s]
            m = factorial(n)
            for x in ks:
                m //= factorial(x)
            tot += m * m
    return tot


def U_from_V(f, Vf):
    """U_n = ceil(2^SC V_f(n) / f^(2n))."""
    out = []
    for n, v in enumerate(Vf):
        den = f ** (2 * n)
        out.append(-((-(v << SC)) // den))
    return out


# ------------------------------------------------------------------------------------ analytic tails
def tail_split4(N1, eps=Fr(1, 20)):
    """rational >= sum_{n >= N1} u_n(0) for f = 4 (Lemma R.1)."""
    eps = Fr(eps)
    assert 0 < eps < Fr(1, 2) and N1 >= 1 / (Fr(1, 2) - eps)
    pi32_inv = 1 / (PI_LO * sqrt_lo(PI_LO))                      # >= pi^(-3/2)
    C4 = 2 * pi32_inv / sqrt_lo(1 - 4 * eps * eps)
    s32 = 1 / (Fr(N1) * sqrt_lo(N1)) + 2 / sqrt_lo(N1)          # >= sum_{n >= N1} n^(-3/2)
    x = 2 * eps * eps
    geo = 2 * exp_neg_hi(x * N1) * (1 + 1 / x)                  # 1/(1 - e^-x) <= 1 + 1/x
    return C4 * s32 + geo


def robbins_tail(f, A1):
    """rational >= sum_{n >= f A1} u_n(0), f >= 5 (the R0 bound, global.md Lemma 4.2(ii))."""
    s2 = f - 1
    assert s2 > 2

    def pow_neg(x, e2):
        v = Fr(1) / Fr(x) ** (e2 // 2)
        if e2 % 2:
            v /= sqrt_lo(x)
        return v
    pref = sqrt_hi(f) * pow_neg(2 * PI_LO, s2)
    ssum = pow_neg(A1, s2) + pow_neg(A1, s2 - 2) / (Fr(s2, 2) - 1)
    return f * pref * ssum


EPS_GRID = [Fr(k, 40) for k in range(1, 16)]          # eps in {1/40, ..., 3/8}


def tail_from(f, N1):
    """rational >= sum_{n >= N1} u_n(0)  (f = 4: best split bound over EPS_GRID)."""
    if f == 4:
        return min(tail_split4(N1, e) for e in EPS_GRID if N1 >= 1 / (Fr(1, 2) - e))
    assert N1 % f == 0
    return robbins_tail(f, N1 // f)


class TightTails:
    """loads a cache built by main(); T_u(N) = sum_{N < n <= NEX} U_n/2^SC + tail_from(f, NEX+1)."""

    def __init__(self, path):
        o = pickle.load(open(path, "rb"))
        self.f, self.NEX, self.U = o['f'], o['NEX'], o['U']
        assert o['SC'] == SC and len(self.U) == self.NEX + 1
        self.tail = tail_from(self.f, self.NEX + 1)
        self.u0 = [Fr(u, 1 << SC) for u in self.U]

    def T_u(self, N):
        assert N < self.NEX
        return Fr(sum(self.U[N + 1:]), 1 << SC) + self.tail


def main():
    f, NEX, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    if f >= 5:
        assert (NEX + 1) % f == 0
    t0 = time.time()
    Vs = V_table(f, 6)
    for n in range(0, 7):
        assert Vs[n] == V_brute(f, n), n
    Vf = V_table(f, NEX, log=sys.stdout)
    U = U_from_V(f, Vf)
    pickle.dump(dict(f=f, NEX=NEX, SC=SC, U=U), open(out, "wb"))
    tail = tail_from(f, NEX + 1)
    s24 = Fr(sum(U[25:]), 1 << SC)
    s60 = Fr(sum(U[61:]), 1 << SC)
    s1500 = Fr(sum(U[25:1501]), 1 << SC) if NEX >= 1500 else None
    print(f"f={f} NEX={NEX}: brute-force check n<=6 OK; exact partial sums (upper, grid 2^-{SC}): "
          f"sum_(24,NEX] u_n(0) <= {float(s24):.10e}, sum_(60,NEX] <= {float(s60):.10e}"
          + (f", sum_(24,1500] <= {float(s1500):.10e}" if s1500 is not None else "")
          + f"; tail_(>NEX) <= {float(tail):.6e}; T_u(24) <= {float(s24 + tail):.6e}; T_u(60) <= {float(s60 + tail):.6e}"
          f" [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()

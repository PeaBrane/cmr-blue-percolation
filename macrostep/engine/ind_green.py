"""Macrostep route, code base R0: Green-function bounds for the pair chain (independent of the L2 rig code).

The pair chain S = (A, D) has independent components (macrostep write-up Lemma T.3):
  D_n = Y_n, the difference of two independent uniform oriented walks in Z^f   (forward part),
  A_n = X_n - X'_n, X_n a sum of n iid word endpoints                           (lateral part),
so G((A'',D''), (A',D')) = sum_n u_n(D' - D'') l_n(A' - A'').

Forward part: u_n(delta) = f^{-2n} sum_{k, |k| = n} multinom(n; k) multinom(n; k - delta), computed EXACTLY
(fractions) by a coordinate-by-coordinate convolution of the sequences j -> 1/(j! (j - delta_q)!).
Tail: u_n(delta) <= u_n(0) <= max_k P(multinom(n; 1/f) = k) <= sqrt(f) (2 pi floor(n/f))^{-(f-1)/2}
(Cauchy-Schwarz, monotonicity of the largest atom, Robbins' Stirling bounds).

Lateral part: p_n = law of X_n, computed in binary64 by direct convolution (round to nearest) with the step law
rounded UP; l_n(z) = sum_x p_n(x) p_n(x - z).  A-priori error bound (macrostep write-up Sec. 7.2): the exact values are
<= computed * (1 + LAT_INFL) + 1e-300, LAT_INFL = 1e-8 >> (1-gamma_{2E})^{-2n}(1-gamma_K)^{-1} - 1, where E <= 63 is
the number of distinct endpoints and K <= 1e8 the number of terms in an autocorrelation sum.
Tail: l_n(z) <= l_n(0) <= l_{N0}(0) for n >= N0 (Young and Cauchy-Schwarz).
"""
from fractions import Fraction as Fr
from math import factorial, isqrt
import numpy as np
from ind_tables import fup

LAT_INFL = 1e-8


# ------------------------------------------------------------------------------------------- forward, exact
def u_exact_table(f, N, delta):
    """[u_0(delta), ..., u_N(delta)] exactly."""
    fac = [factorial(j) for j in range(N + max(abs(x) for x in delta) + 2)]
    conv = [Fr(1)] + [Fr(0)] * N
    for dq in delta:
        seq = [Fr(0)] * (N + 1)
        for j in range(max(0, dq), N + 1):
            seq[j] = Fr(1, fac[j] * fac[j - dq])
        new = [Fr(0)] * (N + 1)
        for a, ca in enumerate(conv):
            if ca == 0:
                continue
            for b in range(N + 1 - a):
                if seq[b]:
                    new[a + b] += ca * seq[b]
        conv = new
    return [conv[n] * fac[n] ** 2 / Fr(f) ** (2 * n) for n in range(N + 1)]


def u_brute(f, n, delta):
    """brute force P(Y_n = delta) for tiny n (validation only)."""
    import itertools
    cnt = 0
    for a in itertools.product(range(f), repeat=n):
        for b in itertools.product(range(f), repeat=n):
            v = [0] * f
            for x in a:
                v[x] += 1
            for x in b:
                v[x] -= 1
            if tuple(v) == tuple(delta):
                cnt += 1
    return Fr(cnt, f ** (2 * n))


PI_LO = Fr(314159265358979, 10 ** 14)


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


def robbins_tail(f, A1):
    """rational >= sum_{n >= f*A1} u_n(0).  Uses u_n(0) <= Amax_{f floor(n/f)} <= sqrt(f)(2 pi a)^{-s}, s=(f-1)/2,
    and sum_{a >= A1} a^{-s} <= A1^{-s} + A1^{1-s}/(s-1) (needs s > 1, i.e. f >= 4)."""
    s2 = f - 1                      # 2s
    assert s2 > 2
    two_pi_lo = 2 * PI_LO

    def pow_neg(x, e2):             # upper bound for x^{-e2/2}, x >= 1
        v = Fr(1) / Fr(x) ** (e2 // 2)
        if e2 % 2:
            v /= sqrt_lo(x)
        return v
    pref = sqrt_hi(f) * pow_neg(two_pi_lo, s2)
    ssum = pow_neg(A1, s2) + pow_neg(A1, s2 - 2) / (Fr(s2, 2) - 1)
    return f * pref * ssum


def forward_tail(f, N0, A1=40):
    """rational T >= sum_{n > N0} u_n(0) (exact terms N0 < n < f A1, Robbins beyond), and the exact u_n(0) list."""
    K = f * A1
    assert K > N0 + 1
    u0 = u_exact_table(f, K - 1, tuple([0] * f))
    return sum(u0[N0 + 1:K]) + robbins_tail(f, A1), u0


# ------------------------------------------------------------------------------------------- lateral, floats
def lateral_table(fam, N0, zs):
    """dict z -> list of upper bounds for l_n(z), n = 0..N0 (z are orbit representatives).
    p_n = law of X_n on the cube [-R,R]^3, R = c N0, by direct convolution (round to nearest, step law rounded up);
    the support of p_{n-1} lies in |x|_1 <= c(n-1) <= R - c, so shifting by |e|_inf <= c loses no mass."""
    c = fam.c
    R = c * N0
    S = 2 * R + 1
    E = sorted(fam.elaw)
    K = [fup(fam.elaw[e]) for e in E]
    p = np.zeros((S, S, S))
    p[R, R, R] = 1.0
    tab = {tuple(z): [] for z in zs}
    for n in range(0, N0 + 1):
        if n > 0:
            q = np.zeros_like(p)
            for e, Ke in zip(E, K):
                sq, sp = [], []
                for a in range(3):
                    if e[a] >= 0:
                        sq.append(slice(e[a], S)); sp.append(slice(0, S - e[a]))
                    else:
                        sq.append(slice(0, S + e[a])); sp.append(slice(-e[a], S))
                q[tuple(sq)] += Ke * p[tuple(sp)]
            p = q
        for z in zs:
            v = autocorr(p, R, z)
            tab[tuple(z)].append(v * (1 + LAT_INFL) + 1e-300)
    assert p.sum() >= 1.0 - 1e-9
    return tab


def autocorr(p, R, z):
    """float sum_x p(x) p(x - z) over the cube (round to nearest)."""
    S = 2 * R + 1
    s1, s2 = [], []
    for a in range(3):
        if z[a] >= 0:
            s1.append(slice(z[a], S)); s2.append(slice(0, S - z[a]))
        else:
            s1.append(slice(0, S + z[a])); s2.append(slice(-z[a], S))
    return float(np.einsum('ijk,ijk->', p[tuple(s1)], p[tuple(s2)]))

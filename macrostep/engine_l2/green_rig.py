"""Rigorous upper bounds for the pair-chain Green function G((A',D'),(B,D'')) = sum_n u_n(D''-D') l_n(B-A').

Forward part: u_n(z) = P_0(Y_n = z), Y = difference of two independent oriented walks in Z^f (uniform steps),
exact rationals for n <= N0:  u_n(z) = (n!)^2 f^{-2n} [x^n] prod_i A_{z_i}(x),  A_r(x) = sum_k x^k/(k!(k-r)!).
Tail: sum_{n>N0} u_n(0) <= exact sum up to K + Robbins bound (global.md Lemma 4.2(ii) with d -> f).
u_n(z) <= u_n(0) (Cauchy-Schwarz).

Lateral part: l_n(z) = P(D^X_n = z), D^X_n = X_n - X'_n, X_n = sum of n iid word endpoints.
p_n = law of X_n by direct convolution in binary64 with the kernel rounded up and EVERY multiply and add rounded
upward (np.nextafter after round-to-nearest); l_n(z) = sum_x p_n(x) p_n(x-z) likewise with upward rounding.
Hence all tabulated values are upper bounds (no a priori error analysis needed).
Monotonicity: l_n(0) = sum_x p_n(x)^2 is nonincreasing in n (Jensen) and l_n(z) <= l_n(0) (Cauchy-Schwarz),
so the tail sum_{n>N0} u_n(z_D) l_n(z_A) <= l_{N0}(0) * T_u(N0).
"""
import math
import itertools
import numpy as np
from fractions import Fraction as Fr
from numba import njit, prange
from common import fup, add_up, mul_up, IG, INF


# ------------------------------------------------------------------ forward walk
def u_table(f, N0, zs):
    """exact u_n(z) (Fractions) for z in zs (tuples), n = 0..N0."""
    fact = [math.factorial(k) for k in range(2 * N0 + 5)]
    out = {}
    for z in zs:
        prod = [Fr(1)] + [Fr(0)] * N0
        for r in z:
            pc = [Fr(0)] * (N0 + 1)
            for k in range(max(0, r), N0 + 1):
                pc[k] = Fr(1, fact[k] * fact[k - r])
            new = [Fr(0)] * (N0 + 1)
            for i, a in enumerate(prod):
                if a == 0:
                    continue
                for j in range(N0 + 1 - i):
                    if pc[j]:
                        new[i + j] += a * pc[j]
            prod = new
        out[tuple(z)] = [prod[n] * fact[n] ** 2 / Fr(f) ** (2 * n) for n in range(N0 + 1)]
    return out


def robbins_tail(f, A):
    """rational >= sum_{n >= f*A} u_n(0)  (global.md Lemma 4.2(ii) with d -> f: u_k <= A_{f floor(k/f)}
    <= sqrt(f) (2 pi floor(k/f))^{-(f-1)/2}; sum_{a>=A} a^{-s} <= A^{-s} + A^{1-s}/(s-1), s = (f-1)/2 > 1)."""
    s2 = f - 1
    assert Fr(s2, 2) > 1
    pref = IG.sqrt_hi(f) / (IG.SQRT_2PI_LO ** s2)
    Sm_s = IG.pow_neg_half(A, s2) + IG.pow_neg_half(A, s2 - 2) / (Fr(s2, 2) - 1)
    return f * pref * Sm_s


def u0_tail(f, N0, A=40):
    """rational T >= sum_{n > N0} u_n(0): exact terms N0 < n < f*A, Robbins bound for n >= f*A.
    Also returns the exact list u_0..u_{fA-1}."""
    K = f * A
    assert K > N0 + 1
    u = u_table(f, K - 1, [tuple([0] * f)])[tuple([0] * f)]
    return sum(u[N0 + 1:K]) + robbins_tail(f, A), u


def dkey(D):
    return tuple(sorted((int(x) for x in D), reverse=True))


def D_members(Dt):
    return sorted(set(itertools.permutations(Dt)))


# ------------------------------------------------------------------ lateral walk
@njit(cache=True)
def conv_step(p, K, E, R, rad):
    """q(x) = sum_e K[e] p(x - e) on the cube [-R, R]^3 (only |x|_1 <= rad computed); round-to-nearest."""
    S = 2 * R + 1
    q = np.zeros_like(p)
    ne = E.shape[0]
    for x0 in range(-rad, rad + 1):
        for x1 in range(-(rad - abs(x0)), rad - abs(x0) + 1):
            r2 = rad - abs(x0) - abs(x1)
            for x2 in range(-r2, r2 + 1):
                s = 0.0
                for k in range(ne):
                    y0 = x0 - E[k, 0]; y1 = x1 - E[k, 1]; y2 = x2 - E[k, 2]
                    if abs(y0) <= R and abs(y1) <= R and abs(y2) <= R:
                        v = p[y0 + R, y1 + R, y2 + R]
                        if v > 0.0:
                            s = np.nextafter(s + np.nextafter(K[k] * v, np.inf), np.inf)
                q[x0 + R, x1 + R, x2 + R] = s
    return q


@njit(parallel=True, cache=True)
def autocorr(p, R, zs):
    """l(z) = sum_x p(x) p(x - z) for each z (round-to-nearest)."""
    nz = zs.shape[0]
    out = np.zeros(nz)
    S = 2 * R + 1
    for t in prange(nz):
        z0 = zs[t, 0]; z1 = zs[t, 1]; z2 = zs[t, 2]
        s = 0.0
        for i0 in range(max(0, z0), min(S, S + z0)):
            for i1 in range(max(0, z1), min(S, S + z1)):
                for i2 in range(max(0, z2), min(S, S + z2)):
                    a = p[i0, i1, i2]
                    if a > 0.0:
                        b = p[i0 - z0, i1 - z1, i2 - z2]
                        if b > 0.0:
                            s = np.nextafter(s + np.nextafter(a * b, np.inf), np.inf)
        out[t] = s
    return out


def lateral_table(fam, N0, zreps):
    """upper bounds l_n(z) for z in zreps (list of 3-tuples), n = 0..N0; returns dict z -> np.array(N0+1)."""
    assert fam.m == 3
    c = fam.c
    R = c * N0
    E = np.array(sorted(fam.elaw), np.int64)
    K = np.array([fup(fam.elaw[tuple(e)]) for e in E])
    p = np.zeros((2 * R + 1,) * 3)
    p[R, R, R] = 1.0
    zs = np.array(zreps, np.int64)
    L = np.zeros((len(zreps), N0 + 1))
    L[:, 0] = [1.0 if not any(z) else 0.0 for z in zreps]
    for n in range(1, N0 + 1):
        p = conv_step(p, K, E, R, c * n)          # upward-rounded, zero outside the support |x|_1 <= c n
        L[:, n] = autocorr(p, R, zs)              # upward-rounded
    return {tuple(z): L[k] for k, z in enumerate(zreps)}

#!/usr/bin/env python3
"""Local constants of the oriented second-moment criterion (Theorem 4.33 of the first manuscript, basic form, and
the rho_c refinement of Lemma 4.27): the global part of the second of the two implementations of Section 4.7,
written from the statements without importing the first implementation.  All arithmetic is exact rational (fractions.Fraction); no floating point enters a claim.

Differences from the first implementation (../../overlap_revealed/crosscheck/), to make this an independent check:
  * u_k = P(Z_k=0) is computed by a sum over integer PARTITIONS of k (not by polynomial powers);
  * tbar = tanh K' = (g_K-1)/(g_K+1) exactly (the first implementation uses the cruder tbar = K');
  * exponentials are bounded by Taylor polynomials with explicit remainders; sqrt and pi by rational bounds;
  * the rho_c constant uses the exact identities e^{2K'}=g_K, e^{2M}=g_K^{2d-1}/g_h (only sqrt(g_K) is bounded).
Inputs per case: d, p (rational lower bound on p_B), g_K = e^{2K'}, g_h = e^{2h} (rationals), xbar (rational,
verified here to satisfy tanh(|h| + 2dK' xbar) < xbar exactly).
"""
from fractions import Fraction as Fr
from math import comb, factorial, isqrt
import sys, time

# ---------------------------------------------------------------- rational helpers
_D = 10 ** 50
def up(x):
    """rational >= x with denominator 10^50 (outward rounding of an upper bound)."""
    x = Fr(x)
    return Fr(-((-x.numerator * _D) // x.denominator), _D)

def dn(x):
    """rational <= x with denominator 10^50 (outward rounding of a lower bound)."""
    x = Fr(x)
    return Fr((x.numerator * _D) // x.denominator, _D)

def sqrt_lo(x, digits=30):
    """rational r <= sqrt(x)."""
    x = Fr(x); S = 10 ** digits
    r = Fr(isqrt(x.numerator * S * S // x.denominator), S)
    assert r * r <= x
    return r

def sqrt_hi(x, digits=30):
    r = sqrt_lo(x, digits) + Fr(1, 10 ** digits)
    assert r * r >= x
    return r

def exp_hi(x, K=30):
    """rational upper bound for e^x, 0 <= x <= 1: Taylor sum + 3 x^{K+1}/(K+1)!."""
    x = Fr(x); assert 0 <= x <= 1
    return sum(x ** k / factorial(k) for k in range(K + 1)) + 3 * x ** (K + 1) / factorial(K + 1)

PI_LO = Fr(314159265358979, 10 ** 14)         # < pi
SQRT_2PI_LO = sqrt_lo(2 * PI_LO)              # < sqrt(2 pi)

def pow_neg_half(x, s2):
    """rational upper bound for x^{-s2/2}, x >= 1 rational, s2 a positive integer."""
    x = Fr(x)
    val = x ** (-(s2 // 2))
    if s2 % 2:
        val /= sqrt_lo(x)
    return val

# ---------------------------------------------------------------- collision Green function (partition sum)
def partitions(k, maxparts, maxpart=None):
    """yield partitions of k (non-increasing tuples of positive parts) with at most maxparts parts."""
    if maxpart is None:
        maxpart = k
    if k == 0:
        yield ()
        return
    if maxparts == 0:
        return
    for first in range(min(k, maxpart), 0, -1):
        for rest in partitions(k - first, maxparts - 1, first):
            yield (first,) + rest

def u_exact(d, k):
    """u_k = sum_{compositions} (k!/(d^k prod k_i!))^2 via partitions with multiplicities."""
    fk = factorial(k); fd = factorial(d)
    tot = 0
    for lam in partitions(k, d):
        mult = {}
        for v in lam:
            mult[v] = mult.get(v, 0) + 1
        den_c = factorial(d - len(lam))            # zeros
        for c in mult.values():
            den_c *= factorial(c)
        ncomp = fd // den_c                        # compositions that rearrange lam (padded with zeros)
        denom = 1
        for v in lam:
            denom *= factorial(v)
        multi = fk // denom                        # multinomial k!/prod k_i!
        tot += ncomp * multi * multi
    return Fr(tot, d ** (2 * k))

def green(d, A=4):
    N = d * A
    u = [u_exact(d, k) for k in range(N)]
    s2 = d - 1                                     # s = (d-1)/2 = s2/2
    pref = sqrt_hi(d) / (SQRT_2PI_LO ** s2)         # >= sqrt(d) (2 pi)^{-s}
    # sum_{a>=A} a^{-s} <= A^{-s} + A^{1-s}/(s-1);  sum_{a>=A} a^{1-s} <= A^{1-s} + A^{2-s}/(s-2)
    Sm_s = pow_neg_half(A, s2) + pow_neg_half(A, s2 - 2) / (Fr(s2, 2) - 1)
    Sm_s1 = pow_neg_half(A, s2 - 2) + pow_neg_half(A, s2 - 4) / (Fr(s2, 2) - 2)
    tail0 = d * pref * Sm_s                        # >= sum_{k>=dA} u_k
    tail1 = pref * (d * d * Sm_s1 + Fr(d * (d + 1), 2) * Sm_s)   # >= sum_{k>=dA} (k+1) u_k
    G_lo = sum(u); G_hi = G_lo + tail0
    S2_hi = sum((k + 1) * u[k] for k in range(N)) + tail1
    return dict(u=u, N=N, G_lo=G_lo, G_hi=G_hi, S2_hi=S2_hi, tail0=tail0, tail1=tail1)

# ---------------------------------------------------------------- the criterion
def Mmax_multinomial(d, m):
    q, j = divmod(m, d)
    return Fr(factorial(m), factorial(q + 1) ** j * factorial(q) ** (d - j))

def Gam(d, t, m):
    return sum(comb(d, j) / (1 - 2 * t * (d - 2 * j)) ** (m + 1) for j in range(d + 1)) / 2 ** d

def N_radius(d, r):
    return sum(comb(d, k) * comb(d - k, l) * comb(r - 1, k - 1) * comb(r - 1, l - 1)
               for k in range(1, d + 1) for l in range(1, d - k + 1))

def atanh_bounds(z, K=25):
    """rational [lo, hi] containing atanh(z) for 0 <= z < 1 (alternating-free positive series + geometric tail)."""
    z = Fr(z); assert 0 <= z < 1
    lo = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(K + 1))
    hi = lo + z ** (2 * K + 3) / ((2 * K + 3) * (1 - z * z))
    return lo, hi

def mf_ok(d, gK, gh, xbar):
    """exact-rational test of tanh(|h| + 2dK' xbar) < xbar, i.e. |h| + 2dK' xbar < atanh(xbar), with
    |h| = -ln(g_h)/2 = atanh((1-g_h)/(1+g_h)) and K' = ln(g_K)/2 = atanh((g_K-1)/(g_K+1))."""
    H_hi = atanh_bounds((1 - gh) / (1 + gh))[1]
    K_hi = atanh_bounds((gK - 1) / (gK + 1))[1]
    return H_hi + 2 * d * K_hi * xbar < atanh_bounds(xbar)[0]

def evaluate(d, p, gK, gh, xbar, R=40, green_cache={}):
    p, gK, gh, xbar = Fr(p), Fr(gK), Fr(gh), Fr(xbar)
    out = {}
    # hypotheses (33) of the manuscript:  h<0, 0<=K'<=|h|, 2dK'<1, d>=6
    assert gK >= 1 and gh < 1 and gK * gh <= 1 and gK ** (2 * d) < Fr(2718, 1000) and d >= 6
    assert mf_ok(d, gK, gh, xbar)
    rho = (1 - xbar) / 2                                   # rho_- (Lemma 4.22)
    kappa = (1 - rho) / rho
    t = (gK - 1) / (gK + 1)                                # tanh K' exactly
    alpha = 2 * d * t; assert alpha < 1
    qs = d * t / (1 - alpha); assert qs < 1
    # b(m) and psibar(r); upper bounds are rounded UP to 50 decimals (outward) to keep denominators small
    Mm = 2 * R + 6
    b_exact = [t ** m * Mmax_multinomial(d, m) * Gam(d, t, m) for m in range(Mm + 1)]
    assert all(b_exact[m + 1] <= qs * b_exact[m] for m in range(Mm))  # consistency with Lemma 4.26(iv)
    b = [up(x) for x in b_exact]
    btail = lambda S: sum(b[S:]) + b[Mm] * qs / (1 - qs)
    psi = {r: up(sum(b[2 * r - abs(s)] for s in range(-r, r + 1)) + 2 * btail(r + 1)) for r in range(1, R + 1)}
    x1 = up(kappa * psi[1]); assert x1 <= 1
    E1 = up(exp_hi(x1))
    U1 = E1 - 1
    SU = up(sum(N_radius(d, r) * kappa * psi[r] * E1 for r in range(1, R + 1)))
    # tail r>R: N_r <= C(r+d-1,d-1)^2, psibar(r) <= (2r+1+c) b(r), b(r) <= b(R) qs^{r-R}; terms decay with ratio <= rat
    c = 2 * qs / (1 - qs)
    first = Fr(comb(R + d, d - 1) ** 2) * (2 * R + 3 + c) * b[R] * qs
    rat = Fr(R + 1 + d, R + 2) ** 2 * (2 * R + 5 + c) / (2 * R + 3 + c) * qs
    assert rat < 1
    SU = up(SU + kappa * E1 * first / (1 - rat))
    if d not in green_cache:
        green_cache[d] = green(d)
    g = green_cache[d]
    eta = up(sum(min(U1, uk * SU) for uk in g["u"]) + SU * g["tail0"])
    assert eta < 1
    Sh2 = up(g["S2_hi"] / g["G_lo"] ** 2 - 1)
    T1 = up(U1 * Sh2)
    F = up(1 - 1 / g["G_hi"])
    term1 = 1 / (d * rho * p)
    term2 = (F - Fr(1, d)) / rho
    term3 = T1 / ((1 - eta) * rho)
    score = term1 + term2 + term3
    # rho_c refinement (Lemma 4.27): c_cov = 2sinh(2K')/(e^{K'}cosh 2M + e^{-K'}cosh 2M')^2, exact up to sqrt(gK)
    E2M = gK ** (2 * d - 1) / gh                           # e^{2M},  M = |h| + (2d-1)K'
    E2Mp = gK ** (2 * d - 1)                               # e^{2M'}, M' = (2d-1)K'
    ch2M = (E2M + 1 / E2M) / 2; ch2Mp = (E2Mp + 1 / E2Mp) / 2
    sK_hi, sK_lo = sqrt_hi(gK), sqrt_lo(gK)                # e^{K'} bounds
    c_cov = (gK - 1 / gK) / (sK_hi * ch2M + ch2Mp / sK_lo) ** 2
    assert rho ** 2 >= c_cov / 4
    rho_c = rho + c_cov / (4 * rho)
    score_c = 1 / (d * rho_c * p) + term2 + term3
    theta = lambda sc: Fr(99, 100) * rho * (1 - eta) * (1 - sc) / (1 - Fr(1, d)) if sc < 1 else None
    out.update(d=d, rho=rho, kappa=kappa, t=t, alpha=alpha, qs=qs, b1=b[1], psi1=psi[1], U1=U1, SU=SU, eta=eta,
               G_lo=g["G_lo"], G_hi=g["G_hi"], F=F, Sh2=Sh2, T1=T1, term1=term1, term2=term2, term3=term3,
               score=score, c_cov=c_cov, rho_c=rho_c, score_c=score_c, theta=theta(score), theta_c=theta(score_c))
    return out

CASES = [  # label, d, p, g_K, g_h  (the rows of Table 2 and the uncertified row d10 t=3/25; xbar found and verified here)
    ("d12 t=3/25 headline", 12, Fr(2161, 10000), Fr(102634, 100000), Fr(83527, 100000)),
    ("d12 t=1/8", 12, Fr(2224, 10000), Fr(1028469, 10 ** 6), Fr(407843, 500000)),
    ("d11 t=13/100", 11, Fr(2310, 10000), Fr(6443, 6250), Fr(25553, 31250)),
    ("d10 t=13/100", 10, Fr(2336, 10000), Fr(1031101, 10 ** 6), Fr(840843, 10 ** 6)),
    ("d10 t=3/25", 10, Fr(2199, 10000), Fr(1026679, 10 ** 6), Fr(437233, 500000)),
    ("d12 t=11/100", 12, Fr(2023, 10000), Fr(255577, 250000), Fr(873143, 10 ** 6)),
    ("d12 t=23/200", 12, Fr(2094, 10000), Fr(1024283, 10 ** 6), Fr(855173, 10 ** 6)),
    ("d12 t=13/100", 12, Fr(2283, 10000), Fr(515339, 500000), Fr(794291, 10 ** 6)),
    ("d11 t=3/25", 11, Fr(2180, 10000), Fr(513251, 500000), Fr(427717, 500000)),
    ("d10 t=27/200", 10, Fr(2399, 10000), Fr(516713, 500000), Fr(102809, 125000)),
]

def find_xbar(d, gK, gh):
    from math import log, tanh
    Kp = log(float(gK)) / 2; H = -log(float(gh)) / 2; x = 0.0
    for _ in range(5000):
        x = tanh(H + 2 * d * Kp * x)
    xb = Fr(int(x * 10 ** 7) + 1, 10 ** 7)
    while not mf_ok(d, gK, gh, xb):
        xb += Fr(1, 10 ** 7)
    return xb

if __name__ == "__main__":
    sel = sys.argv[1:]
    for lab, d, p, gK, gh in CASES:
        if sel and not any(s in lab for s in sel):
            continue
        t0 = time.time()
        xb = find_xbar(d, gK, gh)
        o = evaluate(d, p, gK, gh, xb)
        # DIRECTED rounding of every printed bound (lo = floor, hi = ceil, exact); logic unchanged.
        from math import floor as _fl, ceil as _ce
        lo = lambda x, k=6: "None" if x is None else f"{float(Fr(_fl(Fr(x) * 10 ** k), 10 ** k)):.{k}f}"
        hi = lambda x, k=6: "None" if x is None else f"{float(Fr(_ce(Fr(x) * 10 ** k), 10 ** k)):.{k}f}"
        print(f"{lab}: xbar={xb} rho_->={lo(o['rho'],7)} kappa<={hi(o['kappa'])} tanhK'<={hi(o['t'],7)} alpha<={hi(o['alpha'])} q*<={hi(o['qs'])}")
        print(f"   G_d in [{lo(o['G_lo'],10)}, {hi(o['G_hi'],10)}]  F_d<={hi(o['F'],10)}  sum h^2<={hi(o['Sh2'],8)}")
        print(f"   b1<={hi(o['b1'],8)} psibar1<={hi(o['psi1'],8)} U1<={hi(o['U1'])} SumU<={hi(o['SU'],4)} eta'<={hi(o['eta'])} T1<={hi(o['T1'],8)}")
        print(f"   terms {hi(o['term1'])} + {hi(o['term2'])} + {hi(o['term3'])}:  SCORE<={hi(o['score'])}  theta_*>={lo(o['theta'])}")
        print(f"   c_cov>={lo(o['c_cov'],7)} rho_c>={lo(o['rho_c'])}  SCORE_c<={hi(o['score_c'])}  theta_c>={lo(o['theta_c'])}   [{time.time()-t0:.1f}s]", flush=True)

"""Rigorous evaluation of the oriented Ising-site / Bernoulli-bond second-moment criterion (global.md, Thm G).

Inputs (rationals): d, Kp (= K' >= 0), H (= -h > 0), pB (lower bound for the bond probability).
Optionally rho_override: a rigorous lower bound for the true one-site marginal to use instead of the
mean-field bound (must be <= the mean-field value to be meaningful).

All arithmetic is exact rational (fractions.Fraction) except exp/tanh, which use mpmath interval
arithmetic (outward rounding) and are converted to rational upper/lower bounds.

Score = 1/(d rho pB) + (F_d - 1/d)/rho + T1/((1-eta) rho) ;  criterion: Score < 1 and eta < 1.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import mpmath as mp
from green import green_bounds
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

mp.mp.prec = 200
iv = mp.iv
iv.prec = 200


def _mpf_to_fr(v):
    """Exact rational value of a raw mpf tuple (sign, man, exp, bc) or an mpf."""
    if isinstance(v, tuple):
        sign, man, exp, _bc = v
    else:
        sign, man, exp, _bc = mp.mpf(v)._mpf_
    val = Fr(man) * Fr(2) ** exp if exp >= 0 else Fr(man, 2 ** (-exp))
    return -val if sign else val


def q_up(x_iv):
    """Exact rational value of the (outward-rounded) upper endpoint of an mp interval."""
    return _mpf_to_fr(x_iv._mpi_[1])


def q_lo(x_iv):
    return _mpf_to_fr(x_iv._mpi_[0])


def ivq(fr):
    return iv.mpf(fr.numerator) / iv.mpf(fr.denominator)


def exp_up(fr):
    return q_up(iv.exp(ivq(fr)))


def tanh_up(fr):
    x = ivq(fr)
    return q_up(1 - 2 / (iv.exp(2 * x) + 1))


def mean_field_mbar(d, Kp, H):
    """Rational mbar with tanh(H + 2dK' mbar) <= mbar (rigorously), found by bisection."""
    a = 2 * d * Kp
    assert a < 1 and 0 <= Kp <= H
    lo, hi = Fr(0), Fr(1)
    for _ in range(60):
        mid = (lo + hi) / 2
        if tanh_up(H + a * mid) <= mid:
            hi = mid
        else:
            lo = mid
    # round hi up to 10 decimal places, re-verify
    mbar = Fr(int((hi * 10 ** 10).__ceil__()), 10 ** 10)
    assert tanh_up(H + a * mbar) <= mbar
    return mbar


def multinomial_max(d, m):
    q, j = divmod(m, d)
    return Fr(factorial(m), factorial(q + 1) ** j * factorial(q) ** (d - j))


def Gamma(d, t, m):
    s = Fr(0)
    for j in range(d + 1):
        s += comb(d, j) / (1 - 2 * t * (d - 2 * j)) ** (m + 1)
    return s / 2 ** d


def N_radius(d, r):
    """#{z in Z^d : sum z = 0, |z|_1 = 2r}."""
    tot = 0
    for k in range(1, d + 1):
        for l in range(1, d - k + 1):
            tot += comb(d, k) * comb(d - k, l) * comb(r - 1, k - 1) * comb(r - 1, l - 1)
    return tot


def evaluate(d, Kp, H, pB, R=60, rho_override=None, verbose=True):
    Kp, H, pB = Fr(Kp), Fr(H), Fr(pB)
    t = Kp  # rigorous upper bound for tanh(K')
    alpha = 2 * d * t
    assert alpha < 1
    mbar = mean_field_mbar(d, Kp, H)
    rho_mf = (1 - mbar) / 2
    rho = Fr(rho_override) if rho_override is not None else rho_mf
    kappa = (1 - rho) / rho
    q = d * t / (1 - alpha)  # ratio bound b(m+1)/b(m) <= q
    assert q < 1
    Mmax = 2 * R + 4
    b = [t ** m * multinomial_max(d, m) * Gamma(d, t, m) for m in range(Mmax + 1)]
    # sanity: b nonincreasing (also implied analytically by the ratio bound)
    assert all(b[m + 1] <= b[m] for m in range(Mmax))

    def btail(S):  # upper bound for sum_{s>=S} b(s)
        if S <= Mmax:
            return sum(b[S:Mmax + 1]) + b[Mmax] * q / (1 - q)
        raise ValueError

    def psibar(r):
        tot = sum(b[2 * r - abs(s)] for s in range(-r, r + 1))
        return tot + 2 * btail(r + 1)

    psi = [None] + [psibar(r) for r in range(1, R + 1)]
    assert all(psi[r + 1] <= psi[r] for r in range(1, R))
    E1 = exp_up(kappa * psi[1])
    U1 = E1 - 1
    # sum over the lattice A_{d-1}\{0} of U(z) <= sum_r N_r * kappa*psi_r*e^{kappa psi_1}
    SU = sum(N_radius(d, r) * kappa * psi[r] * E1 for r in range(1, R + 1))
    # tail r > R: N_r <= C(r+d-1,d-1)^2, psi_r <= (2r+1+2q/(1-q)) b(r), b(r) <= b(R) q^{r-R}
    c = 2 * q / (1 - q)
    first = Fr(comb(R + d, d - 1) ** 2) * (2 * (R + 1) + 1 + c) * b[R] * q
    ratio = Fr(R + 1 + d, R + 2) ** 2 * (2 * R + 5 + c) / (2 * R + 3 + c) * q
    assert ratio < 1
    SU += kappa * E1 * first / (1 - ratio)
    g = green_bounds(d)
    u, N = g["u"], g["N"]
    eta = sum(min(U1, u[k] * SU) for k in range(N)) + SU * g["tail0"]
    Sh2 = g["S2_hi"] / g["G_lo"] ** 2 - 1  # >= sum_{y != 0} h(y)^2
    T1 = U1 * Sh2
    F = 1 - 1 / g["G_hi"]
    term1 = 1 / (d * rho * pB)
    term2 = (F - Fr(1, d)) / rho
    term3 = T1 / ((1 - eta) * rho) if eta < 1 else None
    score = term1 + term2 + term3
    theta = rho * (1 - eta) * (1 - score) / (1 - Fr(1, d)) if score < 1 else None
    # refinement (Sec. 4.8): shared edges cost 1/rho_c, rho_c = rho + c_cov/(4 rho),
    # c_cov = 2 sinh(2K') / (e^{K'} cosh(2M) + e^{-K'} cosh(2M'))^2, M = |h| + (2d-1)K', M' = (2d-1)K'
    Kiv, Miv, M2iv = ivq(Kp), ivq(H + (2 * d - 1) * Kp), ivq((2 * d - 1) * Kp)
    ch = lambda x: (iv.exp(x) + iv.exp(-x)) / 2
    sh = lambda x: (iv.exp(x) - iv.exp(-x)) / 2
    den = iv.exp(Kiv) * ch(2 * Miv) + iv.exp(-Kiv) * ch(2 * M2iv)
    c_cov = q_lo(2 * sh(2 * Kiv) / den ** 2)
    assert rho ** 2 >= c_cov / 4  # rho -> rho + c/(4 rho) increasing on [sqrt(c)/2, 1]
    rho_c = rho + c_cov / (4 * rho)
    term1c = 1 / (d * rho_c * pB)
    score_c = term1c + term2 + term3
    theta_c = rho * (1 - eta) * (1 - score_c) / (1 - Fr(1, d)) if score_c < 1 else None
    out = dict(d=d, Kp=Kp, H=H, pB=pB, mbar=mbar, rho_mf=rho_mf, rho=rho, kappa=kappa, t=t, alpha=alpha,
               q=q, b1=b[1], b2=b[2], psi1=psi[1], psi2=psi[2], U1=U1, SU=SU, eta=eta, Sh2=Sh2, T1=T1,
               G_lo=g["G_lo"], G_hi=g["G_hi"], F=F, term1=term1, term2=term2, term3=term3, score=score,
               theta=theta, c_cov=c_cov, rho_c=rho_c, term1c=term1c, score_c=score_c, theta_c=theta_c)
    if verbose:
        f = lambda x: f"{float(x):.6g}" if x is not None else "None"
        print(f"d={d} K'={f(Kp)} H={f(H)} pB={f(pB)}: mbar={f(mbar)} rho_MF>={f(rho_mf)} rho_used={f(rho)} "
              f"kappa={f(kappa)} alpha={f(alpha)} q={f(q)}")
        print(f"   b1={f(b[1])} b2={f(b[2])} psi1={f(psi[1])} psi2={f(psi[2])} U1={f(U1)} SumU={f(SU)} "
              f"eta'={f(eta)} sum h^2<={f(Sh2)} T1={f(T1)}")
        print(f"   G in [{f(g['G_lo'])},{f(g['G_hi'])}] F_d<={f(F)} | terms: {f(term1)} + {f(term2)} + {f(term3)} "
              f"= SCORE {f(score)}  theta_*>={f(theta)} (before e^-0.01)")
        print(f"   refined: c_cov>={f(c_cov)} rho_c>={f(rho_c)} term1c={f(term1c)} SCORE_c={f(score_c)} theta_c>={f(theta_c)}")
    return out


if __name__ == "__main__":
    import sys
    d = int(sys.argv[1]); Kp = Fr(sys.argv[2]); H = Fr(sys.argv[3]); pB = Fr(sys.argv[4])
    rho_o = Fr(sys.argv[5]) if len(sys.argv) > 5 else None
    evaluate(d, Kp, H, pB, rho_override=rho_o)

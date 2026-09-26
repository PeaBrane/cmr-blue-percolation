#!/usr/bin/env python3
"""Track B part 1 -- exact rational certificates for the local inputs (see ../local.md).

Headline case: d=12 (m=24), unit fair couplings, t=tanh(beta)=3/25.
  (1) p_B >= 2161/10000                         [signed-residual vertex certificate, 24 vertex classes]
  (2) Holley line for Ising(K',h), Gamma=e^{2K'}=102634/100000, H=e^{2h}=83527/100000,
      checked in all 25 environments S=-24..24  [decoupled tangent-plane / vertex certificate]
  (3) mean-field root < 649/5000, hence Ising(K',h) density rho > 4351/10000
  (4) exact witness values (upper bounds on the class infima) for context.
All inequalities are decided in exact rational arithmetic (fractions.Fraction). No floating point value enters a
claimed inequality; floats appear only in printed summaries.  Runtime: a few seconds.
"""
from fractions import Fraction as Fr
from math import comb, log
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

T = Fr(3, 25); D = 12; M = 2*D

def base(t):
    w = (1+t)/(1-t)                   # e^{2 beta}
    a = w*w/(1+w*w)                   # e^{2beta}/(2cosh 2beta)
    C = (w+1/w)/2                     # cosh(2 beta)
    return w, a, C

def ch2(w, z):                        # 2cosh(beta z), z even (exact)
    assert z % 2 == 0
    l = z//2
    return w**l + w**(-l)

# ------------------------------------------------------------------ (1) p_B certificate
RHO_U = {1: Fr(1080881, 10**6), 3: Fr(1042249, 10**6), 5: Fr(984864, 10**6), 7: Fr(916142, 10**6),
         9: Fr(842612, 10**6), 11: Fr(768993, 10**6), 13: Fr(698220, 10**6), 15: Fr(631872, 10**6),
         17: Fr(570633, 10**6), 19: Fr(514653, 10**6), 21: Fr(463784, 10**6), 23: Fr(417731, 10**6)}

def pB_certificate(t=T, d=D, target=Fr(2161, 10000), rho_U=RHO_U):
    m = 2*d; n = m-1
    w, a, C = base(t)
    pre = 1-1/(w*w)                                    # 1-e^{-4beta}
    r_low = target/pre
    lam = w*w*(1-r_low)/r_low                          # (a/(1-a)) (1-r_low)/r_low
    K = lambda z: ch2(w, 1+z)                          # K(z)=2cosh(beta(1+z)), z odd
    gam = []
    for j in range(n+1):
        tot = Fr(0)
        for al in range(j+1):
            for b in range(n-j+1):
                A = 2*al-j; B = 2*b-(n-j); U = A+B; z = A-B
                P = comb(j, al)*comb(n-j, b)*a**(al+b)*(1-a)**(n-al-b)
                kap = w**(-2*U)-lam
                Kz = K(z)
                if kap >= 0:
                    phi = 1/(Kz*Kz)                    # inverse Jensen branch
                else:
                    c = rho_U[U]*K(U)                  # tangent branch, c_U = rho_U K(U)
                    phi = 3/(c*c)-2*Kz/(c*c*c)
                tot += P*kap*phi
        gam.append(tot)
    return all(g <= 0 for g in gam), gam, lam, pre

def pB_aligned_witness(t=T, d=D):
    m = 2*d; n = m-1
    w, a, C = base(t); pre = 1-1/(w*w)
    X = sum(comb(n, b)*a**b*(1-a)**(n-b)/ch2(w, 1+(2*b-n))**2 for b in range(n+1))
    Y = sum(comb(n, b)*a**b*(1-a)**(n-b)/ch2(w, (2*b-n)-1)**2 for b in range(n+1))
    return pre*a*X/(a*X+(1-a)*Y)

# ------------------------------------------------------------------ (2) odds / Holley certificate
W_K = {k: Fr(0) for k in range(16)}
W_K.update({16: Fr(442, 10**4), 17: Fr(964, 10**4), 18: Fr(1519, 10**4), 19: Fr(2104, 10**4), 20: Fr(2718, 10**4),
            21: Fr(3359, 10**4), 22: Fr(4027, 10**4), 23: Fr(4719, 10**4), 24: Fr(5435, 10**4)})

def G_table(w, k, m):
    km = m-k
    return {z: sum(Fr(comb(km, l), 2**km)*ch2(w, z+2*l-km) for l in range(km+1))
            for z in range(-k-2, k+3) if (z-k) % 2 == 0}

def numerator_bound(w, a, k, m, wk):
    """min over vertex classes j of the tangent-plane lower bound; all j computed (symmetry j<->k-j is re-checked)."""
    G = G_table(w, k, m)
    c = []
    for i in range(k+1):
        S = 2*i-k
        c.append(G[0] if k == 0 else (1-wk)*G[S] + wk*(Fr(i, k)*G[S-2] + Fr(k-i, k)*G[S+2]))
    vals = []
    for j in range(k+1):
        tot = Fr(0)
        for al in range(j+1):
            for b in range(k-j+1):
                A = 2*al-j; B = 2*b-(k-j); ci = c[al+b]
                tot += comb(j, al)*comb(k-j, b)*a**(al+b)*(1-a)**(k-al-b)*(3/(ci*ci)-2*G[A-B]/(ci*ci*ci))
        vals.append(tot)
    assert all(vals[j] == vals[k-j] for j in range(k+1))
    return min(vals), vals

def denominator_exact(w, a, k, m):
    """max over vertex classes j of E_{lambda_-}[(2cosh beta Z)^{-2}]; returns max, all values."""
    km = m-k; vals = []
    for j in range(km+1):
        P = [Fr(1)]
        def add(P, p):
            Q = [Fr(0)]*(len(P)+1)
            for i, x in enumerate(P):
                Q[i] += x*(1-p); Q[i+1] += x*p
            return Q
        for _ in range(j): P = add(P, a)
        for _ in range(km-j): P = add(P, 1-a)
        for _ in range(k): P = add(P, Fr(1, 2))
        vals.append(sum(x/ch2(w, 2*i-m)**2 for i, x in enumerate(P)))
    return max(vals), vals

def odds_certificate(t=T, d=D, Gam=Fr(102634, 100000), H=Fr(83527, 100000), w_k=W_K, verbose=True):
    m = 2*d
    w, a, C = base(t)
    rows = []; ok = True
    for k in range(m+1):
        N, _ = numerator_bound(w, a, k, m, w_k[k])
        Dk, dv = denominator_exact(w, a, k, m)
        bal = dv[(m-k)//2]
        assert Dk == bal                               # balanced-maximum lemma (proved analytically), re-checked
        S = 2*k-m
        lower = C**S*N/Dk
        good = lower >= H*Gam**S
        ok &= good
        rows.append((k, S, N, Dk, lower, good))
        if verbose:
            print(f"  k+={k:2d} S={S:+3d} w_k={float(w_k[k]):.4f}  N_k={float(N):.9f}  D_k={float(Dk):.9f}  "
                  f"log(C^S N/D)={log(lower):+.6f}  2h+2K'S={log(H*Gam**S):+.6f}  slack={log(lower)-log(H*Gam**S):+.6f}  "
                  f"{'ok' if good else 'FAIL'}")
    return ok, rows

# ------------------------------------------------------------------ (3) mean-field density
def density_certificate(m=M, Gam=Fr(102634, 100000), H=Fr(83527, 100000), xbar=Fr(649, 5000)):
    # side conditions: K'<=|h| <=> Gam*H<=1 ; m K'<1 <=> Gam^m < e^2, and e^2>7
    side = (Gam*H <= 1) and (Gam**m < 7) and (H < 1)
    zx = m*xbar; p, q = zx.numerator, zx.denominator
    root_ok = (1/H)**q*Gam**p < ((1+xbar)/(1-xbar))**q   # tanh(|h|+m K' xbar) < xbar
    return side and root_ok, (1-xbar)/2

if __name__ == "__main__":
    print("== (1) p_B certificate, d=12, t=3/25")
    ok1, gam, lam, pre = pB_certificate()
    print(f"  lambda={float(lam):.10f}; Gamma_j:", " ".join(f"{float(g):+.2e}" for g in gam))
    print(f"  max Gamma_j = {float(max(gam)):+.4e} -> p_B >= 0.2161 {'CERTIFIED' if ok1 else 'FAILED'}")
    print(f"  aligned frozen witness (exact) p_B <= {float(pB_aligned_witness()):.7f}; crude closed form (1-e^-4b)/2 = {float(pre/2):.7f}")
    print("== (2) Holley line, all 25 environments")
    ok2, rows = odds_certificate()
    print(f"  -> {'CERTIFIED' if ok2 else 'FAILED'}  (K'=(1/2)log(102634/100000)=0.0129995..., h=(1/2)log(83527/100000)=-0.0900001...)")
    print("== (3) density")
    ok3, rho = density_certificate()
    print(f"  mean-field root < 649/5000 and side conditions: {'CERTIFIED' if ok3 else 'FAILED'}; rho > {float(rho)}")
    print("ALL CERTIFIED" if (ok1 and ok2 and ok3) else "SOMETHING FAILED")

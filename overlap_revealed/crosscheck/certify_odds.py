"""Exact rational certificate for the overlap-field odds bound (manuscript, Section 4.3) and the
mean-field density bound (Lemma 4.22).

For k=k_+ plus-neighbours (k_-=m-k), S=2k-m, the q-field single-site odds satisfy (Proposition 4.9 and Lemmas 4.10, 4.11)
   Phi_+/Phi_- >= C^S * N_k / D_k ,   C=cosh(2beta)=(w+1/w)/2,
 N_k = min_j Xt_k(j)  (tangent-plane lower bound for inf_nu E_{lambda_+}[R_nu^-2]),
 D_k = max_j Y_k(j)   (exact sup_nu E_{lambda_-}[R_nu^-2], attained at frozen patterns),
with the tangent point c_i (i = # of + among the k tilted signs)
   c_i = (1-w_k) G(S_i) + w_k [ (i/k) G(S_i-2) + ((k-i)/k) G(S_i+2) ],  S_i=2i-k,
   G(z) = E_{F ~ sum of k_- fair signs} 2cosh(beta(z+F))   (= cosh(beta)^{k_-} 2cosh(beta z)),
   Xt_k(j) = sum_{alpha,b} C(j,alpha)C(k-j,b) a^{alpha+b}(1-a)^{k-alpha-b} [3/c^2 - 2 G(A-B)/c^3],  c=c_{alpha+b},
   Y_k(j)  = E[ 1/(4 cosh^2(beta Z)) ],  Z = (j a-biased) + (k_- - j (1-a)-biased) + (k fair) signs.
The Holley condition for Ising(K',h) with Gamma=e^{2K'}, H=e^{2h} rational is  C^S N_k/D_k >= H Gamma^S  for all k.
Everything is exact (fractions.Fraction); floats are used only to choose the mixing weights w_k.
"""
from fractions import Fraction as Fr
from math import comb, atanh, log, cosh, exp
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

def make(d, t):
    m = 2*d; w = (1+t)/(1-t); a = w*w/(1+w*w); C = (w+1/w)/2
    def ch2(z):  # 2cosh(beta z) for even z, exact
        assert z % 2 == 0
        l = z//2
        return w**l + w**(-l)
    return m, w, a, C, ch2

def G_table(k, m, ch2):
    km = m-k
    tab = {}
    for z in range(-k-2, k+3):
        if (z-k) % 2: continue
        tab[z] = sum(Fr(comb(km, l), 2**km)*ch2(z+2*l-km) for l in range(km+1))
    return tab

def numerator(k, m, a, ch2, wk):
    G = G_table(k, m, ch2)
    c = []
    for i in range(k+1):
        S = 2*i-k
        if k == 0:
            c.append(G[0]); continue
        c.append((1-wk)*G[S] + wk*(Fr(i, k)*G[S-2] + Fr(k-i, k)*G[S+2]))
    pa = [a**i for i in range(k+1)]; qa = [(1-a)**i for i in range(k+1)]
    vals = []
    for j in range(k//2+1):
        tot = Fr(0)
        for al in range(j+1):
            for b in range(k-j+1):
                A = 2*al-j; B = 2*b-(k-j); ci = c[al+b]
                P = comb(j, al)*comb(k-j, b)*pa[al+b]*qa[k-al-b]
                tot += P*(3/(ci*ci) - 2*G[A-B]/(ci*ci*ci))
        vals.append(tot)
    return min(vals), vals

def denominator(k, m, a, ch2):
    km = m-k
    vals = []
    for j in range(km+1):
        # distribution of Z as polynomial coefficients over values -m..m (step 2)
        P = [Fr(1)]
        def conv(P, p):  # add one sign with P(+)=p
            Q = [Fr(0)]*(len(P)+1)
            for i, x in enumerate(P):
                Q[i] += x*(1-p); Q[i+1] += x*p
            return Q
        for _ in range(j): P = conv(P, a)
        for _ in range(km-j): P = conv(P, 1-a)
        for _ in range(k): P = conv(P, Fr(1, 2))
        tot = Fr(0)
        for i, x in enumerate(P):
            Z = 2*i-m; v = ch2(Z)
            tot += x/(v*v)
        vals.append(tot)
    return max(vals), vals

# float helper to choose w_k (mixing weight of the one-off layer in the tangent point)
def choose_wk(d, tf, k):
    import numpy as np
    from scipy.optimize import minimize_scalar
    sys.path.insert(0, __file__.rsplit('/', 1)[0])
    from odds_float import num_tangent
    beta = atanh(tf); m = 2*d
    if k < 2: return 0.0
    f = lambda x: -num_tangent(beta, m, k, x)[0]
    res = minimize_scalar(f, bounds=(0, 0.9), method='bounded', options={'xatol': 1e-10})
    return res.x if -res.fun > -f(0.0) else 0.0

def run(d, t, Gam, H, xbar, verbose=True):
    m, w, a, C, ch2 = make(d, t)
    rows = []; ok = True
    for k in range(m+1):
        wk = Fr(round(choose_wk(d, float(t), k)*10**4), 10**4)
        N, nv = numerator(k, m, a, ch2, wk)
        D, dv = denominator(k, m, a, ch2)
        jbal = (m-k)//2
        assert D == dv[jbal], "balanced-maximum lemma violated?"
        S = 2*k-m
        lhs = C**S*N/D; rhs = H*Gam**S
        good = lhs >= rhs
        ok &= good
        rows.append((k, S, wk, N, D, lhs, rhs, good))
        if verbose:
            print(f"k+={k:2d} S={S:+3d} w_k={float(wk):.4f} N={float(N):.9f} D={float(D):.9f} "
                  f"log(C^S N/D)={log(lhs):+.6f} line={log(rhs):+.6f} slack={log(lhs)-log(rhs):+.6f} {'OK' if good else 'FAIL'}")
    # mean-field density bound: need H^{-q} Gam^{p} < ((1+xbar)/(1-xbar))^q with z*xbar=p/q, z=m
    zx = m*xbar; p, q = zx.numerator, zx.denominator
    lhs = (1/H)**q*Gam**p; rhs = ((1+xbar)/(1-xbar))**q
    mf_ok = lhs < rhs
    Kp = log(Gam)/2; h = log(H)/2
    cond = (Kp <= abs(h)) and (m*Kp < 1)
    if verbose:
        print(f"Holley line: Gamma={Gam} (K'={Kp:.7f}), H={H} (h={h:.7f}); all 25 environments: {'CERTIFIED' if ok else 'FAILED'}")
        print(f"mean-field root < xbar={xbar} ({float(xbar)}): {'CERTIFIED' if mf_ok else 'FAILED'}; K'<=|h| and zK'<1: {cond}")
        print(f"=> Ising(K',h) density rho > (1-xbar)/2 = {float((1-xbar)/2):.6f}")
    return ok and mf_ok and cond, rows

if __name__ == "__main__":
    d = 12; t = Fr(3, 25)
    run(d, t, Gam=Fr(102634, 100000), H=Fr(83527, 100000), xbar=Fr(649, 5000))

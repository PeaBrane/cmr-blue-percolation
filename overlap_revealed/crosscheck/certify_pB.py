"""Exact rational certificate for the q-revealed bond floor p_B (track B, local input L-b).

Setting: unit fair couplings, degree m=2d, t=tanh(beta) rational, w=e^{2beta}=(1+t)/(1-t), a=w^2/(1+w^2), n=m-1.
After the partition-extremality lemma only k'=m-1 (all m star signs a-biased) remains.  For a target r_low define
lam=w^2(1-r_low)/r_low.  Signed residual with tangent points c_U>0 (U odd, used where w^{-2U}<lam):
   Gamma_j = sum_{alpha,b} C(j,alpha)C(n-j,b) a^{alpha+b}(1-a)^{n-alpha-b} (w^{-2U}-lam) phi_U(z),
   A=2alpha-j, B=2b-(n-j), U=A+B, z=A-B, K(z)=2cosh(beta(1+z)) = w^{(1+z)/2}+w^{-(1+z)/2},
   phi_U(z)=K(z)^{-2} if w^{-2U}>=lam, else 3/c_U^2-2K(z)/c_U^3.
If Gamma_j<=0 for j=0..n then every flip-invariant cavity law gives mu_{R^-2}(s_u=+) >= r_low, hence
p_B >= (1-w^{-2}) r_low.  All arithmetic below is exact (fractions.Fraction); floats only choose c_U.
"""
from fractions import Fraction as Fr
from math import comb, cosh, atanh, log
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

def certify(d, t, pB_target, kappa=None, lamexp=None, digits=10, verbose=True, cU=None):
    # Repository copy: explicit exact tangent points cU (frozen in ../params.json) may be passed
    # instead of (kappa, lamexp); the certificate below is unchanged.
    m = 2*d; n = m-1
    w = (1+t)/(1-t); a = w*w/(1+w*w); pre = 1-1/(w*w)
    r_low = pB_target/pre
    lam = w*w*(1-r_low)/r_low
    def Kf(z):  # exact K(z)=w^{l}+w^{-l}, l=(1+z)/2 integer
        l = (1+z)//2
        return w**l + w**(-l)
    if cU is None:
        # tangent points (floats rounded to rationals: any positive value is admissible)
        beta = atanh(float(t))
        K0 = 2*cosh(beta)**m
        cU = {}
        for U in range(-n, n+1, 2):
            if w**(-2*U) < lam:
                val = kappa*(2*cosh(beta*(1+U)))**lamexp*K0**(1-lamexp)
                cU[U] = Fr(round(val*10**digits), 10**digits)
    pa = [a**i for i in range(n+1)]; qa = [(1-a)**i for i in range(n+1)]
    Kcache = {z: Kf(z) for z in range(-n, n+1, 2)}
    gam = []
    for j in range(n+1):
        tot = Fr(0)
        for al in range(j+1):
            for b in range(n-j+1):
                A = 2*al-j; B = 2*b-(n-j); U = A+B; z = A-B
                P = comb(j, al)*comb(n-j, b)*pa[al+b]*qa[n-al-b]
                kap = w**(-2*U)-lam
                Kz = Kcache[z]
                if kap >= 0:
                    phi = 1/(Kz*Kz)
                else:
                    c = cU[U]; phi = 3/(c*c) - 2*Kz/(c*c*c)
                tot += P*kap*phi
        gam.append(tot)
    ok = all(g <= 0 for g in gam)
    if verbose:
        print(f"d={d} t={t} ({float(t)}): target p_B >= {pB_target} = {float(pB_target):.6f}; r_low={float(r_low):.8f}; lambda={float(lam):.8f}")
        print("  tangent points c_U (U: c_U/K(U)):", {U: round(float(c/Kcache[U]),6) for U, c in cU.items()})
        print("  Gamma_j (float display of exact rationals):", " ".join(f"{float(g):+.3e}" for g in gam))
        print("  max_j Gamma_j =", f"{float(max(gam)):+.6e}", "->", "CERTIFIED" if ok else "FAILED")
    return ok, gam

def aligned_witness(d, t):
    """exact frozen-aligned value (upper bound on the class infimum): r_al = a X/(a X+(1-a)Y)."""
    m = 2*d; n = m-1
    w = (1+t)/(1-t); a = w*w/(1+w*w); pre = 1-1/(w*w)
    X = Fr(0); Y = Fr(0)
    for b in range(n+1):
        U = 2*b-n; P = comb(n, b)*a**b*(1-a)**(n-b)
        lp = (1+U)//2; lm = (U-1)//2
        X += P/(w**lp+w**(-lp))**2; Y += P/(w**lm+w**(-lm))**2
    r = a*X/(a*X+(1-a)*Y)
    return pre*r

if __name__ == "__main__":
    cases = [  # (d, t, certified target, kappa, lambda) ; kappa/lambda from the float optimisation (pB_family.py)
        (12, Fr(3, 25),   Fr(2160, 10000), 1.01439, 0.56291),
        (12, Fr(11, 100), Fr(2022, 10000), 1.01092, 0.60797),
        (12, Fr(13, 100), Fr(2282, 10000), 1.01842, 0.52692),
        (11, Fr(13, 100), Fr(2309, 10000), 1.01635, 0.56865),
    ]
    for d, t, tgt, kap, lx in cases:
        ok, _ = certify(d, t, tgt, kap, lx)
        wal = aligned_witness(d, t)
        w = (1+t)/(1-t)
        print(f"  aligned frozen witness (exact): p_B <= {float(wal):.7f};  tree value tanh(2beta) = {float((w*w-1)/(w*w+1)):.7f}")

"""Per-configuration EXACT check of both kernels (the a priori error program in ../engine and the directed-rounding
program in ../engine_l2). Needs NumPy and Numba. The recorded runs used
  9 623/2500 1036137/1000000 83033/100000 676921/5000000 3 1/10 40 6 201
  8 2707/10000 1039407/1000000 860671/1000000 539947/5000000 3 1/10 40 6 202
  7 1401/5000 104253/100000 175269/200000 463903/5000000 3 3/25 40 6 203]

For random states z = (A, D) and random word pairs (w, w'), this script evaluates, for all (i, i'), the defining
formula of the per-macrostep weight upper bound (Proposition 5.8 and Lemma 5.7 of the first manuscript)
    W = rho^{-N_V} beta_e^{N_E + fwd} exp(kappa B),
in exact rational arithmetic: B from the exact dyadic table values of ind_tables (integer B[r] / 2^SC, i.e. exact
upper bounds of b(r)), and exp(kappa B) by a rational Taylor upper bound.  It then runs each kernel on a
one-pair "family" (only w, w' carry mass) and checks that the kernel's per-(i, i') excess is >= the exact excess
(both kernels return UPPER bounds, so the ratio kernel/exact must be >= 1 and close to 1).
Usage: python check_configs_exact.py d p gK gh xbar c y nstates npairs seed
"""
import sys, os, math, random

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'engine'))
from ind_tables import BoostTables, fup, SC
from ind_family import Family
from ind_kernel import state_kernel, NTAY
L2DIR = os.path.join(HERE, '..', 'engine_l2')
os.environ['IG_DIR'] = L2DIR
sys.path.insert(1, L2DIR)
import indep_global as IG
import family as L2f, enum_rig as L2e, common as L2c


def exp_up_exact(x):
    return sum(x ** k / math.factorial(k) for k in range(40)) + 3 * x ** 40 / math.factorial(40)


def exact_excess(fam, T, f, c, A, D, w, wp, i, ip, rho, beta_e, kap):
    """exact rational upper bound W - 1 (formula of Proposition 5.8), with exact table values."""
    Bi = T.B
    bq = lambda r: Fr(Bi[r], 1 << SC)
    D = np.array(D); A = np.array(A)
    D1 = int(np.abs(D).sum())
    P1 = [A + fam.PTS[w, k] for k in range(fam.LEN[w] + 1)]
    P2 = [fam.PTS[wp, k] for k in range(fam.LEN[wp] + 1)]
    m1 = {}
    if D1 == 0:
        for k, u in enumerate(P1):
            for kp, v in enumerate(P2):
                if np.array_equal(u, v):
                    m1[k] = kp
    NV = len(m1)
    NE = sum(1 for k in range(1, len(P1)) if k in m1 and k - 1 in m1 and abs(m1[k] - m1[k - 1]) == 1)
    fwd = int(D1 == 0 and not np.any(A + fam.END[w] - fam.END[wp]) and i == ip)
    U1 = [k for k in range(len(P1)) if k not in m1]
    U2 = [kp for kp in range(len(P2)) if kp not in set(m1.values())]
    e = np.eye(f, dtype=np.int64)
    F1a = int(np.abs(D - e[ip]).sum()); F2a = min(int(np.abs(D - e[ip] - e[a]).sum()) for a in range(f))
    F1b = int(np.abs(D + e[i]).sum()); F2b = min(int(np.abs(D + e[i] + e[a]).sum()) for a in range(f))
    phi2 = lambda r, F: Fr(T.phi2_int(r, F), 1 << SC)

    def psi3(r0):
        tot = Fr(0)
        S = 260
        for s in range(3, S + 1):
            tot += phi2(max(0, r0 - c * (s - 1)), max(s, D1 - s))
        q = T.q
        g = 1 + 2 * q / (1 - q)
        return tot + g * bq(S + 1) / (1 - q)
    B = Fr(0)
    for k in U1:
        for kp in U2:
            B += bq(D1 + int(np.abs(P1[k] - P2[kp]).sum()))
    for k in U1:
        r1 = int(np.abs(P1[k] - fam.END[wp]).sum()); r0 = int(np.abs(P1[k]).sum())
        B += sum(bq(F1a + abs(r1 - kk)) for kk in range(c + 1)) + phi2(max(0, r1 - c), F2a) + psi3(r0)
    for kp in U2:
        r1 = int(np.abs(P2[kp] - A - fam.END[w]).sum()); r0 = int(np.abs(P2[kp] - A).sum())
        B += sum(bq(F1b + abs(r1 - kk)) for kk in range(c + 1)) + phi2(max(0, r1 - c), F2b) + psi3(r0)
    x = kap * B
    assert x <= 1
    W = (1 / rho) ** NV * beta_e ** (NE + fwd) * exp_up_exact(x)
    return W - 1, NV, NE, fwd


def main():
    a = sys.argv[1:]
    d = int(a[0]); p, gK, gh, xbar = Fr(a[1]), Fr(a[2]), Fr(a[3]), Fr(a[4]); c = int(a[5]); y = Fr(a[6])
    nst, npairs, seed = int(a[7]), int(a[8]), int(a[9])
    f = d - 3
    rng = random.Random(seed)
    o = IG.evaluate(d, p, gK, gh, xbar)
    rho, rho_c, kap, t = o['rho'], o['rho_c'], o['kappa'], o['t']
    beta_e = rho / (rho_c * p)
    fam = Family(f, c, y); T = BoostTables(d, t, c)
    coef = np.array([fup(Fr(1, math.factorial(k))) for k in range(NTAY + 1)])
    pa = np.array([fup((1 / rho) ** n) for n in range(2 * c + 4)])
    pb = np.array([fup(beta_e ** n) for n in range(2 * c + 4)])
    kapf = fup(kap)
    # inputs of the directed-rounding kernel
    Lf = L2f.Family(d, f, c, y); Lt = L2f.Tables(d, t, T.q, c)
    KE = 16
    lcoef = np.array([L2c.fup(Fr(1, math.factorial(k))) for k in range(KE + 1)]); rem = L2c.fup(Fr(3, math.factorial(KE + 1)))
    lpa = np.array([L2c.fup((1 / rho) ** n) for n in range(2 * c + 4)])
    lpb = np.array([L2c.fup(beta_e ** n) for n in range(2 * c + 4)])
    # map word index of ../engine -> word index of ../engine_l2 (same words, possibly different order)
    key = lambda pts, ln: tuple(tuple(int(v) for v in pts[k]) for k in range(ln + 1))
    L2idx = {key(Lf.PT[k], Lf.LN[k]): k for k in range(Lf.nW)}
    worst = [np.inf, 0.0, np.inf, 0.0]
    nchk = 0
    for s in range(nst):
        if s % 2 == 0:
            A = np.array([rng.randint(-4, 4) for _ in range(3)], np.int64); D = np.zeros(f, np.int64)
        else:
            A = np.array([rng.randint(-6, 6) for _ in range(3)], np.int64)
            D = np.zeros(f, np.int64)
            for _ in range(rng.randint(1, 3)):
                D[rng.randrange(f)] += 1; D[rng.randrange(f)] -= 1
        for _ in range(npairs):
            # bias towards intersecting pairs at D = 0
            w = rng.randrange(len(fam.W)); wp = rng.randrange(len(fam.W))
            if s % 2 == 0 and rng.random() < 0.7:        # force a shared point: A = P_k'(w') - P_k(w)
                A = (fam.PTS[wp, rng.randint(0, fam.LEN[wp])] - fam.PTS[w, rng.randint(0, fam.LEN[w])]).astype(np.int64)
            PM = np.zeros_like(fam.PAIRMASS); PM[w, wp] = 1.0
            acc, _ = state_kernel(A, D, f, c, fam.LEN, fam.PTS, fam.END, PM, fam.DAIDX, len(fam.dAs), T.bt, T.T1,
                                  T.PHI2, T.PSI3, pa, pb, kapf, coef)
            di = fam.DAIDX[w, wp]
            # directed-rounding kernel on the 2-word family {w, wp}
            lw = L2idx[key(fam.PTS[w], fam.LEN[w])]; lwp = L2idx[key(fam.PTS[wp], fam.LEN[wp])]
            sel = [lw] if lw == lwp else [lw, lwp]
            PT = Lf.PT[sel]; LN = Lf.LN[sel]; EN = Lf.EN[sel]
            mass = np.ones((c + 1, c + 1))
            nW2 = len(sel)
            PIDX = np.arange(nW2 * nW2, dtype=np.int64).reshape(nW2, nW2)
            accL, _ = L2e.enum_state(PT, LN, EN, mass, A, D, f, lpa, lpb, L2c.fup(kap), Lt.bt, Lt.Phi2, Lt.psi3, c,
                                     PIDX, nW2 * nW2, lcoef, rem)
            pidx = PIDX[0, 1] if lw != lwp else PIDX[0, 0]
            for i in range(f):
                for ip in range(f):
                    ex, NV, NE, fw = exact_excess(fam, T, f, c, A, D, w, wp, i, ip, rho, beta_e, kap)
                    kR = acc[di, i, ip]; kL = accL[pidx, i, ip]
                    if ex > 0:
                        r1 = Fr(kR) / ex; r2 = Fr(kL) / ex
                        assert r1 >= 1 and r2 >= 1, (A, D, w, wp, i, ip, float(r1), float(r2))
                        worst[0] = min(worst[0], float(r1)); worst[1] = max(worst[1], float(r1))
                        worst[2] = min(worst[2], float(r2)); worst[3] = max(worst[3], float(r2))
                    nchk += 1
    print(f"exact per-configuration check d={d}: {nchk} configurations; a priori error kernel/exact in [{worst[0]:.12f}, {worst[1]:.12f}],"
          f" directed-rounding kernel/exact in [{worst[2]:.12f}, {worst[3]:.12f}]  (all >= 1: both kernels are upper bounds)")


if __name__ == "__main__":
    main()

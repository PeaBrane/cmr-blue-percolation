"""Macrostep route, the a priori error program with the TIGHTENED forward tails of tail_tight.py.
The ONLY changes relative to ind_certify.py are: T_u(N0), T_u(NF) and the u_n(0) list used in eta_out come
from tail_tight.TightTails (exact u_n(0) for n <= NEX, analytic tail beyond) instead of ind_green.forward_tail
(exact u_n(0) for n < f*40, atom bound beyond).  Everything else (kernel, lateral laws, far bands, exact check)
is unchanged.  Env TT_DIR: directory holding tt_f{f}.pkl built by tail_tight.py.

Original docstring:

Usage:  python ind_certify.py d p gK gh xbar f c y RA RD N0 out.pkl [--lam LAM]
  (p, gK, gh, xbar, y exact rationals, e.g. 623/2500).

Pipeline (every step independent of ../engine_l2/; only the local constants rho_-, rho_c, kappa, t are taken from
indep_global.evaluate):
  1. exact boost tables (ind_tables), family and states (ind_family);
  2. near kernel Qbar(z, .) and etabar(z) for every orbit representative z of C' (ind_kernel, a-priori bounds);
  3. Green bounds Gbar(z'', O') (ind_green: exact forward u_n, float lateral l_n with a-priori bound, tails);
  4. Mbar = Qbar Gbar;
  5. far region: upsilon bounds and etaF (forward Khas'minskii step, strong Markov);
  6. test vector w = (lam0 - Mbar)^{-1}(etabar + delta) and an EXACT rational check of (i), (ii).
"""
import sys, os, time, math, pickle, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("ENGINE_DIR", HERE))   # the modules of ind_certify.py live next to this file
sys.path.insert(0, HERE)
from ind_tables import BoostTables, fup, check_monotone
from ind_family import Family, akey, dkey, a_orbit, d_orbit, a_types, d_types
from ind_kernel import state_kernel, NTAY, INFL
from ind_green import u_exact_table, lateral_table
from tail_tight import TightTails
TT_DIR = os.environ.get("TT_DIR", HERE)
_TT = {}


def get_tt(f):
    if f not in _TT:
        _TT[f] = TightTails(os.path.join(TT_DIR, f"tt_f{f}.pkl"))
    return _TT[f]

IG_DIR = os.environ.get("IG_DIR", HERE)          # indep_global.py lives next to this file
sys.path.insert(1, IG_DIR)
import indep_global as IG     # local constants only

_G = {}


def exp_up(x):
    """rational upper bound for e^x, 0 <= x <= 1 (Taylor to degree 30 plus 3 x^31/31!), returned as a float >= it."""
    x = Fr(x)
    assert 0 <= x <= 1
    return fup(sum(x ** k / math.factorial(k) for k in range(31)) + 3 * x ** 31 / math.factorial(31))


def _job(args):
    A, D = args
    g = _G
    acc, xm = state_kernel(np.array(A, np.int64), np.array(D, np.int64), g['f'], g['c'], g['LEN'], g['PTS'],
                           g['END'], g['PAIRMASS'], g['DAIDX'], g['nDA'], g['bt'], g['T1'], g['PHI2'], g['PSI3'],
                           g['pa'], g['pb'], g['kap'], g['coef'])
    return acc, xm


def next_state_masses(A, D, acc, dAs, f):
    """aggregate kernel bins by the ORBIT KEY of the next state; returns dict key -> float (with INFL)."""
    out = defaultdict(float)
    for kk, dA in enumerate(dAs):
        An = akey([A[q] + dA[q] for q in range(3)])
        for i in range(f):
            for ip in range(f):
                v = acc[kk, i, ip]
                if v == 0.0:
                    continue
                Dn = list(D)
                Dn[i] += 1
                Dn[ip] -= 1
                out[(An, dkey(Dn))] += v
    return {k: v * (1 + INFL) for k, v in out.items()}


def main():
    a = sys.argv[1:]
    d = int(a[0]); p, gK, gh, xbar = Fr(a[1]), Fr(a[2]), Fr(a[3]), Fr(a[4])
    f, c, y, RA, RD, N0, out = int(a[5]), int(a[6]), Fr(a[7]), int(a[8]), int(a[9]), int(a[10]), a[11]
    nproc = int(os.environ.get("NPROC", "8"))
    RS = int(os.environ.get("RS", "4"))            # forward shells |D|_1 <= 2 RS treated exactly in etaF
    R3 = int(os.environ.get("R3", "6"))            # exact lateral band for the D-types outside the C' range
    NF = int(os.environ.get("NF", "60"))
    assert d == f + 3
    T0 = time.time()
    # ---------------------------------------------------------------- 1. constants and tables
    o = IG.evaluate(d, p, gK, gh, xbar)
    rho, rho_c, kap, t, alpha = o['rho'], o['rho_c'], o['kappa'], o['t'], o['alpha']
    beta_e = rho / (rho_c * p)
    assert beta_e > 1 and rho < 1                     # W >= 1
    fam = Family(f, c, y)
    T = BoostTables(d, t, c)
    check_monotone(T)
    coef = np.array([fup(Fr(1, math.factorial(k))) for k in range(NTAY + 1)])
    pa = np.array([fup((1 / rho) ** n) for n in range(2 * c + 4)])
    pb = np.array([fup(beta_e ** n) for n in range(2 * c + 4)])
    kapf = fup(kap)
    _G.update(f=f, c=c, LEN=fam.LEN, PTS=fam.PTS, END=fam.END, PAIRMASS=fam.PAIRMASS, DAIDX=fam.DAIDX,
              nDA=len(fam.dAs), bt=T.bt, T1=T.T1, PHI2=T.PHI2, PSI3=T.PSI3, pa=pa, pb=pb, kap=kapf, coef=coef)
    print(f"[engine-t] d={d} f={f} c={c} y={y} RA={RA} RD={RD} N0={N0}: rho_->={float(rho):.7f} rho_c>={float(rho_c):.7f} "
          f"kappa<={float(kap):.6f} t={float(t):.7f} alpha={float(alpha):.5f} q*={float(T.q):.5f} "
          f"#words={len(fam.W)} Z={fam.Z} [{time.time()-T0:.0f}s]", flush=True)
    # ---------------------------------------------------------------- 2. near kernel
    Dts = d_types(f, 2 * RD)
    Ats = a_types(RA)
    states = [(Ak, Dt) for Dt in Dts for Ak in Ats]
    nS = len(states)
    sidx = {s: k for k, s in enumerate(states)}
    import multiprocessing as mp
    with mp.get_context("fork").Pool(nproc) as pool:
        res = pool.map(_job, states)
    Qbar = []
    eta = np.zeros(nS)
    xmax = 0.0
    for k, ((Ak, Dt), (acc, xm)) in enumerate(zip(states, res)):
        xmax = max(xmax, xm)
        ns = next_state_masses(Ak, Dt, acc, fam.dAs, f)
        Qbar.append(ns)
        eta[k] = sum(ns.values()) * (1 + INFL)
    assert xmax <= 1.0
    print(f"[engine-t] near kernel: {nS} states, max kappa*B = {xmax:.4f}, eta(0) = {eta[0]:.6f} [{time.time()-T0:.0f}s]",
          flush=True)
    # ---------------------------------------------------------------- 3. Green bounds
    nexts = sorted({k for ns in Qbar for k in ns} | set(states))
    Amem = {Ak: a_orbit(Ak) for Ak in Ats}
    Dmem = {Dt: d_orbit(Dt) for Dt in Dts}
    zA = set()
    for (An, Dn) in nexts:
        for Ak in Ats:
            for B in Amem[Ak]:
                zA.add(akey([B[q] - An[q] for q in range(3)]))
    zA = sorted(zA)
    zD = set()
    for (An, Dn) in nexts:
        for Dt in Dts:
            for Dm in Dmem[Dt]:
                zD.add(dkey([Dm[q] - Dn[q] for q in range(f)]))
    zD = sorted(zD)
    lat = lateral_table(fam, N0, zA)
    uD = {z: u_exact_table(f, N0, z) for z in zD}
    TT = get_tt(f)
    Tu = TT.T_u(N0)
    l0N = lat[(0, 0, 0)][N0]
    tail_unit = fup(Fr(l0N) * Tu) * (1 + 1e-12)
    print(f"[engine-t] Green: #next reps={len(nexts)} #zA={len(zA)} #zD={len(zD)} l_N0(0)<={l0N:.4e} T_u<={float(Tu):.4e} "
          f"(tight tails: NEX={TT.NEX}, tail beyond NEX <= {float(TT.tail):.4e}) [{time.time()-T0:.0f}s]", flush=True)
    LA = {}
    UD = {}

    def Gbar_row(An, Dn):
        row = np.zeros(nS)
        for j, (Ak, Dt) in enumerate(states):
            if (An, Ak) not in LA:
                v = np.zeros(N0 + 1)
                for B in Amem[Ak]:
                    v += np.array(lat[akey([B[q] - An[q] for q in range(3)])])
                LA[(An, Ak)] = v
            if (Dn, Dt) not in UD:
                UD[(Dn, Dt)] = np.array([fup(sum(uD[dkey([Dm[q] - Dn[q] for q in range(f)])][n] for Dm in Dmem[Dt]))
                                         for n in range(N0 + 1)])
            g = float(np.dot(LA[(An, Ak)], UD[(Dn, Dt)])) + tail_unit * len(Amem[Ak]) * len(Dmem[Dt])
            row[j] = g * (1 + 1e-9)          # a-priori: 25-term nonnegative dot products and sums of <= 48 floats
        return row
    Grow = {z: Gbar_row(*z) for z in nexts}
    M = np.zeros((nS, nS))
    for k in range(nS):
        for key, q in Qbar[k].items():
            M[k] += q * Grow[key]
    M *= (1 + 1e-9)
    rhoM = max(abs(np.linalg.eigvals(M)))
    print(f"[engine-t] Mbar assembled: float rho(Mbar) = {rhoM:.6f} [{time.time()-T0:.0f}s]", flush=True)
    # ---------------------------------------------------------------- 5. far region
    far = far_region(d, f, c, fam, T, kapf, pa, pb, coef, RA, RD, RS, R3, NF, nproc, T0)
    etaF = far['etaF']
    # ---------------------------------------------------------------- 6. certificate
    Gnear = np.array([Grow[s] for s in states])
    cert = certificate(M, eta, Gnear, etaF, rho, kapf, T, T0)
    pickle.dump(dict(d=d, p=p, gK=gK, gh=gh, xbar=xbar, f=f, c=c, y=y, RA=RA, RD=RD, N0=N0, states=states, M=M,
                     eta=eta, Gnear=Gnear, far=far, cert=cert, rho=rho, rho_c=rho_c, kappa=kap, t=t,
                     Qbar=Qbar, xmax=xmax), open(out, "wb"))
    print(f"[engine-t] done [{time.time()-T0:.0f}s]", flush=True)


def crude_upsilon(f, c, fam, T, kapf, D, r):
    """upper bound, valid for EVERY lateral offset A with |A|_1 >= r and no shared points (r >= 2c+1 if D = 0),
    for E_z[e^{kappa B} - 1] <= kappa E[B] exp(kappa B_sup) (Lemma 5.10)."""
    D = np.array(D)
    D1 = int(np.abs(D).sum())
    if D1 == 0:
        assert r >= 2 * c + 1
    Fa, Fb = [], []
    for i in range(f):
        e = np.zeros(f, np.int64); e[i] = 1
        Fa.append((int(np.abs(D - e).sum()), min(int(np.abs(D - e - np.eye(f, dtype=np.int64)[a]).sum()) for a in range(f))))
        Fb.append((int(np.abs(D + e).sum()), min(int(np.abs(D + e + np.eye(f, dtype=np.int64)[a]).sum()) for a in range(f))))
    P2 = lambda rr, F: T.PHI2[min(max(rr, 0), T.Rmax), min(F, T.Fmax)]
    P3 = lambda rr: T.PSI3[min(max(rr, 0), T.Rmax), min(D1, T.Dmax)]
    lenlaw = defaultdict(Fr)
    for w, pw in zip(fam.W, fam.P):
        lenlaw[len(w)] += pw
    EB, Bsup = 0.0, 0.0
    for l1 in range(c + 1):
        for l2 in range(c + 1):
            pr = fup(lenlaw[l1] * lenlaw[l2] / (f * f))
            for ip in range(f):
                for i in range(f):
                    B = 0.0
                    for k in range(l1 + 1):                    # walk-1 vertices (lateral distance to x' >= r - k)
                        if D1 == 0:
                            B += sum(T.bt[max(1, r - k - kp)] for kp in range(l2 + 1))
                        else:
                            B += P2(r - k, D1)                 # same level, all walk-2 points
                        B += P2(r - k - l2, Fa[ip][0]) + P2(r - k - l2 - c, Fa[ip][1]) + P3(r - k)
                    for kp in range(l2 + 1):                   # walk-2 vertices
                        B += P2(r - kp - l1, Fb[i][0]) + P2(r - kp - l1 - c, Fb[i][1]) + P3(r - kp)
                    B *= (1 + 1e-13)
                    EB += pr * B
                    Bsup = max(Bsup, B)
    x = kapf * Bsup * (1 + 1e-15)
    assert x <= 1
    return kapf * EB * exp_up(x) * (1 + 1e-12), Bsup


def shell_upsilon(c, T, kapf, D1):
    """sup over all D with |D|_1 = D1 >= 4 and all A of E[e^{kappa B} - 1] (Lemma 5.10(c))."""
    P2 = lambda F: T.PHI2[0, min(F, T.Fmax)]
    per = P2(D1) + P2(D1 - 1) + P2(D1 - 2) + T.PSI3[0, min(D1, T.Dmax)]
    B = 2 * (c + 1) * per * (1 + 1e-13)
    x = kapf * B * (1 + 1e-15)
    assert x <= 1
    return x * exp_up(x) * (1 + 1e-12)


def far_region(d, f, c, fam, T, kapf, pa, pb, coef, RA, RD, RS, R3, NF, nproc, T0):
    R2 = RA + 6
    inner = d_types(f, 2 * RS)
    ups = {}
    jobs = []
    for Dt in inner:
        D1 = sum(map(abs, Dt))
        if D1 <= 2 * RD:
            band = [Ak for Ak in a_types(R2) if RA + 1 <= sum(Ak)]
            ups[Dt] = crude_upsilon(f, c, fam, T, kapf, Dt, R2 + 1)[0]
        else:
            band = a_types(R3)
            ups[Dt] = crude_upsilon(f, c, fam, T, kapf, Dt, R3 + 1)[0]
        jobs += [(Ak, Dt) for Ak in band]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(nproc) as pool:
        res = pool.map(_job, jobs)
    for (Ak, Dt), (acc, xm) in zip(jobs, res):
        assert xm <= 1
        v = float(acc.sum()) * (1 + INFL) * (1 + 1e-9)
        ups[Dt] = max(ups[Dt], v)
    # forward Green function between inner types, exact to NF steps plus tail
    TT = get_tt(f)
    TuF = TT.T_u(NF)
    u0 = TT.u0
    mem = {Dt: d_orbit(Dt) for Dt in inner}
    need = set()
    for Da in inner:
        for Db in inner:
            for Dm in mem[Db]:
                need.add(dkey([Dm[q] - Da[q] for q in range(f)]))
    Usum = {z: sum(u_exact_table(f, NF, z)) for z in need}
    hin = {}
    for Da in inner:
        tot = Fr(0)
        for Db in inner:
            g = sum(Usum[dkey([Dm[q] - Da[q] for q in range(f)])] for Dm in mem[Db]) + len(mem[Db]) * TuF
            tot += g * Fr(ups[Db])
        hin[Da] = fup(tot)
    h_in = max(hin.values())
    # outer shells R > RS
    Rcut = 60
    Sout = Fr(0)
    umax = 0.0
    for Rr in range(RS + 1, Rcut + 1):
        us = shell_upsilon(c, T, kapf, 2 * Rr)
        umax = max(umax, us)
        Sout += IG.N_radius(f, Rr) * Fr(us)
    # tail R > Rcut: per-vertex boost at |D|_1 = 2R is <= g(3 b(2R-2) + 2 b(R)/(1-q)) <= g(3 + 2/(1-q)) b(R) and
    # b(R) <= b(Rb) q^{R-Rb}; N_R <= C(R+f-1, f-1)^2; ratio test on N_R x e^x.
    q = T.q
    Rb = 20
    g = 1 + 2 * q / (1 - q)
    bRb = Fr(T.B[Rb], 1 << 300)
    C = Fr(kapf) * 2 * (c + 1) * g * (3 + 2 / (1 - q)) * bRb / q ** Rb
    xR = lambda Rr: C * q ** Rr
    R1 = Rcut + 1
    assert xR(R1) <= 1
    e1 = Fr(3)                                   # e^x <= 3 for x <= 1
    first = Fr(math.comb(R1 + f - 1, f - 1) ** 2) * xR(R1) * e1
    rat = Fr(R1 + f, R1 + 1) ** (2 * (f - 1)) * q
    assert rat < 1
    Sout += first / (1 - rat)
    eta_out = Fr(0)
    for n, un in enumerate(u0):
        eta_out += min(Fr(umax), un * Sout)
    assert len(u0) == TT.NEX + 1
    eta_out += Sout * TT.tail
    etaF = fup(Fr(h_in) + eta_out)
    print(f"[engine-t] far: etaF <= {etaF:.4e} (inner {h_in:.4e} max over {len(inner)} forward types, outer "
          f"{float(eta_out):.3e}; Sout={float(Sout):.3e}) [{time.time()-T0:.0f}s]", flush=True)
    return dict(etaF=etaF, ups=ups, h_in=h_in, eta_out=float(eta_out), hin=hin)


def certificate(M, eta, Gnear, etaF, rho, kapf, T, T0):
    nS = len(eta)
    rhoM = max(abs(np.linalg.eigvals(M)))
    delta = 1e-9

    def trial(lam0):
        try:
            w = np.linalg.solve(lam0 * np.eye(nS) - M, eta + delta)
        except np.linalg.LinAlgError:
            return None
        if not (w > 0).all():
            return None
        lam = lam0 * (1 + 1e-7)
        Gam = float(np.max(Gnear @ w)) * (1 + 1e-7)
        if lam <= etaF:
            return None
        H = Gam / (lam - etaF) * (1 + 1e-7)
        pre = (np.all(M @ w + H * etaF * eta <= lam * w * (1 - 1e-12)) and float(np.max(Gnear @ w)) + H * etaF
               <= lam * H * (1 - 1e-12) and lam < 1)
        ok = pre and exact_check(M, eta, Gnear, etaF, w, lam, H)
        return dict(lam0=lam0, lam=lam, w=w, H=H, Gam=Gam, ok=ok)

    hi = trial(0.999)
    if hi is None or not hi['ok']:
        print(f"[engine-t] NOT CERTIFIED: rho(Mbar)={rhoM:.6f}, etaF={etaF:.4e}", flush=True)
        return dict(ok=False, rhoM=rhoM)
    lo_, hi_ = rhoM * (1 + 1e-6), 0.999
    best = hi
    for _ in range(30):
        mid = (lo_ + hi_) / 2
        r = trial(mid)
        if r is not None and r['ok']:
            hi_, best = mid, r
        else:
            lo_ = mid
    rho_f = float(rho)
    Cfin = (1 / rho_f) * exp_up(kapf * (T.bt[1] + 2 * T.PSI3[0, 0]) * (1 + 1e-15)) * (1 + 1e-12)

    def theta(V):
        Lam = max(float(np.max(eta / V['w'])), 1 / V['H'])
        bound = Cfin * (1 + Lam * V['lam'] * V['H'] / (1 - V['lam']))
        return 0.99 / bound * (1 - 1e-9), Lam, bound
    th = theta(best)
    bestth = (th, best)
    for fr in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        l0 = best['lam0'] + fr * (0.999 - best['lam0'])
        r = trial(l0)
        if r is not None and r['ok']:
            t2 = theta(r)
            if t2[0] > bestth[0][0]:
                bestth = (t2, r)
    print(f"[engine-t] CERTIFIED (exact rational check of (i),(ii)): lam = {best['lam']:.6f}; rho(Mbar)={rhoM:.6f}, "
          f"etaF<={etaF:.4e}, Gamma_w<={best['Gam']:.4f}, H={best['H']:.4f}; theta_* >= {th[0]:.6f}", flush=True)
    print(f"[engine-t]   best theta over the lam0 scan: lam = {bestth[1]['lam']:.6f}, sup_n E W_n^2 <= {bestth[0][2]:.4f}, "
          f"theta_* >= {bestth[0][0]:.6f}; C_fin <= {Cfin:.4f}", flush=True)
    return dict(ok=True, lam=best['lam'], w=best['w'], H=best['H'], Gam=best['Gam'], theta=th[0], rhoM=rhoM,
                lam_th=bestth[1]['lam'], theta_best=bestth[0][0], bound_best=bestth[0][2], w_th=bestth[1]['w'],
                H_th=bestth[1]['H'], Cfin=Cfin)


def exact_check(M, eta, Gnear, etaF, w, lam, H):
    """exact rational verification of (i) Mbar w + H etaF etabar <= lam w (every row) and
    (ii) max_z (Gbar_C' w)(z) + H etaF <= lam H, lam < 1, all floats read as exact dyadic rationals."""
    lamF, HF, eF = Fr(lam), Fr(H), Fr(etaF)
    if not lamF < 1:
        return False
    wF = [Fr(float(x)) for x in w]
    nS = len(wF)
    for s in range(nS):
        lhs = sum(Fr(float(M[s, j])) * wF[j] for j in range(nS)) + HF * eF * Fr(float(eta[s]))
        if not lhs <= lamF * wF[s]:
            return False
    Gam = max(sum(Fr(float(Gnear[s, j])) * wF[j] for j in range(nS)) for s in range(nS))
    return Gam + HF * eF <= lamF * HF


if __name__ == "__main__":
    main()

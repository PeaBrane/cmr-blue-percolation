"""Rigorous certificate, part 2: far region constants and the Collatz-Wielandt verification.

Conditions of Theorem 5.9 of the first manuscript:
  (i)  for every near state z in C':  (Mbar w)(z) + H * etaF * eta_bar(z) <= lam * w(z)
  (ii) Gamma_w + H * etaF <= lam * H,   Gamma_w >= max_{z in C'} sum_{z' in C'} G(z,z') w(z')
with etaF >= sup_{z''} sum_{far z'} G(z'',z') Ubar(z').  Then E W_n^2 <= e^{kappa eps} C_fin (1 + Lam lam H/(1-lam)),
Lam = max(max_z eta_bar(z)/w(z), 1/H).
"""
import sys, pickle, math, itertools, time
import numpy as np
from fractions import Fraction as Fr
from common import fup, add_up, mul_up, sum_up, div_up, local_constants, POINTS, INF, IG
from family import Family, Tables
from green_rig import u_table, u0_tail, robbins_tail, dkey, D_members

RS = 4          # forward shells |D|_1 <= 2*RS are treated exactly in etaF
NF = 60         # forward Green function exact up to NF steps


def fwd_types(f, R):
    """S_f-types (sorted desc tuples) of root-lattice vectors with |D|_1 = 2R."""
    out = set()
    def parts(n, maxp):
        if n == 0:
            yield ()
            return
        for k in range(min(n, maxp), 0, -1):
            for rest in parts(n - k, k):
                yield (k,) + rest
    for pos in parts(R, R):
        for neg in parts(R, R):
            if len(pos) + len(neg) <= f:
                out.add(tuple(sorted(list(pos) + [0] * (f - len(pos) - len(neg)) + [-x for x in neg], reverse=True)))
    return sorted(out)


def F_values(D, f):
    eyes = np.eye(f, dtype=np.int64)
    D = np.array(D)
    res = []
    for i in range(f):
        for ip in range(f):
            F1a = int(np.abs(D - eyes[ip]).sum()); F1b = int(np.abs(D + eyes[i]).sum())
            F2a = min(int(np.abs(D - eyes[ip] - eyes[q]).sum()) for q in range(f))
            F2b = min(int(np.abs(D + eyes[i] + eyes[q]).sum()) for q in range(f))
            res.append((F1a, F1b, F2a, F2b))
    return res


def upsilon(tab, fam, kap_up, D, r):
    """upper bound for Ubar(A, D) = E[e^{kappa B} - 1] valid for all A with |A|_1 >= r (no shared points):
    kappa * E[B_lin] * exp(kappa * B_sup), with every distance bounded below monotonically in r."""
    c, f = fam.c, fam.f
    bt, Phi2, psi3 = tab.bt, tab.Phi2, tab.psi3
    D1 = int(np.abs(np.array(D)).sum())
    rP = Phi2.shape[0] - 1; r3 = psi3.shape[0] - 1
    Fv = F_values(D, f)
    def same_level(i1, l2):
        # pairs of walk-1 vertex i1 with the l2+1 points of walk 2: sum_j b(D1 + |rho - j|), rho >= r - i1
        if D1 == 0:
            return sum_up(bt[max(1, r - i1 - j)] for j in range(l2 + 1))   # far A: r - i1 - j >= r - 2c >= 1
        return Phi2[min(max(0, r - i1), rP), D1]
    def vert_terms(i1, l2, F1, F2):
        return sum_up([Phi2[min(max(0, r - i1 - l2), rP), F1], Phi2[min(max(0, r - i1 - l2 - c), rP), F2],
                       psi3[min(max(0, r - i1), r3), D1]])
    Blin = 0.0
    Bsup = 0.0
    for l1 in range(c + 1):
        for l2 in range(c + 1):
            pr = fup(fam.len_law[l1] * fam.len_law[l2])
            for (F1a, F1b, F2a, F2b) in Fv:
                B = 0.0
                for i1 in range(l1 + 1):
                    B = add_up(B, same_level(i1, l2))
                    B = add_up(B, vert_terms(i1, l2, F1a, F2a))
                for j2 in range(l2 + 1):
                    B = add_up(B, vert_terms(j2, l1, F1b, F2b))
                Blin = add_up(Blin, mul_up(pr, div_up(B, float(f * f))))
                Bsup = max(Bsup, B)
    x = mul_up(kap_up, Bsup)
    assert x <= 1
    ex = fup(IG.exp_hi(Fr(x)))
    return mul_up(mul_up(kap_up, Blin), ex), Bsup


def upsilon_shell(tab, fam, kap_up, D1):
    """sup over all D with |D|_1 = D1 >= 4 and all A: kappa*B_sup*exp(kappa*B_sup), F1 >= D1-1, F2 >= D1-2."""
    c = fam.c
    Phi2, psi3 = tab.Phi2, tab.psi3
    per_vertex = sum_up([Phi2[0, D1], Phi2[0, D1 - 1], Phi2[0, D1 - 2], psi3[0, min(D1, psi3.shape[1] - 1)]])
    B = mul_up(float(2 * (c + 1)), per_vertex)
    x = mul_up(kap_up, B)
    assert x <= 1
    return mul_up(x, fup(IG.exp_hi(Fr(x))))


_FAR_ARGS = None


def _far_job(job):
    from enum_rig import enum_state
    Ak, Dt = job
    PT, LN, EN, mass_up, f, ones, kap_up, bt, Phi2, psi3, c, PAIRIDX0, coef, rem = _FAR_ARGS
    acc, xm = enum_state(PT, LN, EN, mass_up, np.array(Ak, np.int64), np.array(Dt, np.int64), f, ones, ones, kap_up,
                         bt, Phi2, psi3, c, PAIRIDX0, 1, coef, rem)
    assert xm <= 1
    return sum_up(acc.ravel())


def main(pkl):
    t0 = time.time()
    R = pickle.load(open(pkl, "rb"))
    prm = R["params"]
    d, f, c, y = prm["d"], prm["f"], prm["c"], prm["y"]
    K = R["K"]
    fam = Family(d, f, c, y)
    tab = Tables(d, K["t"], K["qs"], c)
    kap_up = fup(K["kappa"])
    states, M, eta_bar, Gn = R["states"], R["M"], R["eta_bar"], R["Gn"]
    nS = len(states)
    # ---------------- far constants
    # (F1): D-types of C' (|D|_1 <= 2 RD), far lateral |A|_1 >= RA+1: exact kernel averages for RA+1 <= |A|_1 <= R2
    #       (upward-rounded, boost only: no shared points there) and the crude averaged bound at r = R2+1 beyond.
    RA = prm["RA"]; RD = prm.get("RD", 1); R2 = RA + 6
    from enum_rig import enum_state
    import itertools as _it
    from certify import akey
    ones = np.ones(len(R["pa"]))
    KE = 16
    coef = np.array([fup(Fr(1, math.factorial(k))) for k in range(KE + 1)])
    rem = fup(Fr(3, math.factorial(KE + 1)))
    PAIRIDX0 = np.zeros((fam.nW, fam.nW), np.int64)
    ups = {}
    far_orbits = sorted({akey(A) for A in _it.product(range(-R2, R2 + 1), repeat=fam.m)
                         if RA + 1 <= sum(map(abs, A)) <= R2})
    import multiprocessing as mp, os as _os
    cache = pkl + ".ups.pkl"
    if _os.path.exists(cache):
        ups_far = pickle.load(open(cache, "rb"))
    else:
        ups_far = None
    jobs = [(Ak, Dt) for Dt in R["Dtypes"] for Ak in far_orbits] if ups_far is None else []
    global _FAR_ARGS
    _FAR_ARGS = (fam.PT, fam.LN, fam.EN, fam.mass_up, f, ones, kap_up, tab.bt, tab.Phi2, tab.psi3, c, PAIRIDX0, coef, rem)
    vals = []
    if jobs:
        with mp.get_context("fork").Pool(int(__import__("os").environ.get("NPROC", "6"))) as pool:
            vals = pool.map(_far_job, jobs)
    if ups_far is None:
        for Dt in R["Dtypes"]:
            ups[Dt] = upsilon(tab, fam, kap_up, Dt, R2 + 1)[0]
        for (Ak, Dt), v in zip(jobs, vals):
            ups[Dt] = max(ups[Dt], v)
        pickle.dump({Dt: ups[Dt] for Dt in R["Dtypes"]}, open(cache, "wb"))
    else:
        ups.update(ups_far)
    # (F2): all types with 2 RD < |D|_1 <= 2 RS, any A
    for Rr in range(RD + 1, RS + 1):
        for Dt in fwd_types(f, Rr):
            ups[Dt], _ = upsilon(tab, fam, kap_up, Dt, 0)
    types_in = [Dt for Rr in range(0, RS + 1) for Dt in fwd_types(f, Rr)]
    print("upsilon: " + ", ".join(f"{Dt}:{ups[Dt]:.2e}" for Dt in types_in[:10]) + " ...", flush=True)
    # forward Green function between S-types, exact to NF steps, tail T_u(NF)
    zs = set()
    mem = {Dt: D_members(Dt) for Dt in types_in}
    for Da in types_in:
        for Db in types_in:
            for Dm in mem[Db]:
                zs.add(dkey(np.array(Dm) - np.array(Da)))
    Ut = u_table(f, NF, sorted(zs))
    TuF, u0list = u0_tail(f, NF)
    Usum = {z: sum(v) for z, v in Ut.items()}
    h = {}
    for Da in types_in:
        tot = 0.0
        for Db in types_in:
            g = sum(Usum[dkey(np.array(Dm) - np.array(Da))] for Dm in mem[Db]) + len(mem[Db]) * TuF
            tot = add_up(tot, mul_up(fup(g), ups[Db]))
        h[Da] = tot
    h_in = max(h.values())
    # outside shells: sum_{R > RS} N_R * ups_shell(2R), geometric tail beyond Rcut
    Rcut = 60
    Sout = 0.0; umax_out = 0.0
    for Rr in range(RS + 1, Rcut + 1):
        us = upsilon_shell(tab, fam, kap_up, 2 * Rr)
        umax_out = max(umax_out, us)
        Sout = add_up(Sout, mul_up(float(IG.N_radius(f, Rr)), us))
    # tail R > Rcut (rigorous): per-vertex boost at |D|_1 = 2R is <= Cq * q^R with
    #   Phi2[0,F] <= b(F)(1+2q/(1-q)),  psi3[0,2R] <= (1+2q/(1-q)) 2 b(R)/(1-q),  b(F) <= b(Rb) q^{F-Rb};
    #   so ups_shell(2R) <= x e^x with x = kappa*2(c+1)*Cq*q^R <= 1;  N_R <= C(R+f-1, f-1)^2.
    q = K["qs"]
    Rb = 10
    bRb = __import__("common").b_exact(d, K["t"], Rb)                 # exact; b(F) <= b(Rb) q^{F-Rb}, F >= Rb
    g = 1 + 2 * q / (1 - q)
    Cq = (3 * g * bRb / q ** Rb + g * 2 * bRb / q ** Rb / (1 - q))   # per-vertex bound / q^R (uses 2R-2 >= R)
    xR = lambda Rr: Fr(fup(kap_up)) * 2 * (c + 1) * Cq * q ** Rr
    assert xR(Rcut + 1) <= 1
    first = Fr(math.comb(Rcut + 1 + f - 1, f - 1) ** 2) * xR(Rcut + 1) * IG.exp_hi(min(Fr(1), xR(Rcut + 1)))
    rat = Fr(Rcut + 1 + f, Rcut + 2) ** (2 * (f - 1)) * q
    assert rat < 1
    Sout = add_up(Sout, fup(first / (1 - rat)))
    eta_out = 0.0
    for n, un in enumerate(u0list):
        eta_out = add_up(eta_out, min(umax_out, mul_up(fup(un), Sout)))
    eta_out = add_up(eta_out, mul_up(Sout, fup(robbins_tail(f, len(u0list) // f))))
    etaF = add_up(h_in, eta_out)
    print(f"etaF <= {etaF:.3e}  (inner {h_in:.3e}, outer {eta_out:.3e}; Sout={Sout:.2e}) [{time.time()-t0:.0f}s]", flush=True)
    # ---------------- test vector and verification
    # w = (lam0 I - Mbar)^{-1} (eta_bar + delta): then (Mbar w) = lam0 w - eta_bar - delta and condition (i) reduces
    # to H etaF <= 1 (checked below with directed rounding anyway).  Scan lam0 upward from rho(Mbar).
    ev = np.linalg.eigvals(M)
    rhoM = max(abs(ev))
    delta = 1e-9
    def verify(lam0):
        w = np.linalg.solve(lam0 * np.eye(nS) - M, eta_bar + delta)
        if not (w > 0).all():
            return None
        Mw = np.zeros(nS)
        for s in range(nS):
            acc = 0.0
            for t in range(nS):
                acc = add_up(acc, mul_up(M[s, t], w[t]))
            Mw[s] = acc
        Gw = []
        for (Ak, Dt) in states:
            row = Gn[(Ak, Dt)]
            acc = 0.0
            for t in range(nS):
                acc = add_up(acc, mul_up(row[t], w[t]))
            Gw.append(acc)
        Gamma_w = max(Gw)
        lam = lam0 * (1 + 1e-9)
        if lam <= etaF:
            return None
        H = math.nextafter(div_up(Gamma_w, math.nextafter(lam - etaF, -INF)), INF)
        ok_i = all(add_up(Mw[s], mul_up(mul_up(H, etaF), eta_bar[s])) <= math.nextafter(lam * w[s], -INF)
                   for s in range(nS))
        ok_ii = add_up(Gamma_w, mul_up(H, etaF)) <= math.nextafter(lam * H, -INF)
        return dict(lam=lam, w=w, Mw=Mw, Gamma_w=Gamma_w, H=H, ok=ok_i and ok_ii and lam < 1, ok_i=ok_i, ok_ii=ok_ii)
    best = None
    lo_, hi_ = rhoM * (1 + 1e-6), 0.999999
    res_hi = verify(hi_)
    if res_hi is None or not res_hi["ok"]:
        print(f"d={d} t={prm['tkey']}: NOT certified (rho(Mbar)={rhoM:.6f}, etaF={etaF:.3e})", flush=True)
        return dict(ok=False, rhoM=rhoM, etaF=etaF)
    for _ in range(40):
        mid = (lo_ + hi_) / 2
        r_ = verify(mid)
        if r_ is not None and r_["ok"]:
            hi_ = mid; res_hi = r_
        else:
            lo_ = mid
    rho = K["rho"]
    beta_fin = add_up(tab.bt[1], mul_up(2.0, tab.psi3[0, 0]))
    Cfin = mul_up(fup(1 / rho), fup(IG.exp_hi(Fr(mul_up(kap_up, beta_fin)))))
    def theta_of(V):
        lam, w, H = V["lam"], V["w"], V["H"]
        Lam = max(max(div_up(eta_bar[s], w[s]) for s in range(nS)), div_up(1.0, H))
        bound = mul_up(Cfin, add_up(1.0, div_up(mul_up(mul_up(Lam, lam), H), math.nextafter(1 - lam, -INF))))
        return math.nextafter(0.99 / bound, -INF), Lam, bound
    lam_min = res_hi["lam"]
    best_V = res_hi; best_th = theta_of(res_hi)
    for frac in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        l0 = lam_min + frac * (0.999 - lam_min)
        V_ = verify(l0)
        if V_ is not None and V_["ok"]:
            th_ = theta_of(V_)
            if th_[0] > best_th[0]:
                best_V, best_th = V_, th_
    V = res_hi
    lam, w, H, Gamma_w = V["lam"], V["w"], V["H"], V["Gamma_w"]
    theta, Lam, bound = theta_of(V)
    print(f"d={d} t={prm['tkey']} f={f} c={c} y={y} RA={RA} N0={prm['N0']}: rho(Mbar) = {rhoM:.6f}; etaF <= {etaF:.3e}", flush=True)
    print(f"   CERTIFIED lam = {lam:.6f} < 1  [(i) {V['ok_i']}, (ii) {V['ok_ii']}];  Gamma_w <= {Gamma_w:.4f}, H = {H:.4f}, H*etaF = {H*etaF:.4f}", flush=True)
    print(f"   Lam = {Lam:.4f}; C_fin <= {Cfin:.4f}; sup_n E W_n^2 <= {bound:.4f}; theta_* >= {theta:.6f}", flush=True)
    print(f"   best theta over the lam0 scan: lam = {best_V['lam']:.6f} (also certified: (i) {best_V['ok_i']}, (ii) {best_V['ok_ii']}),"
          f" sup_n E W_n^2 <= {best_th[2]:.4f}, theta_* >= {best_th[0]:.6f}", flush=True)
    # exact re-verification of the final linear inequalities (floats are exact dyadic rationals)
    def exact_check(V):
        lamF, HF, eF = Fr(V["lam"]), Fr(V["H"]), Fr(etaF)
        wF = [Fr(float(x)) for x in V["w"]]
        ok1 = True
        for s in range(nS):
            lhs = sum(Fr(float(M[s, t])) * wF[t] for t in range(nS)) + HF * eF * Fr(float(eta_bar[s]))
            if not lhs <= lamF * wF[s]:
                ok1 = False; break
        Gam = max(sum(Fr(float(Gn[(Ak, Dt)][t])) * wF[t] for t in range(nS)) for (Ak, Dt) in states)
        ok2 = Gam + HF * eF <= lamF * HF
        return ok1 and ok2 and lamF < 1
    ex1 = exact_check(V); ex2 = exact_check(best_V)
    print(f"   exact-rational re-check of (i),(ii): minimal-lam vector {ex1}; best-theta vector {ex2}", flush=True)
    out = dict(lam=lam, theta=theta, bound=bound, lam_best=best_V["lam"], theta_best=best_th[0], bound_best=best_th[2],
               etaF=etaF, rhoM=rhoM, Gamma_w=Gamma_w, H=H, Cfin=Cfin, exact_ok=(ex1, ex2), ups=ups,
               w=V["w"], w_best=best_V["w"], H_best=best_V["H"])
    pickle.dump(out, open(pkl + ".cert.pkl", "wb"))
    return out


if __name__ == "__main__":
    main(sys.argv[1])

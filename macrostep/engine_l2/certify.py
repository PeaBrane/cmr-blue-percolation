"""Rigorous certificate for the macrostep second-moment criterion (the directed-rounding program), part 1: the near matrix.

States of the pair chain: z = (A, D), A in Z^3 lateral offset, D in the forward root lattice of Z^f.
Near set C' = {D in {0, roots}, |A|_1 <= RA}, orbit-reduced under B_3 (lateral signed permutations) x S_f.
Outputs (pickled): Qbar, eta_bar, Gbar tables, Mbar, Perron test vector w, (Mbar w), Gamma_w.
"""
import sys, time, pickle, itertools, math
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from common import (fup, add_up, mul_up, sum_up, local_constants, POINTS, INF)
from family import Family, Tables
from enum_rig import enum_state
from green_rig import u_table, u0_tail, dkey, D_members, lateral_table


def akey(A):
    return tuple(sorted((abs(int(x)) for x in A), reverse=True))


def A_members(Ak):
    m = len(Ak)
    return sorted({tuple(sg * v for sg, v in zip(sgns, perm)) for perm in itertools.permutations(Ak)
                   for sgns in itertools.product((1, -1), repeat=m)})


def dot_up(u, v):
    s = 0.0
    for a, b in zip(u, v):
        s = math.nextafter(s + math.nextafter(a * b, INF), INF)
    return s


def build(d, tkey, f, c, y, RA=8, N0=24, RD=1, verbose=True):
    t0 = time.time()
    P = POINTS[(d, tkey)]
    K = local_constants(d, P["p"], P["gK"], P["gh"], P["xbar"])
    rho, rho_c, p, kap = K["rho"], K["rho_c"], K["p"], K["kappa"]
    import os as _os0
    if _os0.environ.get("P_OVERRIDE"):          # robustness variant: a smaller (hence still valid) bond floor p
        p = Fr(_os0.environ["P_OVERRIDE"]); assert p <= K["p"]
        K = dict(K); K["p"] = p
    if _os0.environ.get("NO_RHOC") == "1":       # robustness variant: no shared-edge refinement (rho_e = rho_-)
        rho_c = rho
        K = dict(K); K["rho_c"] = rho
    fam = Family(d, f, c, y)
    m = fam.m
    tab = Tables(d, K["t"], K["qs"], c)
    a = 1 / rho
    bb = rho / (rho_c * p)
    pa = np.array([fup(a ** n) for n in range(2 * c + 4)])
    pb = np.array([fup(bb ** n) for n in range(2 * c + 4)])
    kap_up = fup(kap)
    KE = 16
    coef = np.array([fup(Fr(1, math.factorial(k))) for k in range(KE + 1)])
    rem = fup(Fr(3, math.factorial(KE + 1)))
    # states
    from certify2 import fwd_types
    Dtypes = [Dt for Rr in range(0, RD + 1) for Dt in fwd_types(f, Rr)]
    Akeys = sorted({akey(A) for A in itertools.product(range(-RA, RA + 1), repeat=m) if sum(map(abs, A)) <= RA},
                   key=lambda k: (sum(k), k))
    states = [(Ak, Dt) for Dt in Dtypes for Ak in Akeys]
    sidx = {s: i for i, s in enumerate(states)}
    # dA index
    dAs = sorted({tuple(int(x) for x in (fam.EN[w1] - fam.EN[w2])) for w1 in range(fam.nW) for w2 in range(fam.nW)})
    dAidx = {v: k for k, v in enumerate(dAs)}
    PAIRIDX = np.array([[dAidx[tuple(int(x) for x in (fam.EN[w1] - fam.EN[w2]))] for w2 in range(fam.nW)]
                        for w1 in range(fam.nW)], np.int64)
    eyes = np.eye(f, dtype=np.int64)
    Qbar = defaultdict(float)
    eta_bar = np.zeros(len(states))
    xmax_all = 0.0
    for si, (Ak, Dt) in enumerate(states):
        A = np.array(Ak, np.int64); D = np.array(Dt, np.int64)
        acc, xmax = enum_state(fam.PT, fam.LN, fam.EN, fam.mass_up, A, D, f, pa, pb, kap_up, tab.bt, tab.Phi2,
                               tab.psi3, c, PAIRIDX, len(dAs), coef, rem)
        xmax_all = max(xmax_all, xmax)
        loc = defaultdict(float)
        for kk, dA in enumerate(dAs):
            Ap = akey([A[q] + dA[q] for q in range(m)])
            for i in range(f):
                for ip in range(f):
                    v = acc[kk, i, ip]
                    if v == 0.0:
                        continue
                    Dn = dkey(D + eyes[i] - eyes[ip])
                    loc[(Ap, Dn)] = add_up(loc[(Ap, Dn)], v)
        for key, v in loc.items():
            Qbar[(si,) + key] = v
            eta_bar[si] = add_up(eta_bar[si], v)
        if verbose and si % 10 == 0:
            print(f"  state {si}/{len(states)} {Ak} {Dt}: eta={eta_bar[si]:.5f} xmax={xmax:.4f} [{time.time()-t0:.0f}s]", flush=True)
    assert xmax_all <= 1.0
    # Green function tables
    nexts = sorted({(k[1], k[2]) for k in Qbar})
    starts = [(Ak, Dt) for (Ak, Dt) in states]
    allAstart = sorted({x[0] for x in nexts} | {s[0] for s in starts})
    allDstart = sorted({x[1] for x in nexts} | {s[1] for s in starts})
    Amem = {Ak: A_members(Ak) for Ak in Akeys}
    Dmem = {Dt: D_members(Dt) for Dt in Dtypes}
    zA = set()
    for An in allAstart:
        for Bk in Akeys:
            for B in Amem[Bk]:
                zA.add(akey([B[q] - An[q] for q in range(m)]))
    zA = sorted(zA, key=lambda k: (sum(k), k))
    zD = set([tuple([0] * f)])
    for Dn in allDstart:
        for Dt in Dtypes:
            for Dm in Dmem[Dt]:
                zD.add(dkey(np.array(Dm) - np.array(Dn)))
    if verbose:
        print(f"  #states={len(states)} #next-types={len(nexts)} #zA={len(zA)} #zD={len(zD)} [{time.time()-t0:.0f}s]", flush=True)
    import os as _os
    cache = _os.environ.get("LATCACHE")
    Lt = None
    if cache and _os.path.exists(cache):
        LC = pickle.load(open(cache, "rb"))
        if LC["c"] == c and LC["y"] == Fr(y) and LC["N0"] >= N0 and all(z in LC["L"] for z in zA):
            Lt = {z: LC["L"][z][:N0 + 1] for z in zA}
            if verbose:
                print(f"  lateral table loaded from cache {cache}", flush=True)
    if Lt is None:
        Lt = lateral_table(fam, N0, zA)
    Ut = u_table(f, N0, sorted(zD))
    Tu, _ = u0_tail(f, N0)
    L0N = Lt[tuple([0] * m)][N0]
    tail_unit = mul_up(L0N, fup(Tu))
    if verbose:
        print(f"  lateral+forward tables done: l_N0(0)<={L0N:.3e} T_u<={float(Tu):.3e} tail_unit={tail_unit:.3e} [{time.time()-t0:.0f}s]", flush=True)
    LA = {}
    for An in allAstart:
        for Bk in Akeys:
            vec = np.zeros(N0 + 1)
            for B in Amem[Bk]:
                lv = Lt[akey([B[q] - An[q] for q in range(m)])]
                for n in range(N0 + 1):
                    vec[n] = add_up(vec[n], lv[n])
            LA[(An, Bk)] = vec
    UD = {}
    for Dn in allDstart:
        for Dt in Dtypes:
            ex = [sum(Ut[dkey(np.array(Dm) - np.array(Dn))][n] for Dm in Dmem[Dt]) for n in range(N0 + 1)]
            UD[(Dn, Dt)] = [fup(x) for x in ex]
    def Gbar(An, Dn, Bk, Dt):
        g = dot_up(UD[(Dn, Dt)], LA[(An, Bk)])
        return add_up(g, mul_up(tail_unit, float(len(Amem[Bk]) * len(Dmem[Dt]))))
    nS = len(states)
    Gn = {}
    for (An, Dn) in set(nexts) | set(starts):
        Gn[(An, Dn)] = np.array([Gbar(An, Dn, Bk, Dt) for (Bk, Dt) in states])
    M = np.zeros((nS, nS))
    for (si, An, Dn), q in Qbar.items():
        row = Gn[(An, Dn)]
        for sj in range(nS):
            M[si, sj] = add_up(M[si, sj], mul_up(q, row[sj]))
    if verbose:
        print(f"  Mbar assembled [{time.time()-t0:.0f}s]", flush=True)
    return dict(K=K, fam=fam, tab=tab, states=states, Qbar=dict(Qbar), eta_bar=eta_bar, M=M, Gn=Gn,
                Akeys=Akeys, Dtypes=Dtypes, Amem={k: len(v) for k, v in Amem.items()},
                Dmem={k: len(v) for k, v in Dmem.items()}, N0=N0, tail_unit=tail_unit, L0N=L0N, Tu=Tu,
                params=dict(d=d, tkey=tkey, f=f, c=c, y=Fr(y), RA=RA, N0=N0, RD=RD), xmax=xmax_all,
                pa=pa, pb=pb, kap_up=kap_up)


if __name__ == "__main__":
    d, tkey, f, c = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    y = Fr(sys.argv[5]); RA = int(sys.argv[6]); N0 = int(sys.argv[7]); out = sys.argv[8]
    RD = int(sys.argv[9]) if len(sys.argv) > 9 else 1
    res = build(d, tkey, f, c, y, RA=RA, N0=N0, RD=RD)
    M = res["M"]
    ev, evec = np.linalg.eig(M)
    i0 = np.argmax(abs(ev))
    print(f"float Perron root of Mbar: {abs(ev[i0]):.6f}", flush=True)
    res.pop("fam"); res.pop("tab")
    with open(out, "wb") as fh:
        pickle.dump(res, fh)

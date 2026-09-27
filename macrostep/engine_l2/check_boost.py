"""Independent checks (not part of the certificate's arithmetic):
 (1) pure-Python re-implementation of the per-macrostep weight from the definitions, summed over all choices at a
     state, compared with the numba kernel's eta_bar(state);
 (2) Monte Carlo over random pairs of n-macrostep paths on Z^d: the actual shared counts and the actual pair sum
     sum_{x in V'\\V, y in V\\V'} b(|x-y|_1) are compared with the macrostep decomposition
     (sum_J nV_J + 1{S_n=0}, sum_J nE_J + fwd_J, sum_J B_J + beta_fin); the boost inequality must never fail.
"""
import sys, math, random, itertools
import numpy as np
from fractions import Fraction as Fr
from common import fup, local_constants, POINTS
from family import Family, Tables


def B_ref(fam, tab, A, D, w1, w2, i, ip):
    """per-macrostep (nV, nE, fwd, B) from the definitions (float, round-to-nearest)."""
    c, f, m = fam.c, fam.f, fam.m
    bt, Phi2, psi3 = tab.bt, tab.Phi2, tab.psi3
    P1 = [tuple(fam.PT[w1, k]) for k in range(fam.LN[w1] + 1)]
    P2 = [tuple(fam.PT[w2, k]) for k in range(fam.LN[w2] + 1)]
    e1, e2 = tuple(fam.EN[w1]), tuple(fam.EN[w2])
    D1 = sum(abs(x) for x in D)
    # walk-2 lateral points relative to walk-1 start: P2 - A
    Q2 = [tuple(p[q] - A[q] for q in range(m)) for p in P2]
    shared1 = [D1 == 0 and (p in Q2) for p in P1]
    shared2 = [D1 == 0 and (p in P1) for p in Q2]
    nV = sum(shared1)
    E1 = {frozenset((P1[k], P1[k + 1])) for k in range(len(P1) - 1)}
    E2 = {frozenset((Q2[k], Q2[k + 1])) for k in range(len(Q2) - 1)}
    nE = len(E1 & E2) if D1 == 0 else 0
    Ap = tuple(A[q] + e1[q] - e2[q] for q in range(m))
    fwd = 1 if (D1 == 0 and i == ip and not any(Ap)) else 0
    ei = [1 if q == i else 0 for q in range(f)]; eip = [1 if q == ip else 0 for q in range(f)]
    F1a = sum(abs(D[q] - eip[q]) for q in range(f)); F1b = sum(abs(D[q] + ei[q]) for q in range(f))
    F2a = min(sum(abs(D[q] - eip[q] - (1 if q == a else 0)) for q in range(f)) for a in range(f))
    F2b = min(sum(abs(D[q] + ei[q] + (1 if q == a else 0)) for q in range(f)) for a in range(f))
    B = 0.0
    for k1, p in enumerate(P1):
        if shared1[k1]: continue
        for k2, qq in enumerate(Q2):
            if shared2[k2]: continue
            B += bt[D1 + sum(abs(p[r] - qq[r]) for r in range(m))]
    # walk-1 vertices: r1 = dist to walk-2 next start (L' + e2 = L - A + e2), r0 = dist to walk-2 start (L - A)
    n2start = tuple(-A[q] + e2[q] for q in range(m)); c2start = tuple(-A[q] for q in range(m))
    for k1, p in enumerate(P1):
        if shared1[k1]: continue
        r1 = sum(abs(p[q] - n2start[q]) for q in range(m)); r0 = sum(abs(p[q] - c2start[q]) for q in range(m))
        B += sum(bt[F1a + abs(r1 - k)] for k in range(c + 1))
        B += Phi2[min(max(0, r1 - c), Phi2.shape[0] - 1), F2a]
        B += psi3[min(r0, psi3.shape[0] - 1), D1]
    n1start = e1; c1start = tuple([0] * m)
    for k2, qq in enumerate(Q2):
        if shared2[k2]: continue
        r1 = sum(abs(qq[q] - n1start[q]) for q in range(m)); r0 = sum(abs(qq[q] - c1start[q]) for q in range(m))
        B += sum(bt[F1b + abs(r1 - k)] for k in range(c + 1))
        B += Phi2[min(max(0, r1 - c), Phi2.shape[0] - 1), F2b]
        B += psi3[min(r0, psi3.shape[0] - 1), D1]
    return nV, nE, fwd, B


def eta_ref(fam, tab, K, A, D):
    f, c = fam.f, fam.c
    a = float(1 / K["rho"]); bb = float(K["rho"] / (K["rho_c"] * K["p"])); kap = float(K["kappa"])
    tot = 0.0
    for w1 in range(fam.nW):
        for w2 in range(fam.nW):
            pr = float(fam.Pl[fam.LN[w1]] * fam.Pl[fam.LN[w2]]) / f ** 2
            for i in range(f):
                for ip in range(f):
                    nV, nE, fwd, B = B_ref(fam, tab, A, D, w1, w2, i, ip)
                    tot += pr * (a ** nV * bb ** (nE + fwd) * math.exp(kap * B) - 1)
    return tot


def mc_check(fam, tab, n, trials, seed=1, copyp=0.0):
    """random path pairs; verify shared counts and the boost inequality."""
    rng = random.Random(seed)
    d, f, m, c = fam.d, fam.f, fam.m, fam.c
    bt = tab.bt
    Pw = [float(fam.Pl[fam.LN[k]]) for k in range(fam.nW)]
    worst = 0.0
    for tr in range(trials):
        paths = []
        choices = []
        for walk in range(2):
            pos = [0] * d
            V = []   # list of (point, macrostep index)
            E = set()
            ch = []
            for J in range(1, n + 1):
                if walk == 1 and rng.random() < copyp:
                    w, i = choices[0][J - 1]
                    if rng.random() < 0.3:   # perturb: same word, other forward direction / nearby word
                        i = rng.randrange(f)
                    if rng.random() < 0.3:
                        w = rng.choices(range(fam.nW), weights=Pw)[0]
                else:
                    w = rng.choices(range(fam.nW), weights=Pw)[0]
                    i = rng.randrange(f)
                ch.append((w, i))
                start = list(pos)
                prev = tuple(pos)
                V.append((prev, J))
                for k in range(1, fam.LN[w] + 1):
                    pt = list(start)
                    for q in range(m):
                        pt[f + q] = start[f + q] + int(fam.PT[w, k, q])
                    pt = tuple(pt)
                    E.add(frozenset((prev, pt))); V.append((pt, J)); prev = pt
                pos = list(prev); pos[i] += 1
                E.add(frozenset((prev, tuple(pos))))
            V.append((tuple(pos), n + 1))
            paths.append((V, E)); choices.append(ch)
        (V1, E1), (V2, E2) = paths
        S1 = {p for p, _ in V1}; S2 = {p for p, _ in V2}
        shV = len(S1 & S2); shE = len(E1 & E2)
        actual = 0.0
        for p, _ in V1:
            if p in S2: continue
            for q_, _ in V2:
                if q_ in S1: continue
                actual += bt[sum(abs(p[r] - q_[r]) for r in range(d))]
        # macrostep decomposition
        A = [0] * m; D = [0] * f
        sV = 0; sE = 0; sB = 0.0
        for J in range(n):
            (w1, i), (w2, ip) = choices[0][J], choices[1][J]
            nV, nE, fwd, B = B_ref(fam, tab, tuple(A), tuple(D), w1, w2, i, ip)
            sV += nV; sE += nE + fwd; sB += B
            A = [A[q] + int(fam.EN[w1, q]) - int(fam.EN[w2, q]) for q in range(m)]
            D = [D[q] + (1 if q == i else 0) - (1 if q == ip else 0) for q in range(f)]
        if not any(A) and not any(D):
            sV += 1
        beta_fin = bt[1] + 2 * tab.psi3[0, 0]
        assert sV == shV, (tr, sV, shV)
        assert sE == shE, (tr, sE, shE)
        ratio = actual / (sB + beta_fin)
        worst = max(worst, ratio)
        assert actual <= sB + beta_fin + 1e-15, (tr, actual, sB + beta_fin)
    return worst


if __name__ == "__main__":
    d, tkey, f, c, y = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), Fr(sys.argv[5])
    mode = sys.argv[6]
    Pt = POINTS[(d, tkey)]
    K = local_constants(d, Pt["p"], Pt["gK"], Pt["gh"], Pt["xbar"])
    fam = Family(d, f, c, y)
    tab = Tables(d, K["t"], K["qs"], c)
    if mode == "eta":
        import pickle
        R = pickle.load(open(sys.argv[7], "rb"))
        for si in [int(x) for x in sys.argv[8].split(",")]:
            Ak, Dt = R["states"][si]
            e = eta_ref(fam, tab, K, Ak, Dt)
            print(f"state {si} {Ak} {Dt}: eta_ref = {e:.12f}  kernel eta_bar = {R['eta_bar'][si]:.12f}  rel diff = {R['eta_bar'][si]/e - 1:.2e}", flush=True)
    else:
        n, trials = int(sys.argv[7]), int(sys.argv[8])
        copyp = float(sys.argv[9]) if len(sys.argv) > 9 else 0.0
        worst = mc_check(fam, tab, n, trials, copyp=copyp)
        print(f"MC check passed: {trials} path pairs of {n} macrosteps; shared counts exact; max actual/bound = {worst:.4f}")

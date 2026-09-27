#!/usr/bin/env python3
"""Noisy-cavity route -- EXACT (rational) enumeration checks of Lemma N (noisy-cavity reduction), noisy-cavity write-up Sec. 3.1.
[Packaged copy of the research script r1_lemmaN_exact.py; changes: none in the computation; an optimisation guard; exit status 1 on failure.]

(1) Torus T_4^2 (d = 2, L = 4, m = 4, n = 3; 15 cavity spins, 2^15 configurations), random off-star couplings s_off.
    LHS: R_nu(s) = E_nu[2cosh(beta s.rho)], nu = marginal on N(v) of the zero-field Gibbs law on V minus v with
         couplings s_off, by brute force over all 2^15 configurations.
    RHS: c0 * E_{nu~}[2cosh(beta_c' s.eps)], c0 = (cosh beta / cosh beta_c')^m, with the law nu~ constructed in the
         proof of Lemma N: rho_Z from its Gibbs marginal, then eps_w independent with P(eps_w = +) = (1 + m_w/b')/2,
         m_w = tanh(beta h_w), h_w = sum_{z ~ w, z != v} s_wz rho_z, b' = tanh(beta_c')/t in [b, 1], b = tanh(n beta).
    Checked: LHS == RHS exactly for all 2^m star-sign vectors s, for several rational t and b' in {b, (1+b)/2, 1};
    lambda_w in [0,1]; nu~ is flip invariant (it is, by the global spin-flip symmetry; the proof does not need it).
    Everything is rational: 2cosh(beta y) = w^{y/2} + w^{-y/2} (y even), Gibbs weight w^{#satisfied edges},
    tanh(beta h) = (w^h - 1)/(w^h + 1), e^{2 beta_c'} = (1 + t b')/(1 - t b'), cosh^2 x = 1/(1 - tanh^2 x).
(2) A synthetic graph with the star structure of Lemma N for m = 6, n = 5 (vertex v, 6 pairwise non-adjacent neighbours,
    each with exactly 5 further neighbours in Z, plus random Z-Z edges; 14 cavity spins): the same exact identity.
(3) Geometry used by Lemma N on T_L^d: for L >= 4, N(v) is an independent set and every w in N(v) has exactly 2d-1
    neighbours outside {v} u N(v) (checked for d = 2..9, L = 4..7 via the neighbourhood of the origin); for L = 3 N(v)
    is not independent (v + e_1 ~ v - e_1).
"""
import itertools, random, sys, time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
import numpy as np


def tanh_mult(t, n):
    t = Fr(t); P, Q = (1 + t) ** n, (1 - t) ** n
    return (P - Q) / (P + Q)


def check_graph(label, nodes, edges, v, t, bprime_list, rng):
    """nodes: list; edges: list of (x, y); v: the cavity vertex. Returns True iff the identity holds exactly."""
    t = Fr(t); w = (1 + t) / (1 - t)
    N = sorted({y for (x, y) in edges if x == v} | {x for (x, y) in edges if y == v})
    m = len(N)
    assert m % 2 == 0
    rest = [x for x in nodes if x != v]; pos = {x: i for i, x in enumerate(rest)}
    Z = [x for x in rest if x not in N]
    off = [(x, y) for (x, y) in edges if v not in (x, y)]
    assert not any(x in N and y in N for (x, y) in off), "N(v) must be independent"
    nb = {u: [] for u in N}
    for (x, y) in off:
        if x in N: nb[x].append(y)
        if y in N: nb[y].append(x)
    deg = {len(nb[u]) for u in N}; assert len(deg) == 1; n = deg.pop()
    assert all(z in Z for u in N for z in nb[u])
    s = {e: rng.choice([-1, 1]) for e in off}
    nR = len(rest)
    cfg = ((np.arange(2 ** nR)[:, None] >> np.arange(nR)[None, :]) & 1).astype(np.int64)
    rho = 1 - 2 * cfg                                  # bit 1 -> spin -1
    sat = np.zeros(2 ** nR, dtype=np.int64)
    for (x, y), sv in s.items():
        sat += (sv * rho[:, pos[x]] * rho[:, pos[y]] == 1)
    nE = len(off)
    Pw, Qw = w.numerator, w.denominator
    Wint = [Pw ** k * Qw ** (nE - k) for k in range(nE + 1)]   # proportional to w^{#sat}
    Nidx = [pos[u] for u in N]; Zidx = [pos[z] for z in Z]
    rhoN = rho[:, Nidx]; rhoZ = rho[:, Zidx]
    # group by (sat, rho_N) and by (rho_Z, sat) with integer counts
    keyN = (rhoN < 0).astype(np.int64) @ (1 << np.arange(m))
    keyZ = (rhoZ < 0).astype(np.int64) @ (1 << np.arange(len(Z)))
    sigmas = list(itertools.product([1, -1], repeat=m))
    f_beta = lambda y: w ** (abs(y) // 2) + 1 / w ** (abs(y) // 2)
    # LHS
    from collections import Counter
    cN = Counter(zip(sat.tolist(), keyN.tolist()))
    Ztot = sum(cnt * Wint[k] for (k, _), cnt in cN.items())
    patN = {kk: [1 - 2 * ((kk >> i) & 1) for i in range(m)] for kk in range(2 ** m)}
    LHS = {}
    for sg in sigmas:
        tot = 0
        for (k, kk), cnt in cN.items():
            y = sum(a * b for a, b in zip(sg, patN[kk]))
            tot += cnt * Wint[k] * f_beta(y)
        LHS[sg] = tot / Ztot
    # Gibbs marginal of rho_Z
    cZ = Counter(zip(keyZ.tolist(), sat.tolist()))
    PZ = {}
    for (kz, k), cnt in cZ.items():
        PZ[kz] = PZ.get(kz, 0) + cnt * Wint[k]
    ok_all = True
    b = tanh_mult(t, n)
    for bp in bprime_list(b):
        tc = t * bp; wc = (1 + tc) / (1 - tc)
        c0 = ((1 - tc * tc) / (1 - t * t)) ** (m // 2)
        f_c = lambda y: wc ** (abs(y) // 2) + 1 / wc ** (abs(y) // 2)
        RHS = {sg: Fr(0) for sg in sigmas}
        nut = {}
        lam_ok = True
        for kz, wz in PZ.items():
            zval = {Z[i]: 1 - 2 * ((kz >> i) & 1) for i in range(len(Z))}
            lam = []
            for u in N:
                h = sum((s[(u, z)] if (u, z) in s else s[(z, u)]) * zval[z] for z in nb[u])
                mw = (w ** h - 1) / (w ** h + 1)
                lw = (1 + mw / bp) / 2
                lam_ok &= 0 <= lw <= 1
                lam.append(lw)
            pz = Fr(wz, Ztot)
            for eps in itertools.product([1, -1], repeat=m):
                pe = Fr(1)
                for lw, e in zip(lam, eps):
                    pe *= lw if e == 1 else 1 - lw
                nut[eps] = nut.get(eps, Fr(0)) + pz * pe
        for sg in sigmas:
            RHS[sg] = c0 * sum(pe * f_c(sum(a * e for a, e in zip(sg, eps))) for eps, pe in nut.items())
        flip_inv = all(nut[eps] == nut[tuple(-e for e in eps)] for eps in nut)
        ident = all(LHS[sg] == RHS[sg] for sg in sigmas)
        ok = lam_ok and ident
        ok_all &= ok
        print(f"  {label}: t={t}, n={n}, b={float(b):.6f}, b'={float(bp):.6f} (beta_c'={float(np_atanh(tc)):.6f}): "
              f"R_nu == c0 R'_nu~ exactly for all {len(sigmas)} s: {ident}; lambda in [0,1]: {lam_ok}; nu~ flip-invariant: {flip_inv}",
              flush=True)
    return ok_all


def np_atanh(x):
    from math import atanh
    return atanh(float(x))


def torus(d, L):
    nodes = list(itertools.product(range(L), repeat=d))
    E = set()
    for x in nodes:
        for i in range(d):
            y = list(x); y[i] = (y[i] + 1) % L; y = tuple(y)
            E.add(tuple(sorted((x, y))))
    return nodes, sorted(E)


def geometry(dmax=9):
    ok = True
    for d in range(2, dmax + 1):
        for L in range(3, 8):
            add = lambda x, i, s: tuple((x[j] + (s if j == i else 0)) % L for j in range(d))
            o = tuple([0] * d)
            Nv = [add(o, i, s) for i in range(d) for s in (1, -1)]
            distinct = len(set(Nv)) == 2 * d
            indep = all(add(u, i, s) not in Nv for u in Nv for i in range(d) for s in (1, -1))
            outside = all(len({add(u, i, s) for i in range(d) for s in (1, -1)} - set(Nv) - {o}) == 2 * d - 1 for u in Nv)
            good = distinct and indep and outside
            if L >= 4:
                ok &= good
            elif d == 2 or d == 9:
                print(f"  L=3, d={d}: N(v) independent: {indep} (expected False)")
                ok &= not indep
    print(f"  geometry: for d=2..{dmax}, L=4..7: N(v) has 2d distinct points, is independent, and each w in N(v) has "
          f"exactly 2d-1 neighbours outside {{v}} u N(v): {ok}")
    return ok


if __name__ == "__main__":
    T0 = time.time(); rng = random.Random(20260926)
    allok = geometry()
    nodes, E = torus(2, 4)
    v = (0, 0)
    blist = lambda b: [b, (1 + b) / 2, Fr(1)]
    for t in (Fr(1, 5), Fr(3, 20), Fr(2, 5)):
        allok &= check_graph("T_4^2", nodes, E, v, t, blist, rng)
    # synthetic m = 6 star: v, w0..w5, Z = z0..z8
    for trial in range(2):
        Zs = [f"z{i}" for i in range(8)]; Ws = [f"w{i}" for i in range(6)]
        edges = [("v", u) for u in Ws]
        for u in Ws:
            for z in rng.sample(Zs, 5):
                edges.append((u, z))
        zz = [(Zs[i], Zs[j]) for i in range(8) for j in range(i + 1, 8)]
        edges += rng.sample(zz, 6)
        allok &= check_graph(f"synthetic m=6 #{trial}", ["v"] + Ws + Zs, edges, "v", Fr(3, 20), blist, rng)
    # negative control (not part of the certificate): b' = b/2 < b.  The algebraic identity survives with a SIGNED
    # 'law' (lambda_w outside [0,1]), which is exactly why Lemma N needs b' >= b, i.e. beta_c' >= atanh(t b).
    print("  negative control (b' = b/2, expected: identity True but lambda in [0,1] False):")
    check_graph("T_4^2 control", nodes, E, v, Fr(2, 5), lambda b: [b / 2], rng)
    print(f"LEMMA N EXACT CHECKS: {'ALL PASSED' if allok else 'FAILED'}  [{time.time() - T0:.0f}s]")
    sys.exit(0 if allok else 1)

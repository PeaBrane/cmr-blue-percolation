"""Exact finite-volume check of the Gibbs factorization (optional; needs NumPy).

Research-note source: theorem-A.md, Secs. 3 and 6 (Lemma 6.1, Corollary 3.3). Not a proof input.

Free-boundary n1 x n2 boxes of Z^2 with iid couplings in {0, +1, -1} (P(J != 0) = p).
Box-level isolated diamonds (non-4-cycle), box-activated = frustrated, I = internal vertices,
B = rest, G = (B, open edges within B).  Prediction (Lemma 6.1 / Prop. 3.1 finite version):
for all sites x, y with disjoint relevant-component sets R(x), R(y)
(R(x) = {C(x)} for x in B, R(x) = {C(u), C(w)} for x internal with terminals u, w),
<sigma_x sigma_y> = 0 EXACTLY.  We verify this by exact enumeration with integer arithmetic
at x = e^beta in {2, 3} (weights x^{sum J s s}; the zero test is an exact integer identity).
We also record how many predicted-zero pairs would be NONZERO if the diamond cancellation were
absent (i.e. pairs connected in the plain open graph but not in G), to show the check is not vacuous.
"""
import itertools
import random
import sys

import numpy as np

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")


def edges_of_box(n1, n2):
    V = [(i, j) for i in range(n1) for j in range(n2)]
    idx = {v: k for k, v in enumerate(V)}
    E = []
    for (i, j) in V:
        if i + 1 < n1:
            E.append((idx[(i, j)], idx[(i + 1, j)]))
        if j + 1 < n2:
            E.append((idx[(i, j)], idx[(i, j + 1)]))
    plaq = []
    for i in range(n1 - 1):
        for j in range(n2 - 1):
            c = [idx[(i, j)], idx[(i + 1, j)], idx[(i + 1, j + 1)], idx[(i, j + 1)]]  # cyclic
            plaq.append(c)
    return V, idx, E, plaq


def structure(nv, E, J, plaq):
    Jd = {}
    for (a, b), j in zip(E, J):
        Jd[(a, b)] = Jd[(b, a)] = j
    deg = [0] * nv
    for (a, b), j in zip(E, J):
        if j != 0:
            deg[a] += 1
            deg[b] += 1
    activated = []   # (internal pair, terminal pair)
    for c in plaq:
        es = [(c[k], c[(k + 1) % 4]) for k in range(4)]
        if any(Jd[e] == 0 for e in es):
            continue
        q02 = deg[c[0]] == 2 and deg[c[2]] == 2
        q13 = deg[c[1]] == 2 and deg[c[3]] == 2
        if q02 == q13:          # neither pair (not a diamond) or both (isolated 4-cycle)
            continue
        prod = 1
        for e in es:
            prod *= Jd[e]
        if prod == -1:
            internal = (c[0], c[2]) if q02 else (c[1], c[3])
            terminal = (c[1], c[3]) if q02 else (c[0], c[2])
            activated.append((internal, terminal))
    I = set(v for ip, _ in activated for v in ip)
    # components of G = (B, open edges within B)
    parent = list(range(nv))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for (a, b), j in zip(E, J):
        if j != 0 and a not in I and b not in I:
            parent[find(a)] = find(b)
    term_of = {}
    for ip, tp in activated:
        for v in ip:
            term_of[v] = tp
    R = []
    for v in range(nv):
        if v in I:
            R.append(frozenset((find(tp) for tp in term_of[v])))
        else:
            R.append(frozenset([find(v)]))
    # plain open-graph components (no cancellation) for the non-vacuity count
    parent2 = list(range(nv))

    def find2(v):
        while parent2[v] != v:
            parent2[v] = parent2[parent2[v]]
            v = parent2[v]
        return v
    for (a, b), j in zip(E, J):
        if j != 0:
            parent2[find2(a)] = find2(b)
    plain = [find2(v) for v in range(nv)]
    return R, plain, len(activated)


def exact_correlations(nv, E, J, xval):
    states = np.array(list(itertools.product((-1, 1), repeat=nv)), dtype=np.int64)
    En = np.zeros(len(states), dtype=np.int64)
    for (a, b), j in zip(E, J):
        if j != 0:
            En += j * states[:, a] * states[:, b]
    emin = int(En.min())
    levels = sorted(set(En.tolist()))
    # integer weights xval^(E - emin) (exact, python ints)
    wlev = {e: xval ** (e - emin) for e in levels}
    Z = sum(wlev[e] * int((En == e).sum()) for e in levels)
    corr_num = {}
    for i in range(nv):
        for k in range(i + 1, nv):
            prod = states[:, i] * states[:, k]
            num = 0
            for e in levels:
                m = int(prod[En == e].sum())
                if m:
                    num += wlev[e] * m
            corr_num[(i, k)] = num
    return Z, corr_num


if __name__ == '__main__':
    rng = random.Random(20260925)
    tests = 0
    zero_pred = zero_ok = nonvacuous = acts = 0
    for trial in range(80):
        n1, n2 = rng.choice([(4, 4), (4, 3), (5, 3), (3, 3)])
        V, idx, E, plaq = edges_of_box(n1, n2)
        p = rng.choice([0.6, 0.7, 0.8, 0.9])
        J = [0 if rng.random() > p else rng.choice((-1, 1)) for _ in E]
        if trial % 2 == 0:
            # plant a frustrated isolated diamond: open plaquette c, close all other edges at the
            # internal pair (c[1], c[3]), random frustrated signs
            c = rng.choice(plaq)
            es = [(c[k], c[(k + 1) % 4]) for k in range(4)]
            sg = [rng.choice((-1, 1)) for _ in range(3)]
            sg.append(-sg[0] * sg[1] * sg[2])
            for t, (a, b) in enumerate(E):
                if (a, b) in es or (b, a) in es:
                    k = [i for i, e in enumerate(es) if e in ((a, b), (b, a))][0]
                    J[t] = sg[k]
                elif a in (c[1], c[3]) or b in (c[1], c[3]):
                    J[t] = 0
            # make sure the terminals are not both of degree 2 (else it is a 4-cycle): open one more edge at c[0]
            for t, (a, b) in enumerate(E):
                if (a == c[0] or b == c[0]) and (a, b) not in es and (b, a) not in es and J[t] == 0:
                    other = b if a == c[0] else a
                    if other not in (c[1], c[3]):
                        J[t] = rng.choice((-1, 1))
                        break
        R, plain, na = structure(len(V), E, J, plaq)
        acts += na
        for xval in (2, 3):
            Z, cn = exact_correlations(len(V), E, J, xval)
            for (i, k), num in cn.items():
                if not (R[i] & R[k]):
                    zero_pred += 1
                    zero_ok += (num == 0)
                    if plain[i] == plain[k]:
                        nonvacuous += 1
            tests += 1
    print(f'boxes x temperatures tested: {tests}; activated diamonds seen: {acts}')
    print(f'pairs predicted to have <s_x s_y> = 0: {zero_pred}; exactly zero: {zero_ok}')
    print(f'of these, pairs connected in the plain open graph (zero only because of cancellation): {nonvacuous}')
    assert zero_ok == zero_pred
    assert (tests, acts, zero_pred, nonvacuous) == (160, 44, 2794, 812)   # counts displayed in the manuscript
    print('PASS')

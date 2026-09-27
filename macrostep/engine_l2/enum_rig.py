"""Rigorous enumeration of the excess kernel Q(z, .) at a state z = (A, D) of the pair chain.

For every macrostep choice omega = (w1, w2, i, i') (probability P(w1)P(w2)/f^2) at state z = (A, D) we compute
  nV   = shared lateral points (only if D = 0),   nE = shared lateral edges,
  fwd  = 1{D = 0, i = i', A + e(w1) - e(w2) = 0}  (shared forward edge),
  B    = boost bound of Lemma E (sum of four groups of b-terms, see report),
  W    = a^nV * bb^(nE+fwd) * exp(kappa B),
and accumulate mass*(W - 1) into acc[index of e(w1)-e(w2), i, i'].  Every floating operation is rounded UPWARD
(np.nextafter after round-to-nearest), all table entries are upper bounds, so acc is an upper bound of the exact
excess mass.  Returns also the largest kappa*B seen (must be <= 1 for the exp bound).
"""
import numpy as np
from numba import njit

INF = np.inf


@njit(cache=True)
def aup(x, y):
    return np.nextafter(x + y, INF)


@njit(cache=True)
def mup(x, y):
    return np.nextafter(x * y, INF)


@njit(cache=True)
def exp_up(x, coef, rem):
    """upper bound for e^x, 0 <= x <= 1: Horner on sum_{k<=K} x^k/k! (coef = upper bounds of 1/k!) + rem*x^{K+1}."""
    K = coef.shape[0] - 1
    s = coef[K]
    for k in range(K - 1, -1, -1):
        s = aup(mup(s, x), coef[k])
    xp = 1.0
    for k in range(K + 1):
        xp = mup(xp, x)
    return aup(s, mup(rem, xp))


@njit(cache=True)
def enum_state(PT, LN, EN, mass, A, D, f, pa, pb, kap, bt, Phi2, psi3, c, PAIRIDX, ndA, coef, rem):
    nW = PT.shape[0]; m = PT.shape[2]
    rmaxP = Phi2.shape[0] - 1; rmax3 = psi3.shape[0] - 1
    acc = np.zeros((ndA, f, f))
    D1 = 0
    for q in range(f):
        D1 += abs(D[q])
    # forward quantities per (i, ip)
    F1a = np.zeros((f, f), np.int64); F1b = np.zeros((f, f), np.int64)
    F2a = np.zeros((f, f), np.int64); F2b = np.zeros((f, f), np.int64)
    for i in range(f):
        for ip in range(f):
            s1 = 0; s2 = 0
            for q in range(f):
                s1 += abs(D[q] - (1 if q == ip else 0))
                s2 += abs(D[q] + (1 if q == i else 0))
            F1a[i, ip] = s1; F1b[i, ip] = s2
            m1 = 1 << 30; m2 = 1 << 30
            for qq in range(f):
                t1 = 0; t2 = 0
                for q in range(f):
                    t1 += abs(D[q] - (1 if q == ip else 0) - (1 if q == qq else 0))
                    t2 += abs(D[q] + (1 if q == i else 0) + (1 if q == qq else 0))
                if t1 < m1: m1 = t1
                if t2 < m2: m2 = t2
            F2a[i, ip] = m1; F2b[i, ip] = m2
    xmax = 0.0
    match = np.full(PT.shape[1], -1); used2 = np.zeros(PT.shape[1], np.int64)
    r1v = np.zeros(PT.shape[1], np.int64); r0v = np.zeros(PT.shape[1], np.int64)
    r1w = np.zeros(PT.shape[1], np.int64); r0w = np.zeros(PT.shape[1], np.int64)
    for w1 in range(nW):
        l1 = LN[w1]
        for w2 in range(nW):
            l2 = LN[w2]
            for i in range(l1 + 1): match[i] = -1
            for j in range(l2 + 1): used2[j] = 0
            nV = 0
            if D1 == 0:
                for i in range(l1 + 1):
                    for j in range(l2 + 1):
                        eq = True
                        for q in range(m):
                            if PT[w1, i, q] != PT[w2, j, q] - A[q]:
                                eq = False; break
                        if eq:
                            match[i] = j; used2[j] = 1; nV += 1; break
            nE = 0
            for i in range(l1):
                if match[i] >= 0 and match[i + 1] >= 0 and abs(match[i] - match[i + 1]) == 1:
                    nE += 1
            zeroA = True
            for q in range(m):
                if A[q] + EN[w1, q] - EN[w2, q] != 0:
                    zeroA = False; break
            # (0) same level unshared pairs
            B0 = 0.0
            for i in range(l1 + 1):
                if match[i] >= 0: continue
                for j in range(l2 + 1):
                    if used2[j] == 1: continue
                    dist = D1
                    for q in range(m):
                        dist += abs(PT[w1, i, q] - PT[w2, j, q] + A[q])
                    B0 = aup(B0, bt[dist])
            n1 = 0
            for i in range(l1 + 1):
                if match[i] >= 0: continue
                a1 = 0; a0 = 0
                for q in range(m):
                    a1 += abs(A[q] + PT[w1, i, q] - EN[w2, q])
                    a0 += abs(A[q] + PT[w1, i, q])
                r1v[n1] = a1; r0v[n1] = a0; n1 += 1
            n2 = 0
            for j in range(l2 + 1):
                if used2[j] == 1: continue
                a1 = 0; a0 = 0
                for q in range(m):
                    a1 += abs(PT[w2, j, q] - A[q] - EN[w1, q])
                    a0 += abs(PT[w2, j, q] - A[q])
                r1w[n2] = a1; r0w[n2] = a0; n2 += 1
            ms = mass[l1, l2]
            pidx = PAIRIDX[w1, w2]
            for i in range(f):
                for ip in range(f):
                    B = B0
                    fa = F1a[i, ip]; fb = F1b[i, ip]; ga = F2a[i, ip]; gb = F2b[i, ip]
                    for v in range(n1):
                        for k in range(c + 1):
                            B = aup(B, bt[fa + abs(r1v[v] - k)])
                        B = aup(B, Phi2[min(max(0, r1v[v] - c), rmaxP), ga])
                        B = aup(B, psi3[min(r0v[v], rmax3), D1])
                    for v in range(n2):
                        for k in range(c + 1):
                            B = aup(B, bt[fb + abs(r1w[v] - k)])
                        B = aup(B, Phi2[min(max(0, r1w[v] - c), rmaxP), gb])
                        B = aup(B, psi3[min(r0w[v], rmax3), D1])
                    x = mup(kap, B)
                    if x > xmax: xmax = x
                    fwd = 1 if (D1 == 0 and i == ip and zeroA) else 0
                    Wt = mup(mup(pa[nV], pb[nE + fwd]), exp_up(x, coef, rem))
                    ex = np.nextafter(Wt - 1.0, INF)
                    if ex <= 0.0:
                        continue
                    acc[pidx, i, ip] = aup(acc[pidx, i, ip], mup(ms, ex))
    return acc, xmax

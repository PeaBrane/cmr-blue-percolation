"""Macrostep route, independent implementation (code base R0), part 2: the one-macrostep excess kernel.

For a pair-chain state z = (A, D) (A lateral offset in Z^3, D forward offset in the root lattice of Z^f) and a
macrostep pair omega = (w, i; w', i'), macrostep write-up Sec. 3 defines
    W(z; omega) = rho^-N_V * beta_e^(N_E + fwd) * exp(kappa * B(z; omega)),   beta_e = rho/(rho_e p) >= 1,
with the shared counts N_V, N_E, fwd and the boost bound B = B0 + B1 + B2 + B3 (Lemma E).  This module returns,
for every (w, w') pair class dA = e(w) - e(w') and every (i, i'), an upper bound for
    sum_{w,w' with e(w)-e(w') = dA} P(w) P(w') f^-2 (W - 1).

Arithmetic: binary64, round to nearest, with an a-priori error analysis (macrostep write-up Sec. 7.2):
  * all inputs (tables, masses, pa, pb, kappa, Taylor coefficients) are floats >= the exact values;
  * B is a sum of <= 64 nonnegative floats, e^x (0 <= x <= 1) is a degree-22 Taylor polynomial by Horner
    (nonnegative coefficients) plus the remainder e x^23/23! < 3e-22, and W is a product of three factors;
    hence W_exact <= W_fl (1 + EPS_W) with EPS_W = 1e-12 >> 300 u;
  * the per-configuration upper bound for W - 1 is (W_fl - 1)(1 + 1e-15) + EPS_W W_fl;
  * every bin is a sum of <= 10^7 nonnegative terms; the caller multiplies each bin by (1 + INFL), INFL = 1e-8,
    which dominates gamma_{10^7} + 3u.
The structure (separable i / i' parts, per-pair records) is deliberately different from the L2 rig kernel.
"""
import numpy as np
from numba import njit

EPS_W = 1e-12
INFL = 1e-8
NTAY = 22


@njit(cache=True)
def exp_taylor(x, coef):
    # Horner for sum_{k<=NTAY} coef[k] x^k, coef[k] >= 1/k!; valid upper bound (after EPS_W) for 0 <= x <= 1
    s = coef[NTAY]
    for k in range(NTAY - 1, -1, -1):
        s = s * x + coef[k]
    return s + 3e-22


@njit(cache=True)
def l1(v):
    s = 0
    for q in range(v.shape[0]):
        s += abs(v[q])
    return s


@njit(cache=True)
def state_kernel(A, D, f, c, LEN, PTS, END, PAIRMASS, DAIDX, nDA, BT, T1, PHI2, PSI3, pa, pb, kap, coef):
    nW = LEN.shape[0]
    Fmax = T1.shape[0] - 1
    Rmax = T1.shape[1] - 1
    Dmax = PSI3.shape[1] - 1
    D1 = l1(D)
    Dc = min(D1, Dmax)
    # forward distances: walk-1 vertices vs walk-2 next / next-next starts depend on i' only;
    # walk-2 vertices vs walk-1 next / next-next starts depend on i only.
    F1a = np.zeros(f, np.int64); F2a = np.zeros(f, np.int64)
    F1b = np.zeros(f, np.int64); F2b = np.zeros(f, np.int64)
    tmp = np.zeros(f, np.int64)
    for j in range(f):
        for q in range(f):
            tmp[q] = D[q]
        tmp[j] -= 1
        F1a[j] = l1(tmp)
        best = 1 << 40
        for a in range(f):
            tmp[a] -= 1
            v = l1(tmp)
            if v < best:
                best = v
            tmp[a] += 1
        F2a[j] = best
        for q in range(f):
            tmp[q] = D[q]
        tmp[j] += 1
        F1b[j] = l1(tmp)
        best = 1 << 40
        for a in range(f):
            tmp[a] += 1
            v = l1(tmp)
            if v < best:
                best = v
            tmp[a] -= 1
        F2b[j] = best
    acc = np.zeros((nDA, f, f))
    xmax = 0.0
    sh1 = np.zeros(c + 1, np.int64)   # partner index of walk-1 point k, or -1
    sh2 = np.zeros(c + 1, np.int64)
    X1 = np.zeros(f); X2 = np.zeros(f)
    for w in range(nW):
        lw = LEN[w]
        for wp in range(nW):
            lwp = LEN[wp]
            for k in range(lw + 1):
                sh1[k] = -1
            for k in range(lwp + 1):
                sh2[k] = -1
            nV = 0
            if D1 == 0:
                for k in range(lw + 1):
                    for kp in range(lwp + 1):
                        if (A[0] + PTS[w, k, 0] == PTS[wp, kp, 0] and A[1] + PTS[w, k, 1] == PTS[wp, kp, 1]
                                and A[2] + PTS[w, k, 2] == PTS[wp, kp, 2]):
                            sh1[k] = kp
                            sh2[kp] = k
                            nV += 1
            nE = 0
            for k in range(1, lw + 1):
                if sh1[k - 1] >= 0 and sh1[k] >= 0 and abs(sh1[k] - sh1[k - 1]) == 1:
                    nE += 1
            zeroA = (A[0] + END[w, 0] - END[wp, 0] == 0 and A[1] + END[w, 1] - END[wp, 1] == 0
                     and A[2] + END[w, 2] - END[wp, 2] == 0)
            # B0 (same level, unshared x unshared) and B3 (psi3 at each unshared vertex)
            B03 = 0.0
            for k in range(lw + 1):
                if sh1[k] >= 0:
                    continue
                r0 = abs(A[0] + PTS[w, k, 0]) + abs(A[1] + PTS[w, k, 1]) + abs(A[2] + PTS[w, k, 2])
                B03 += PSI3[min(r0, Rmax), Dc]
                for kp in range(lwp + 1):
                    if sh2[kp] >= 0:
                        continue
                    dist = D1 + (abs(A[0] + PTS[w, k, 0] - PTS[wp, kp, 0]) + abs(A[1] + PTS[w, k, 1] - PTS[wp, kp, 1])
                                 + abs(A[2] + PTS[w, k, 2] - PTS[wp, kp, 2]))
                    B03 += BT[dist]
            for kp in range(lwp + 1):
                if sh2[kp] >= 0:
                    continue
                r0 = abs(PTS[wp, kp, 0] - A[0]) + abs(PTS[wp, kp, 1] - A[1]) + abs(PTS[wp, kp, 2] - A[2])
                B03 += PSI3[min(r0, Rmax), Dc]
            # B1 + B2, walk-1 vertices (depend on i'), walk-2 vertices (depend on i)
            for j in range(f):
                X1[j] = 0.0
                X2[j] = 0.0
            for k in range(lw + 1):
                if sh1[k] >= 0:
                    continue
                r1 = (abs(A[0] + PTS[w, k, 0] - END[wp, 0]) + abs(A[1] + PTS[w, k, 1] - END[wp, 1])
                      + abs(A[2] + PTS[w, k, 2] - END[wp, 2]))
                r2 = r1 - c
                if r2 < 0:
                    r2 = 0
                for j in range(f):
                    X1[j] += T1[min(F1a[j], Fmax), min(r1, Rmax)] + PHI2[min(r2, Rmax), min(F2a[j], Fmax)]
            for kp in range(lwp + 1):
                if sh2[kp] >= 0:
                    continue
                r1 = (abs(PTS[wp, kp, 0] - A[0] - END[w, 0]) + abs(PTS[wp, kp, 1] - A[1] - END[w, 1])
                      + abs(PTS[wp, kp, 2] - A[2] - END[w, 2]))
                r2 = r1 - c
                if r2 < 0:
                    r2 = 0
                for j in range(f):
                    X2[j] += T1[min(F1b[j], Fmax), min(r1, Rmax)] + PHI2[min(r2, Rmax), min(F2b[j], Fmax)]
            ms = PAIRMASS[w, wp]
            di = DAIDX[w, wp]
            for i in range(f):
                for ip in range(f):
                    B = B03 + X1[ip] + X2[i]
                    x = kap * B
                    if x > xmax:
                        xmax = x
                    fw = 1 if (D1 == 0 and zeroA and i == ip) else 0
                    Wf = pa[nV] * pb[nE + fw] * exp_taylor(x, coef)
                    ex = (Wf - 1.0) * (1.0 + 1e-15) + EPS_W * Wf
                    acc[di, i, ip] += ms * ex
    return acc, xmax

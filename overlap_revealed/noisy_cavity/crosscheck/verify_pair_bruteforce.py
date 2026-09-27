#!/usr/bin/env python3
"""Noisy-cavity route -- third, independent check of the pair certificates (Lemma P, noisy-cavity write-up Sec. 3.2).
[Packaged copy of the research script r1_verify_pair.py; changes: the default data path is ../data/CERT_DATA_d9.json relative to this file; an optimisation guard; exit status 1 on failure.]

Written from the statement of Lemma P only (no code shared with code base L1's exact_pair.py or verify_pair_indep.py).
Input: data/CERT_DATA_d9.json (t, p, wc', psi_U(x) as exact rational strings).  All decisions are exact
(fractions.Fraction and Python integers); numpy is used only for integer bit-counting in the enumeration of (C).

Checks, for each temperature in the file:
 (0) model constants:  w = (1+t)/(1-t), a = w^2/(1+w^2), p_A = 1-w^-2, tanh(2 beta) = p_A a = 2t/(1+t^2),
     r = p/p_A, lam = w^2 (1-r)/r;  p_A/2 < p < tanh(2 beta);  alpha_U = w^U - lam w^-U > 0 and
     gamma_U = lam w^U - w^-U > 0 for U = 1,3,...,n (n = 2d-1);  noisy constant e^{2 beta_c} = (1+t b)/(1-t b) with
     b = tanh(n beta) = ((1+t)^n-(1-t)^n)/((1+t)^n+(1-t)^n), and e^{2 beta_c} <= wc' <= w.
 (V) point conditions psi_U(x) > phi_U(p_x) and chord conditions P(theta) > 0 on [0,1] for all x != y, where
     P = L A^2 B^2 - alpha_U B^2 + gamma_U A^2 (degree 5), decided by an own exact Sturm-sequence root count
     (P(0) > 0, P(1) > 0 and no root in (0,1)).
 (C) Psi(j) = sum_{s' in {+-1}^n, U(s') >= 1} psi_{U(s')}(s'.rho'_j) <= 0 for rho'_j = (+^j, -^(n-j)), j = 0..n,
     by BRUTE-FORCE enumeration of all 2^n vectors s' (integer histogram of (U, x)), not by a hypergeometric formula.
 Also printed: the noisy aligned-frozen witness p_A r'(aligned) (exact; an upper bound on the class infimum).
"""
import json, os, sys, time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
from math import comb
import numpy as np


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def poly_mul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def poly_rem(a, b):
    a, b = trim(a), trim(b)
    assert not (len(b) == 1 and b[0] == 0)
    db = len(b) - 1
    while len(a) - 1 >= db and not (len(a) == 1 and a[0] == 0):
        coef = a[-1] / b[-1]; sh = len(a) - 1 - db
        for i in range(len(b)):
            a[i + sh] -= coef * b[i]
        assert a[-1] == 0
        a.pop()
        if not a:
            return [Fr(0)]
        a = trim(a)
    return a


def poly_eval(p, x):
    s = Fr(0)
    for c in reversed(p):
        s = s * x + c
    return s


def sturm_count_open01(P):
    """number of distinct real roots of P in the open interval (0,1); requires P(0) != 0 != P(1)."""
    P = trim(P)
    seq = [P, trim([i * P[i] for i in range(1, len(P))] or [Fr(0)])]
    while not (len(seq[-1]) == 1 and seq[-1][0] == 0):
        r = poly_rem(seq[-2], seq[-1])
        if len(r) == 1 and r[0] == 0:
            break
        seq.append([-c for c in r])

    def var(x):
        vals = [v for v in (poly_eval(p, x) for p in seq) if v != 0]
        return sum(1 for u, v in zip(vals, vals[1:]) if (u > 0) != (v > 0))
    return var(Fr(0)) - var(Fr(1))


def popcount(v):
    v = v.copy(); c = np.zeros_like(v)
    while np.any(v):
        c += v & 1
        v >>= 1
    return c


def check(key, E, d=9):
    t0 = time.time(); n = 2 * d - 1
    t = Fr(E['t']); p = Fr(E['p']); wc = Fr(E['wc_prime'])
    psi = {int(U): {int(x): Fr(v) for x, v in row.items()} for U, row in E['psi'].items()}
    w = (1 + t) / (1 - t); a = w * w / (1 + w * w); pA = 1 - 1 / w ** 2; th2b = pA * a
    assert th2b == 2 * t / (1 + t * t)
    r = p / pA; lam = w * w * (1 - r) / r
    Us = list(range(1, n + 1, 2)); X = list(range(-n, n + 1, 2))
    al = {U: w ** U - lam / w ** U for U in Us}; ga = {U: lam * w ** U - 1 / w ** U for U in Us}
    N1, N2 = (1 + t) ** n, (1 - t) ** n; b = (N1 - N2) / (N1 + N2)
    wcx = (1 + t * b) / (1 - t * b)
    checks0 = dict(p_range=pA / 2 < p < th2b, alpha_gamma_pos=all(al[U] > 0 and ga[U] > 0 for U in Us),
                   wc_range=wcx <= wc <= w, data_shape=set(psi) == set(Us) and all(set(psi[U]) == set(X) for U in Us))
    ok0 = all(checks0.values())
    f = lambda y: wc ** (abs(y) // 2) + 1 / wc ** (abs(y) // 2)       # 2cosh(beta_c' y) for even y
    # (V)
    nfail = 0; npol = 0; npt = 0
    for U in Us:
        for x in X:
            npt += 1
            if not psi[U][x] > al[U] / f(x - 1) ** 2 - ga[U] / f(x + 1) ** 2:
                nfail += 1
        for i in range(len(X)):
            for j in range(i + 1, len(X)):
                x, y = X[i], X[j]
                A = [f(y - 1), f(x - 1) - f(y - 1)]          # first coordinate of theta p_x + (1-theta) p_y
                B = [f(y + 1), f(x + 1) - f(y + 1)]          # second coordinate
                L = [psi[U][y], psi[U][x] - psi[U][y]]        # theta psi(x) + (1-theta) psi(y)
                A2 = poly_mul(A, A); B2 = poly_mul(B, B)
                P = poly_mul(poly_mul(L, A2), B2)
                for k, c in enumerate(B2):
                    P[k] -= al[U] * c
                for k, c in enumerate(A2):
                    P[k] += ga[U] * c
                npol += 1
                if not (poly_eval(P, Fr(0)) > 0 and poly_eval(P, Fr(1)) > 0 and sturm_count_open01(P) == 0):
                    nfail += 1
    # (C): brute force over s' in {+-1}^n; bit i of the integer s is 1 iff s'_i = -1
    s = np.arange(2 ** n, dtype=np.int64)
    tot_minus = popcount(s); Uv = n - 2 * tot_minus; sel = Uv >= 1
    Psi = []
    for j in range(n + 1):
        am = popcount(s & ((1 << j) - 1))           # minus signs of s' on the j plus-coordinates of rho'_j
        bm = tot_minus - am                          # minus signs of s' on the n-j minus-coordinates
        xv = (j - 2 * am) - ((n - j) - 2 * bm)       # s'.rho'_j
        keys = (Uv[sel] - 1) // 2 * (2 * n + 1) + (xv[sel] + n)
        cnt = np.bincount(keys, minlength=((n + 1) // 2) * (2 * n + 1))
        assert int(cnt.sum()) == 2 ** (n - 1)
        tot = Fr(0)
        for kk in np.nonzero(cnt)[0]:
            U = 2 * int(kk // (2 * n + 1)) + 1; x = int(kk % (2 * n + 1)) - n
            tot += int(cnt[kk]) * psi[U][x]
        Psi.append(tot)
    okC = max(Psi) <= 0
    # noisy aligned-frozen witness
    Xw = sum(comb(n, k) * a ** k * (1 - a) ** (n - k) / f(2 * k - n + 1) ** 2 for k in range(n + 1))
    Yw = sum(comb(n, k) * a ** k * (1 - a) ** (n - k) / f(2 * k - n - 1) ** 2 for k in range(n + 1))
    wit = pA * a * Xw / (a * Xw + (1 - a) * Yw)
    ok = ok0 and nfail == 0 and okC and p <= wit
    print(f"{key}: p={p}; noisy aligned witness p_A r'(al) = {float(wit):.9f}; p/witness = {float(p / wit):.7f}; "
          f"checks {checks0}", flush=True)
    print(f"   (V) {npt} point tests + {npol} chord polynomials (own Sturm): {nfail} failures", flush=True)
    print(f"   (C) brute force over 2^{n} vectors s': max_j Psi(j) = {float(max(Psi)):+.6e} (at j={Psi.index(max(Psi))}); "
          f"all <= 0: {okC}  [{time.time() - t0:.0f}s]", flush=True)
    print(f"   ==> {'PAIR CERTIFICATE VERIFIED (noisy-cavity implementation)' if ok else 'NOT VERIFIED'}", flush=True)
    return ok, Psi, wit


if __name__ == "__main__":
    J = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'CERT_DATA_d9.json')))
    sel = sys.argv[2:]
    allok = True
    for key, E in J.items():
        if sel and key not in sel:
            continue
        ok, _, _ = check(key, E)
        allok &= ok
    print("ALL SELECTED PAIR CERTIFICATES VERIFIED" if allok else "SOME PAIR CERTIFICATE FAILED")
    sys.exit(0 if allok else 1)

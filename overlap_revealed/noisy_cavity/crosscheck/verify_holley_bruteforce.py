#!/usr/bin/env python3
"""Noisy-cavity route -- independent check of the Holley-line certificates (Lemma 4.20 with Lemma 4.15; condition (N2)
of Section 4.7 of the manuscript).  Written from the statement of Lemma 4.20 only, independently of the other
implementations of this certificate.  The default data path is ../data/CERT_DATA_d9.json relative to this file; Python
refuses optimized mode (-O), and the exit status is 1 on failure.  Input: data/CERT_DATA_d9.json (t, wc', g_K, g_h, and for every environment k the stored
Lambda'_k and the tangent points c(i+, i-)).  All decisions are exact (fractions.Fraction, Python integers); numpy is
used only for integer counting.

For environment k (N+ = k coordinates, N- = m-k coordinates, S = 2k-m), with a = w^2/(1+w^2), C = (w+1/w)/2:
   lambda+(s) = a^{i+}(1-a)^{k-i+} 2^{-(m-k)},   lambda-(s) = a^{i-}(1-a)^{m-k-i-} 2^{-k},   i+- = #plus of s on N+-,
   Lambda'_k  = g_h g_K^S C^{-S}   (recomputed and compared with the stored value),
   kappa = lambda+ - Lambda' lambda-,   kappa_s(i+, i-) = (kappa(i+, i-) + kappa(k-i+, m-k-i-))/2,
   Gamma(l+, l-) = sum_s kappa_s(s) * [kappa_s > 0 : 3 c^-2 - 2 f(s.rho) c^-3 ;  kappa_s < 0 : f(s.rho)^-2],
   f(y) = 2cosh(beta_c' y) = wc'^{y/2} + wc'^{-y/2} (y even),  rho = (+^{l+} -^{k-l+} | +^{l-} -^{m-k-l-}).
The counts #{s : i+, i-, s.rho = D} are obtained by BRUTE-FORCE enumeration of the 2^k and 2^{m-k} block patterns
(integer histograms) and an integer convolution -- not by binomial/hypergeometric formulas.
Certified: Gamma(l+, l-) >= 0 for every class and every k (then Phi+/Phi- >= g_h g_K^S for every cavity law).
Also checked: e^{2beta_c} <= wc' <= w; c > 0 wherever used; g_h < 1, g_K >= 1, g_K g_h <= 1, g_K^d < 2718/1000,
g_K^{2d} < 2718/1000; the flip symmetry Gamma(l+, l-) = Gamma(k-l+, m-k-l-) (a consistency check).
"""
import json, os, sys, time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
import numpy as np


def popcount(v):
    v = v.copy(); c = np.zeros_like(v)
    while np.any(v):
        c += v & 1
        v >>= 1
    return c


_cache = {}
def block_counts(nb, l):
    """cnt[i, (dot+nb)//2] = #{s in {+-1}^nb : #plus(s) = i, s.rho = dot} with rho = (+^l, -^(nb-l))."""
    if (nb, l) not in _cache:
        s = np.arange(2 ** nb, dtype=np.int64)
        tm = popcount(s); am = popcount(s & ((1 << l) - 1)); bm = tm - am
        i = nb - tm
        dot = (l - 2 * am) - ((nb - l) - 2 * bm)
        cnt = np.zeros((nb + 1, nb + 1), dtype=np.int64)
        np.add.at(cnt, (i, (dot + nb) // 2), 1)
        assert int(cnt.sum()) == 2 ** nb
        _cache[(nb, l)] = cnt
    return _cache[(nb, l)]


def check(key, E, d=9, verbose=True):
    T0 = time.time(); m = 2 * d; n = m - 1
    t = Fr(E['t']); wc = Fr(E['wc_prime']); H = E['holley']; gK = Fr(H['gK']); gh = Fr(H['gh'])
    w = (1 + t) / (1 - t); a = w * w / (1 + w * w); C = (w + 1 / w) / 2
    N1, N2 = (1 + t) ** n, (1 - t) ** n; b = (N1 - N2) / (N1 + N2); wcx = (1 + t * b) / (1 - t * b)
    hyp = dict(wc_range=wcx <= wc <= w, gh_lt_1=gh < 1, gK_ge_1=gK >= 1, H1_gKgh=gK * gh <= 1,
               H2_gKd=gK ** d < Fr(2718, 1000), gK2d=gK ** (2 * d) < Fr(2718, 1000))
    f = {D: wc ** (abs(D) // 2) + 1 / wc ** (abs(D) // 2) for D in range(-m, m + 1, 2)}
    envs = {e['k']: e for e in H['envs']}
    assert sorted(envs) == list(range(m + 1))
    allok = all(hyp.values()); rows = []
    for k in range(m + 1):
        t0 = time.time(); kp, km = k, m - k; S = kp - km
        Lam = gh * gK ** S / C ** S
        e = envs[k]; lam_ok = Fr(e['Lam']) == Lam
        c = {tuple(map(int, kk.split(','))): Fr(v) for kk, v in e['c'].items()}
        lp_ = lambda ip: a ** ip * (1 - a) ** (kp - ip) / 2 ** km
        lm_ = lambda im: a ** im * (1 - a) ** (km - im) / 2 ** kp
        kap = lambda ip, im: lp_(ip) - Lam * lm_(im)
        F = {}; c_ok = True
        for ip in range(kp + 1):
            for im in range(km + 1):
                ks = (kap(ip, im) + kap(kp - ip, km - im)) / 2
                if ks > 0:
                    cc = c[(ip, im)]; c_ok &= cc > 0
                    F[(ip, im)] = [ks * (3 / cc ** 2 - 2 * f[2 * j - m] / cc ** 3) for j in range(m + 1)]
                elif ks < 0:
                    F[(ip, im)] = [ks / f[2 * j - m] ** 2 for j in range(m + 1)]
        G = {}
        for lp in range(kp + 1):
            cP = block_counts(kp, lp)
            for lm in range(km + 1):
                cM = block_counts(km, lm)
                tot = Fr(0)
                for (ip, im), Fv in F.items():
                    Hc = np.convolve(cP[ip], cM[im])            # integer counts of D = 2j - m, j = 0..m
                    for j in np.nonzero(Hc)[0]:
                        tot += int(Hc[j]) * Fv[j]
                G[(lp, lm)] = tot
        sym_ok = all(G[(lp, lm)] == G[(kp - lp, km - lm)] for (lp, lm) in G)
        gmin = min(G.values()); arg = min(G, key=G.get)
        ok = lam_ok and c_ok and gmin >= 0
        allok &= ok and sym_ok
        rows.append((k, S, gmin, arg, len(G)))
        if verbose:
            print(f"  k={k:2d} S={S:+3d}: {len(G):3d} classes, min Gamma = {float(gmin):+.6e} at (l+,l-)={arg}; "
                  f"Lambda' matches: {lam_ok}; c>0: {c_ok}; flip symmetry: {sym_ok} -> {'OK' if ok else 'FAIL'} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
    print(f"{key}: line g_K={gK} g_h={gh}, wc'={wc}; hypotheses {hyp}", flush=True)
    print(f"   ==> {'HOLLEY LINE VERIFIED in all ' + str(m + 1) + ' environments (noisy-cavity implementation)' if allok else 'NOT VERIFIED'}"
          f"; overall min Gamma = {float(min(r[2] for r in rows)):+.6e}  [{time.time() - T0:.0f}s]", flush=True)
    return allok, rows


if __name__ == "__main__":
    J = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'CERT_DATA_d9.json')))
    sel = sys.argv[2:]
    allok = True
    for key, E in J.items():
        if sel and key not in sel:
            continue
        ok, _ = check(key, E)
        allok &= ok
    print("ALL SELECTED HOLLEY LINES VERIFIED" if allok else "SOME HOLLEY LINE FAILED")
    sys.exit(0 if allok else 1)

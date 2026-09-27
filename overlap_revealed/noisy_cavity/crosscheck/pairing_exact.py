#!/usr/bin/env python3
"""Noisy-cavity route -- EXACT checks of the pairing identity (Lemma P, step (i)) and of the full inequality chain of Lemma P.
[Packaged copy of the research script r1_pairing_exact.py; changes: the data path is ../data/CERT_DATA_d9.json relative to this file; an optimisation guard; exit status 1 on failure.]

Part A (small n, exact, all laws random and flip invariant): for n = 2d-1 in {3, 5, 7}, m = n+1,
   r(nu; [m]) from its DEFINITION (sum over all s in {+-1}^m with the tilted product law), equals aX/(aX+(1-a)Y);
   B(s') := E_{nu'} f(s'.rho'+1) = R(+, s') and A(s') := E_{nu'} f(s'.rho'-1) = R(+, -s')  (nu' = nu( . | rho_u=+));
   G := Y - lam X  ==  (w/(1+w^2))^n * sum_{U(s') >= 1} [alpha_U A(s')^-2 - gamma_U B(s')^-2]   (note (a(1-a))^{1/2} = w/(1+w^2));
   and  r(nu;[m]) >= p/p_A  <=>  G <= 0.
Part B (d = 9, exact): for random sparse laws nu' on {+-1}^17 (1..4 atoms, random rational weights, plus a few
   structured ones: aligned, one-off layer, balanced), with the CERTIFIED psi of data/CERT_DATA_d9.json:
   G(nu') <= (w/(1+w^2))^17 E_{nu'} Psi(rho') <= 0, and r(nu') >= p/p_A.  (A spot check of the proved chain; the
   proof does not depend on it.)  s' is enumerated exhaustively (2^17 vectors), grouped by the integer key
   (U, s'.rho'_1, ..., s'.rho'_K) before the exact rational arithmetic.
"""
import itertools, json, os, random, sys, time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
from math import comb
import numpy as np


def model(t, wc, p, n):
    w = (1 + t) / (1 - t); a = w * w / (1 + w * w); pA = 1 - 1 / w ** 2
    r = p / pA; lam = w * w * (1 - r) / r
    f = lambda y: wc ** (abs(y) // 2) + 1 / wc ** (abs(y) // 2)
    al = {U: w ** U - lam / w ** U for U in range(1, n + 1, 2)}
    ga = {U: lam * w ** U - 1 / w ** U for U in range(1, n + 1, 2)}
    return w, a, pA, r, lam, f, al, ga


def part_a(rng):
    ok_all = True
    for (d, t, wc) in ((2, Fr(1, 3), Fr(3, 2)), (3, Fr(1, 5), Fr(13, 10)), (4, Fr(3, 20), Fr(6, 5)), (3, Fr(2, 5), Fr(2))):
        n = 2 * d - 1; m = n + 1
        w0 = (1 + t) / (1 - t); pA0 = 1 - 1 / w0 ** 2; th2b = 2 * t / (1 + t * t)
        for trial in range(6):
            p = pA0 / 2 + (th2b - pA0 / 2) * Fr(rng.randint(1, 99), 100)
            w, a, pA, r, lam, f, al, ga = model(t, wc, p, n)
            assert (a * (1 - a)) == (w / (1 + w * w)) ** 2
            pats = list(itertools.product([1, -1], repeat=m))
            K = [1, 2, 5, 2 ** m][trial % 4]
            atoms = rng.sample(pats, min(K, len(pats)))
            nu = {}
            for rho in atoms:
                q = Fr(rng.randint(1, 30))
                nu[rho] = nu.get(rho, 0) + q; neg = tuple(-x for x in rho); nu[neg] = nu.get(neg, 0) + q
            Z = sum(nu.values()); nu = {k: v / Z for k, v in nu.items()}
            R = {s: sum(pr * f(sum(x * y for x, y in zip(s, rho))) for rho, pr in nu.items()) for s in pats}
            mu = lambda s: a ** sum(1 for x in s if x > 0) * (1 - a) ** sum(1 for x in s if x < 0)
            num = sum(mu(s) / R[s] ** 2 for s in pats if s[0] == 1); den = sum(mu(s) / R[s] ** 2 for s in pats)
            r_def = num / den
            sp = list(itertools.product([1, -1], repeat=n))
            X = sum(mu(s) / R[(1,) + s] ** 2 for s in sp); Y = sum(mu(s) / R[(-1,) + s] ** 2 for s in sp)
            r_xy = a * X / (a * X + (1 - a) * Y)
            nup = {rho[1:]: 2 * pr for rho, pr in nu.items() if rho[0] == 1}
            assert sum(nup.values()) == 1
            A = {s: sum(pr * f(sum(x * y for x, y in zip(s, rp)) - 1) for rp, pr in nup.items()) for s in sp}
            B = {s: sum(pr * f(sum(x * y for x, y in zip(s, rp)) + 1) for rp, pr in nup.items()) for s in sp}
            ab_ok = all(B[s] == R[(1,) + s] and A[s] == R[(1,) + tuple(-x for x in s)] for s in sp)
            G = Y - lam * X
            paired = (w / (1 + w * w)) ** n * sum(al[sum(s)] / A[s] ** 2 - ga[sum(s)] / B[s] ** 2 for s in sp if sum(s) >= 1)
            eq = (r_def == r_xy) and ab_ok and (G == paired) and ((r_def >= r) == (G <= 0))
            ok_all &= eq
            print(f"  A: n={n} t={t} wc'={wc} p={float(p):.5f} |supp nu|={len(nu):3d}: r_def=r_XY: {r_def == r_xy}; "
                  f"A,B = R(+,-s'),R(+,s'): {ab_ok}; G == paired form: {G == paired}; "
                  f"[r>=p/p_A] == [G<=0]: {(r_def >= r) == (G <= 0)} (G={float(G):+.3e})", flush=True)
    return ok_all


def part_b(E, rng, ntrials=24, d=9):
    n = 2 * d - 1
    t = Fr(E['t']); p = Fr(E['p']); wc = Fr(E['wc_prime'])
    psi = {int(U): {int(x): Fr(v) for x, v in row.items()} for U, row in E['psi'].items()}
    w, a, pA, r, lam, f, al, ga = model(t, wc, p, n)
    # Psi(j) by the hypergeometric/multiplicity formula (as in the proof)
    Psi = []
    for j in range(n + 1):
        tot = Fr(0)
        for U in range(1, n + 1, 2):
            kp = (n + U) // 2
            for h in range(max(0, kp - j), min(kp, n - j) + 1):
                tot += comb(n - j, h) * comb(j, kp - h) * psi[U][U - 2 * (2 * h - (n - j))]
        Psi.append(tot)
    c0 = (w / (1 + w * w)) ** n
    s = np.arange(2 ** n, dtype=np.int64)
    bits = ((s[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
    S = 1 - 2 * bits                                   # s' vectors, shape (2^n, n)
    U = S.sum(1)
    laws = []
    al_vec = tuple([1] * n)
    laws.append(("aligned", {al_vec: Fr(1)}))
    laws.append(("one-off layer", {tuple(-1 if i == k else 1 for i in range(n)): Fr(1, n) for k in range(n)}))
    laws.append(("balanced frozen j=9", {tuple([1] * 9 + [-1] * 8): Fr(1)}))
    for tr in range(ntrials):
        K = rng.randint(1, 4); nu = {}
        for _ in range(K):
            j = rng.randint(0, n)
            rho = tuple(rng.sample([1] * j + [-1] * (n - j), n))
            nu[rho] = nu.get(rho, 0) + Fr(rng.randint(1, 20))
        Z = sum(nu.values()); laws.append((f"random K={K}", {k: v / Z for k, v in nu.items()}))
    ok_all = True; worst = None; minr = None
    for name, nu in laws:
        rhos = list(nu); P = np.array(rhos, dtype=np.int64)       # (K, n)
        dots = S @ P.T                                            # (2^n, K): s'.rho'_i
        keys = np.concatenate([U[:, None], dots], axis=1)
        uniq, cnt = np.unique(keys, axis=0, return_counts=True)
        mu_of = lambda Uv: a ** ((n + Uv) // 2) * (1 - a) ** ((n - Uv) // 2)
        X = Fr(0); Y = Fr(0)
        for row, c in zip(uniq.tolist(), cnt.tolist()):
            Uv, xs = row[0], row[1:]
            Bv = sum(nu[rh] * f(x + 1) for rh, x in zip(rhos, xs))     # R(+, s')
            X += c * mu_of(Uv) / Bv ** 2
            Y += c * mu_of(-Uv) / Bv ** 2                            # Y = sum_s' mu'(s') R(+,-s')^-2 = sum mu'(-s') R(+,s')^-2
        G = Y - lam * X
        bound = c0 * sum(pr * Psi[sum(1 for x in rh if x > 0)] for rh, pr in nu.items())
        rv = a * X / (a * X + (1 - a) * Y)
        ok = G <= bound <= 0 and rv >= r
        ok_all &= ok
        worst = (G - bound) if worst is None else max(worst, G - bound)
        minr = rv if minr is None else min(minr, rv)
        print(f"  B t={t}: {name:20s}: G={float(G):+.4e} <= c0 E Psi = {float(bound):+.4e} <= 0: {G <= bound <= 0}; "
              f"p_A r(nu') = {float(pA * rv):.6f} >= p = {float(p):.4f}: {rv >= r}", flush=True)
    print(f"  B t={t}: {len(laws)} laws, all chain inequalities hold exactly: {ok_all}; max(G - bound) = {float(worst):+.3e}; "
          f"min p_A r = {float(pA * minr):.6f}", flush=True)
    return ok_all


if __name__ == "__main__":
    T0 = time.time(); rng = random.Random(7)
    ok = part_a(rng)
    J = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'CERT_DATA_d9.json')))
    for key in (sys.argv[1:] or ['t=3/20', 't=31/200']):
        ok &= part_b(J[key], rng)
    print(f"PAIRING IDENTITY / CHAIN EXACT CHECKS: {'ALL PASSED' if ok else 'FAILED'}  [{time.time() - T0:.0f}s]")
    sys.exit(0 if ok else 1)

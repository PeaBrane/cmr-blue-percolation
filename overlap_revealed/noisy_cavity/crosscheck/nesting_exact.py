#!/usr/bin/env python3
"""Noisy-cavity route -- exact checks for the nesting of the noisy classes (Remark 4.17 of the manuscript,
Section 4.4). ROOT resolves to ../ , which holds data/; Python refuses optimized mode (-O), and the exit status is 1
on failure.

Nesting.  Let 0 < beta1 <= beta2, t_i = tanh(beta_i), b'' = t1/t2 in (0, 1].  For a law nu of eps on {+-1}^m let
xi = (xi_w) be iid +-1 with E xi_w = b'', independent of eps, and let nu'' be the law of rho = (eps_w xi_w)_w.  Then
    R^{beta2}_{nu''}(s) = c * R^{beta1}_{nu}(s)   for all s in {+-1}^m,   c = (cosh beta2 / cosh beta1)^m,
and nu'' is flip invariant when nu is.  Here R^{beta}_nu(s) = E_nu[2 cosh(beta s.rho)].

(1) The identity is checked as an equality of rationals.  For even m, s.rho is even, 2cosh(beta y) = w^{y/2} + w^{-y/2}
    with w = e^{2 beta} = (1+t)/(1-t), and c = ((1 - t1^2)/(1 - t2^2))^{m/2} (cosh^2 = 1/(1 - tanh^2)).
    Grid: m in {2, 4, 6}; t2 in {1/5, 3/20, 2/5}; b'' in {1/3, 1/2, 9/10, 1}; 3 random flip-invariant laws nu per cell
    (random rational weights on random supports, symmetrised).  Also checked: nu'' is a probability law and is flip
    invariant.
(2) Negative control: for b'' > 1 (i.e. beta1 > beta2) the noise law P(xi = +1) = (1 + b'')/2 exceeds 1.
(3) Strictness at the five certified points (the strict comparison of Remark 4.17): the plain aligned-frozen witness
    p_A r_beta(aligned; [m]) (beta inside R, an exact upper bound on p_A * inf over the PLAIN class) is < p, the certified
    noisy floor.  Also printed: the noisy aligned-frozen witness at beta_c and at the certificate's beta' (w_c'), both
    >= p.  All values exact rationals; witnesses (upper bounds) shown rounded UP to 9 digits, ratios DOWN.
"""
import itertools, json, os, random, sys, time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)


def two_cosh(w, y):
    assert y % 2 == 0
    k = abs(y) // 2
    return w ** k + 1 / w ** k


def R(w, law, m):
    """R^{beta}_law(s) for all s, w = e^{2 beta}; law: dict rho(tuple) -> prob."""
    out = {}
    for s in itertools.product((1, -1), repeat=m):
        out[s] = sum(pr * two_cosh(w, sum(a * b for a, b in zip(s, rho))) for rho, pr in law.items())
    return out


def random_flip_invariant_law(m, rng):
    cube = list(itertools.product((1, -1), repeat=m))
    k = rng.randint(1, min(6, len(cube)))
    supp = rng.sample(cube, k)
    wts = [Fr(rng.randint(1, 50)) for _ in supp]
    tot = sum(wts)
    law = {}
    for rho, wt in zip(supp, wts):
        law[rho] = law.get(rho, Fr(0)) + wt / (2 * tot)
        neg = tuple(-x for x in rho)
        law[neg] = law.get(neg, Fr(0)) + wt / (2 * tot)
    return law


def noise(law, m, bpp):
    """law of rho = eps * xi, xi iid with P(xi=+1) = (1+bpp)/2."""
    pp = (1 + bpp) / 2; pm = 1 - pp
    assert 0 <= pm <= pp <= 1
    out = {}
    for eps, pr in law.items():
        for xi in itertools.product((1, -1), repeat=m):
            wt = pr
            for x in xi:
                wt *= pp if x == 1 else pm
            if wt == 0:
                continue
            rho = tuple(a * b for a, b in zip(eps, xi))
            out[rho] = out.get(rho, Fr(0)) + wt
    return out


def witness(d, t, wR):
    """p_A r(aligned-frozen; [m]) with f(y) = w_R^{y/2} + w_R^{-y/2} inside R (exact; the aligned frozen value of Remark 4.17)."""
    n = 2 * d - 1; t = Fr(t)
    w = (1 + t) / (1 - t); a = w * w / (1 + w * w); pA = 1 - 1 / w ** 2
    f = lambda y: two_cosh(wR, y)
    X = sum(comb(n, k) * a ** k * (1 - a) ** (n - k) / f(2 * k - n + 1) ** 2 for k in range(n + 1))
    Y = sum(comb(n, k) * a ** k * (1 - a) ** (n - k) / f(2 * k - n - 1) ** 2 for k in range(n + 1))
    return pA * a * X / (a * X + (1 - a) * Y)


def up(x, k=9):
    """decimal string of the least k-digit decimal >= x (upper bounds are displayed rounded up)."""
    x = Fr(x); s = 10 ** k; q = -((-x.numerator * s) // x.denominator)
    return f"{q // s}.{q % s:0{k}d}"


def down(x, k=9):
    x = Fr(x); s = 10 ** k; q = (x.numerator * s) // x.denominator
    return f"{q // s}.{q % s:0{k}d}"


def main():
    T0 = time.time(); rng = random.Random(20260926)
    n_id = 0; ok_all = True
    print("(1) nesting identity  R^{beta2}_{nu''}(s) == c R^{beta1}_{nu}(s)  (exact rationals, all s)")
    for m in (2, 4, 6):
        for t2 in (Fr(1, 5), Fr(3, 20), Fr(2, 5)):
            w2 = (1 + t2) / (1 - t2)
            for bpp in (Fr(1, 3), Fr(1, 2), Fr(9, 10), Fr(1)):
                t1 = bpp * t2; w1 = (1 + t1) / (1 - t1)
                c = ((1 - t1 ** 2) / (1 - t2 ** 2)) ** (m // 2)
                for rep in range(3):
                    nu = random_flip_invariant_law(m, rng)
                    nupp = noise(nu, m, bpp)
                    is_law = sum(nupp.values()) == 1 and all(v >= 0 for v in nupp.values())
                    flip_inv = all(nupp.get(tuple(-x for x in r), Fr(0)) == v for r, v in nupp.items())
                    L = R(w2, nupp, m); Rt = R(w1, nu, m)
                    ok = is_law and flip_inv and all(L[s] == c * Rt[s] for s in L)
                    ok_all &= ok; n_id += 1
                print(f"   m={m} t2={t2} b''={bpp}: 3 laws, identity/law/flip-invariance: {'ok' if ok else 'FAIL'}")
    print(f"   {n_id} instances, all exact: {ok_all}")

    print("(2) negative control: b'' = 5/4 > 1 gives P(xi=+1) = 9/8 > 1, not a law ->", end=" ")
    try:
        noise({(1, 1): Fr(1, 2), (-1, -1): Fr(1, 2)}, 2, Fr(5, 4)); print("UNEXPECTED")
    except AssertionError:
        print("rejected (as it must be)")

    print("(3) strictness at the five certified points (d = 9): plain witness < p <= noisy witness (beta_c and beta')")
    J = json.load(open(os.path.join(ROOT, 'data', 'CERT_DATA_d9.json')))
    ok3 = True
    for key in sorted(J, key=lambda k: Fr(J[k]['t'])):
        E = J[key]; d = E['d']; t = Fr(E['t']); p = Fr(E['p']); wcp = Fr(E['wc_prime'])
        n = 2 * d - 1; w = (1 + t) / (1 - t)
        P, Q = (1 + t) ** n, (1 - t) ** n; b = (P - Q) / (P + Q); wc = (1 + t * b) / (1 - t * b)
        plain = witness(d, t, w); noisy_c = witness(d, t, wc); noisy_p = witness(d, t, wcp)
        ok = plain < p <= noisy_c and p <= noisy_p and wc <= wcp <= w
        ok3 &= ok
        print(f"   t={t}: plain {up(plain)} < p = {p} = {float(p):.4f} <= noisy(beta_c) {up(noisy_c)}, "
              f"noisy(beta') {up(noisy_p)}; p/noisy(beta_c) >= {down(p / noisy_c, 6)}; "
              f"w_c <= w_c' <= w: {wc <= wcp <= w} -> {'ok' if ok else 'FAIL'}")
    print(f"ALL OK: {ok_all and ok3}   [{time.time() - T0:.1f}s]")
    sys.exit(0 if ok_all and ok3 else 1)


if __name__ == "__main__":
    main()

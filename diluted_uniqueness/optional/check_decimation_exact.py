"""Exact (symbolic) checks of the isolated-diamond decimation identity (optional; needs SymPy).

Research-note source: theorem-A.md, Lemma 2.4. The lemma is proved by hand; this is a finite check.

All checks are identities of Laurent polynomials / rational functions in x = e^beta,
so they hold for every beta > 0 simultaneously (sympy exact arithmetic, no floating point).

Diamond: terminals u, w; internal a, b; paths u-a-w and u-b-w.
Weight W(s_u, s_w) = sum_{s_a, s_b} x^{ s_a (J_ua s_u + J_aw s_w) + s_b (J_ub s_u + J_bw s_w) }.

(A) unit magnitudes: for all 16 sign patterns, W(s_u,s_w) is independent of (s_u,s_w)
    iff the plaquette is frustrated (product of the four signs = -1); in that case
    W = 2 (x^2 + x^-2) = 4 cosh(2 beta).
(B) the single-path decimation formula: for magnitudes m1, m2 in {1,2,3} and signs,
    sum_s x^{s(e1 m1 s_u + e2 m2 s_w)} = C exp(K s_u s_w) with
    tanh K = tanh(beta e1 m1) tanh(beta e2 m2)   (checked as a rational-function identity
    for e^{2K} = W(+,+)/W(+,-) = (1+t1 t2)/(1-t1 t2)).
(C) magnitudes in {1,2}: W independent of (s_u,s_w) as Laurent polynomials iff
    (sign products of the two paths opposite) and {m_ua,m_aw} = {m_ub,m_bw} as multisets.
"""
import itertools
import sys

import sympy as sp

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

x = sp.symbols('x', positive=True)


def W(Jua, Jaw, Jub, Jbw, su, sw):
    tot = 0
    for sa, sb in itertools.product((-1, 1), repeat=2):
        tot += x ** (sa * (Jua * su + Jaw * sw) + sb * (Jub * su + Jbw * sw))
    return sp.expand(tot)


def check_unit():
    n_frus = n_ok = 0
    for signs in itertools.product((-1, 1), repeat=4):
        Jua, Jaw, Jub, Jbw = signs
        vals = [W(Jua, Jaw, Jub, Jbw, su, sw) for su, sw in itertools.product((-1, 1), repeat=2)]
        const = all(sp.simplify(v - vals[0]) == 0 for v in vals)
        frustrated = (Jua * Jaw * Jub * Jbw == -1)
        assert const == frustrated, (signs, vals)
        if frustrated:
            n_frus += 1
            assert sp.simplify(vals[0] - 2 * (x ** 2 + x ** -2)) == 0
        n_ok += 1
    return n_ok, n_frus


def check_single_path():
    cnt = 0
    for m1, m2 in itertools.product((1, 2, 3), repeat=2):
        for e1, e2 in itertools.product((-1, 1), repeat=2):
            def w1(su, sw):
                return sum(x ** (s * (e1 * m1 * su + e2 * m2 * sw)) for s in (-1, 1))
            ratio = sp.cancel(w1(1, 1) / w1(1, -1))           # = e^{2K}
            t1 = sp.cancel((x ** (2 * e1 * m1) - 1) / (x ** (2 * e1 * m1) + 1))   # tanh(beta e1 m1)
            t2 = sp.cancel((x ** (2 * e2 * m2) - 1) / (x ** (2 * e2 * m2) + 1))
            pred = sp.cancel((1 + t1 * t2) / (1 - t1 * t2))
            assert sp.simplify(ratio - pred) == 0, (m1, m2, e1, e2)
            # symmetry: W(-,-)=W(+,+), W(-,+)=W(+,-)
            assert sp.simplify(w1(-1, -1) - w1(1, 1)) == 0 and sp.simplify(w1(-1, 1) - w1(1, -1)) == 0
            cnt += 1
    return cnt


def check_mixed():
    cnt = cancel = 0
    for mags in itertools.product((1, 2), repeat=4):
        for signs in itertools.product((-1, 1), repeat=4):
            J = [s * m for s, m in zip(signs, mags)]
            vals = [W(*J, su, sw) for su, sw in itertools.product((-1, 1), repeat=2)]
            const = all(sp.simplify(v - vals[0]) == 0 for v in vals)
            pred = (signs[0] * signs[1] * signs[2] * signs[3] == -1) and sorted(mags[:2]) == sorted(mags[2:])
            assert const == pred, (mags, signs)
            cnt += 1
            cancel += const
    return cnt, cancel


if __name__ == '__main__':
    n, nf = check_unit()
    assert (n, nf) == (16, 8)                                  # counts displayed in the manuscript
    print(f'(A) unit magnitudes: {n} sign patterns checked; {nf} frustrated, all with W == 2(x^2+x^-2) '
          f'independent of terminals; all {n - nf} unfrustrated patterns terminal-dependent. PASS')
    c = check_single_path()
    assert c == 36
    print(f'(B) single-path formula e^(2K) = (1+t1 t2)/(1-t1 t2): {c} (magnitude,sign) cases. PASS')
    c, k = check_mixed()
    assert (c, k) == (256, 48)
    print(f'(C) magnitudes {{1,2}}: {c} patterns; {k} cancel identically, exactly those with opposite path '
          f'sign products and equal magnitude multisets. PASS')

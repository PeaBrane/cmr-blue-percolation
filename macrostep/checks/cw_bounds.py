"""Macrostep route: exact Collatz-Wielandt enclosure of the Perron root of a stored computed near matrix Mbar.
[Packaged copy of the research script cw_bounds.py, unchanged apart from this note and an
optimisation guard (needs NumPy). Macrostep write-up Sec. 9.4 and Table 9.5c: exact Collatz-Wielandt enclosures of the
Perron root of a computed near matrix, used only for the robustness statements (older local inputs at d = 7).]

For a nonnegative matrix A and a vector v > 0:  min_i (Av)_i / v_i <= rho(A) <= max_i (Av)_i / v_i.
We take v = the float Perron vector of Mbar (made positive) and evaluate both ratios EXACTLY (integer arithmetic on
the dyadic rationals represented by the stored binary64 entries).  This bounds rho of the COMPUTED Mbar only; the
exact near matrix M satisfies M <= Mbar entrywise, so a lower bound >= 1 excludes a certificate built on this Mbar
(condition (i) of Theorem G forces Mbar w <= lam w, hence rho(Mbar) <= lam < 1), not a certificate for the family.

Usage: python cw_bounds.py file.pkl [...]
"""
import sys, pickle

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
import numpy as np
from recheck_certs import down, up


def dy(x):
    n, den = float(x).as_integer_ratio()
    k = den.bit_length() - 1
    assert den == 1 << k
    return n, k


def scaled(vals):
    pairs = [dy(v) for v in vals]
    K = max(k for _, k in pairs)
    return [n << (K - k) for n, k in pairs], K


def main():
    for path in sys.argv[1:]:
        o = pickle.load(open(path, "rb"))
        M = o['M']
        assert (M >= 0).all()
        ev, V = np.linalg.eig(M)
        k = int(np.argmax(np.abs(ev)))
        v = np.real(V[:, k])
        v = v / v[np.argmax(np.abs(v))]
        v = np.maximum(v, 1e-300)
        W, Kw = scaled(list(v))
        lo = hi = None
        for s in range(len(v)):
            R, Kr = scaled(list(M[s]))
            r = Fr(sum(a * b for a, b in zip(R, W)), 1 << (Kr + Kw)) / Fr(float(v[s]))
            lo = r if lo is None else min(lo, r)
            hi = r if hi is None else max(hi, r)
        print(f"{path}: |C'|={len(v)} float rho(Mbar)={abs(ev[k]):.6f}; exact Collatz-Wielandt enclosure "
              f"{down(lo, 7)} <= rho(Mbar) <= {up(hi, 7)}; rho(Mbar) >= 1 proved: {lo >= 1}", flush=True)


if __name__ == "__main__":
    main()

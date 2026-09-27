"""Exact side results of the single-floor account (standard library only).

None of these is needed for Theorem A.2 of the first manuscript or its corollaries.

(a) The nonnegative-residual floor is not minimised at k = 1 uniformly in
    beta: at (Delta, t) = (16, 1/2), (20, 3/10), (20, 1/3), (32, 1/5) and
    (44, 3/20) the minimum over k is at k = Delta and lies strictly below
    p(1). At seven control points the minimum is at k = 1.
(b) The same-cavity monotonicity p_k(a) >= p_1(a) fails for an explicit law
    a on {-1, 1}^4 at e^beta = 3 (t = 4/5), also after mixing in 1/1000 or
    1/100 of the uniform law. It holds after mixing in 1/10, and for the
    undiluted law at e^beta = 11/10.
(c) Remark A.12: at d = 15, t = 27/250 the overlap counts k >= 2 clear F_15: the signed
    residual at k = 2 with lambda = 1/0.0725, c = 6/5 is negative and the
    nonnegative residual gives min_(k >= 3) p(k) > 0.0725 > F_15. The
    nonnegative residual alone gives p(1) < p(2) < F_15.
"""

from fractions import Fraction as Fr
from itertools import product
from pathlib import Path
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from check_lemmas import nonneg_floor_all_j  # noqa: E402
from verify_single_floor import Star, g_head, nonneg_floor, signed_residual  # noqa: E402

# (Delta, t, argmin k, displayed p(1) and p(argmin), both truncated to 10 decimals)
WITNESSES = [(16, Fr(1, 2), 16, "0.0000131724", "0.0000109283"),
             (20, Fr(3, 10), 20, "0.0060312635", "0.0042399151"),
             (20, Fr(1, 3), 20, "0.0019997080", "0.0011159822"),
             (32, Fr(1, 5), 32, "0.0162538454", "0.0126110043"),
             (44, Fr(3, 20), 44, "0.0254714479", "0.0242009832")]
CONTROLS = [(44, Fr(7, 50)), (16, Fr(3, 10)), (16, Fr(3, 5)), (4, Fr(4, 5)),
            (44, Fr(11, 100)), (44, Fr(1, 8)), (44, Fr(13, 100))]


def truncated(x, shown, digits=10):
    """shown is x truncated to `digits` decimals."""
    return Fr(shown) <= x < Fr(shown) + Fr(1, 10 ** digits)


def floor_row(D, t):
    S = Star(t, D)
    evaluate = nonneg_floor_all_j if D <= 32 else nonneg_floor
    return [evaluate(S, k) for k in range(1, D + 1)]


def p_k_all(Delta, q, law):
    """Exact p_k(a) = E[a_k R^-2] / E[w_k R^-2] at e^beta = q (paper App. A.1, h = 0)."""
    e = Fr(q)
    cosh = lambda m: (e ** m + e ** (-m)) / 2
    alpha = 1 - e ** (-4)
    signs = list(product((-1, 1), repeat=Delta))
    inv_R2 = {}
    for eps in signs:
        R = sum(w * cosh(sum(a * b for a, b in zip(eps, y))) for w, y in law)
        inv_R2[eps] = 1 / (R * R)
    out = []
    for k in range(1, Delta + 1):
        A = D = Fr(0)
        for eps in signs:
            XU, XV = sum(eps[:k]), sum(eps[k:])
            A += alpha / 4 * e ** (2 * eps[0] * XU) * inv_R2[eps]
            D += (cosh(2 * XU) + cosh(2 * XV)) / 2 * inv_R2[eps]
        out.append(A / D)
    return out


def main():
    start = time.monotonic()
    for D, t, k_min, shown_1, shown_min in WITNESSES:
        row = floor_row(D, t)
        assert min(row) == row[k_min - 1] < row[0] and row.index(min(row)) == k_min - 1
        assert truncated(row[0], shown_1) and truncated(row[k_min - 1], shown_min)
        print(f"PASS (a) witness Delta={D} t={t}: argmin k={k_min}, p({k_min}) < p(1) "
              f"({time.monotonic() - start:.1f}s)", flush=True)
    for D, t in CONTROLS:
        row = floor_row(D, t)
        assert min(row) == row[0] and row.index(row[0]) == 0
        print(f"PASS (a) control Delta={D} t={t}: argmin k=1 ({time.monotonic() - start:.1f}s)", flush=True)

    base = [(Fr(75, 200), (1, -1, -1, -1)), (Fr(58, 200), (1, 1, -1, 1)),
            (Fr(58, 200), (1, 1, 1, -1)), (Fr(9, 200), (1, 1, 1, 1))]
    uniform = [(Fr(1, 16), y) for y in product((-1, 1), repeat=4)]
    for delta, bound in [(Fr(0), "0.49344"), (Fr(1, 1000), "0.49961"), (Fr(1, 100), "0.55696"), (Fr(1, 10), None)]:
        law = [((1 - delta) * w, y) for w, y in base] + [(delta * w, y) for w, y in uniform if delta]
        p = p_k_all(4, 3, law)
        if bound is None:
            assert p[1] >= p[0]
        else:
            assert p[1] < p[0] and p[1] / p[0] <= Fr(bound)
    p = p_k_all(4, Fr(11, 10), base)
    assert all(pk >= p[0] for pk in p)
    print(f"PASS (b) same-cavity counterexample and controls ({time.monotonic() - start:.1f}s)", flush=True)

    head, up = g_head(15)
    F_lo, F_up = 1 - 1 / head, 1 - 1 / up
    assert Fr("0.0718626355") <= F_lo and F_up <= Fr("0.0718626462")
    S = Star(Fr(27, 250), 30)
    row = [nonneg_floor_all_j(S, k) for k in range(1, 31)]
    floor = Fr(725, 10000)
    M2, orientations = signed_residual(S, 2, 1 / floor, Fr(6, 5))
    assert orientations == 58 and M2 < 0 and M2 <= Fr("-0.0105469944")
    assert floor > F_up and min(row[2:]) > floor and min(row[2:]) >= Fr("0.072807")
    assert row[0] < row[1] < F_lo
    print(f"PASS (c) d=15 t=27/250: k=2 signed and k>=3 nonnegative clear F_15 <= 0.0718626462 "
          f"({time.monotonic() - start:.1f}s)")


if __name__ == "__main__":
    main()

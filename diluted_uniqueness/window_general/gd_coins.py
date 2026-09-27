"""Axis-star coins for the octahedral rule (Remark 5.27 of the diluted-model manuscript): exact enumeration for specific d.

Coin partition: coin(P_ij(w)) = (w, min(i, j)) for the plaquette w + {0, e_i, e_j, e_i + e_j}; the members of coin (w, i)
are P_ij(w), j > i, which pairwise share the edge {w, w + e_i}.  By Proposition 5.20 the window of Theorem 5.26 holds
with (q, N_Q) replaced by
    q_c(d)  = #{(w, i) : some P_ij(w), j > i, meets O^+}          (O^+ = {|v|_inf <= 2, |v|_1 <= 3})
    N_Qc(d) = d * max_i |V(coin(0, i)) - O^+|
and the same K_d.  Exact rationals, floors; brackets as in gd_table.py.
Usage: python gd_coins.py dmax
"""
import sys
import itertools
import math
from fractions import Fraction as Fr
from gd_table import window, floor_sig, q_NQ


def oplus(d):
    out = set()
    for s in range(4):
        for axes in itertools.combinations(range(d), s):
            for vals in itertools.product((-2, -1, 1, 2), repeat=s):
                if sum(abs(x) for x in vals) <= 3:
                    v = [0] * d
                    for a, x in zip(axes, vals):
                        v[a] = x
                    out.add(tuple(v))
    return out


def counts(d):
    Op = oplus(d)
    coins, plaq = set(), set()
    for v in Op:
        for i, j in itertools.combinations(range(d), 2):
            for a, b in ((0, 0), (1, 0), (0, 1), (1, 1)):
                w = list(v)
                w[i] -= a
                w[j] -= b
                w = bytes(x + 4 for x in w)          # compact key (coordinates lie in [-3, 2])
                plaq.add(w + bytes((i, j)))
                coins.add(w + bytes((i,)))
    assert len(plaq) == q_NQ(d)[0]
    o = (0,) * d
    NQc = 0
    for i in range(d - 1):
        V = {o, tuple(1 if k == i else 0 for k in range(d))}
        for j in range(i + 1, d):
            V.add(tuple(1 if k == j else 0 for k in range(d)))
            V.add(tuple(1 if k in (i, j) else 0 for k in range(d)))
        NQc = max(NQc, d * len({tuple(x - y for x, y in zip(u, w)) for u in V for w in Op}))
    return len(plaq), len(coins), NQc


def main():
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print('d | q | q_c | N_Q | N_Qc | g_d (plaquette coins) >= | g_d (axis-star coins) >= | gain')
    for d in range(3, dmax + 1):
        q, qc, NQc = counts(d)
        q0, NQ, p_lo, p_hi, K, D = window(d)
        tail = Fr(1, 2 ** (qc + 1)) if qc < 400 else Fr(1, 2 ** 400)
        Dc = (1 - tail) / ((qc + 1) * (K + Fr(NQc) / (2 * (1 - p_hi))))
        print(f'{d} | {q} | {qc} | {NQ} | {NQc} | {floor_sig(2 * D)} | {floor_sig(2 * Dc)} | {float(Dc / D):.3f}', flush=True)


if __name__ == '__main__':
    main()

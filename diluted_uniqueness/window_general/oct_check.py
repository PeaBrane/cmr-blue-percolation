"""Clean-room check of the octahedral rule (general d >= 3) in a given dimension d.

Octahedral region O = {v : |v|_inf <= 1, |v|_1 <= 2}, Lambda(e) = z(e) + O.  The rule (types B: G/F1/F2/F3,
N: Nint1/Nint2 and B-with-x'=0, S: S1a/S1b/S2a/S2b/S2c) is implemented from the proof text only; the design
hypotheses are checked with cr_check.check_design (own code).  Checks:
  (1) every configuration of types B, N (all) and S (all at the given L) gets designs satisfying the Design Lemma;
  (2) (|O|, |W*\\O|) equals the value stated for its case;
  (3) the arc-reconstruction consequence: two keys of the rule at the same edge e with the same labelled plaquette
      are incompatible unless both are N-int with arcs at different terminals;
  (4) per labelled plaquette, the realizable-key bound sum_e max_{D at e} W_D c_D (x2 at N-int edges) is at most the
      hand constant K_d = 2d [max(c(7,10d-12)/(2(d-2)), c(6,8d-9)) + c(6,8d-9)];
  (5) for type B in the bulk frame, the exact maximal-clique value (for comparison with K_d).
The optional 5th argument (added after the first audit) selects the p-bracket used in (4) and (5):
  old  -- the first version's bracket (d = 3: [211/1000, 3473/10000]; d >= 4: [1/(2d), 7/20]); default;
  R8   -- the exact Table R8 bracket (p_-, p_+) of gd_table.window(d);
  univ -- [1/(2d), 7/20].
The bracket is printed as exact fractions and log2 values with 6 decimals.
Usage: python oct_check.py d L_S [bulk_clique(0/1)] [store_S_keys(0/1)] [old|R8|univ]
"""
import sys
import itertools
import math
from fractions import Fraction as Fr
from collections import defaultdict
from cr_check import (add, sub, unit, nbrs, edge, linf, l1, Fset, check_design, compatible, max_clique_value, Bad)


def units(d):
    return [unit(d, i, s) for i in range(d) for s in (1, -1)]


def axis(u):
    return next(i for i, x in enumerate(u) if x)


def supp(v):
    d = len(v)
    return [unit(d, i, 1 if x > 0 else -1) for i, x in enumerate(v) if x]


def neg(u):
    return tuple(-x for x in u)


def Pz(z, a, b):
    return [z, add(z, a), add(z, add(a, b)), add(z, b)]


def mk(P, ip, tp, arcs, tag, w=Fr(1)):
    return dict(P=P, intp=list(ip), terp=list(tp), arcs=[list(A) for A in arcs], tag=tag, w=w)


def region(d):
    return frozenset(v for v in itertools.product((-1, 0, 1), repeat=d) if l1(v) <= 2)


def ruleB(d, z, x, y):
    """designs (with weights) for entry x = x', exit y = y' (absolute), region z + O"""
    X, Y = sub(x, z), sub(y, z)
    ab = lambda v: add(z, v)
    for a in supp(X):                                                   # G
        for b in supp(Y):
            if axis(a) != axis(b) and add(a, b) not in (X, Y):
                A1 = [x] if X == a else [x, ab(a)]
                A2 = [y] if Y == b else [y, ab(b)]
                return [mk(Pz(z, a, b), (z, ab(add(a, b))), (ab(a), ab(b)), (A1, A2), 'G')]
    if len(supp(X)) == 1 and len(supp(Y)) == 1:                        # F1: y' = z - a, x' = z + a
        a = X
        assert Y == neg(a)
        c = next(u for u in units(d) if axis(u) != axis(a))
        return [mk(Pz(z, a, c), (z, ab(add(a, c))), (ab(a), ab(c)), ([x], [y, ab(add(neg(a), c)), ab(c)]), 'F1')]
    if len(supp(X)) != len(supp(Y)):                                    # F2: level 1 = z + a, level 2 = z + a + b
        swap = len(supp(X)) == 2
        (p1, P1), (p2, P2) = ((y, Y), (x, X)) if swap else ((x, X), (y, Y))
        a = P1
        assert a in supp(P2)
        b = next(u for u in supp(P2) if u != a)
        cs = [u for u in units(d) if axis(u) not in (axis(a), axis(b))]
        out = []
        for c in cs:
            arc1 = [p1, ab(add(a, c)), ab(c)]
            arc2 = [p2, ab(b)]
            arcs = (arc2, arc1) if swap else (arc1, arc2)
            out.append(mk(Pz(z, b, c), (z, ab(add(b, c))), (ab(b), ab(c)), arcs, 'F2', Fr(1, len(cs))))
        return out
    common = [u for u in supp(X) if u in supp(Y)]                       # F3
    assert len(common) == 1 and {axis(u) for u in supp(X)} == {axis(u) for u in supp(Y)}
    a1 = common[0]
    a2 = next(u for u in supp(X) if u != a1)
    assert Y == add(a1, neg(a2))
    return [mk(Pz(z, a1, a2), (ab(a1), ab(a2)), (z, x), ([x], [y, ab(neg(a2)), z]), 'F3')]


def ruleN(d, z, k, y):
    o = tuple([0] * d)
    if z != o:
        return ruleB(d, z, o, y)
    S = supp(y)
    a = S[0]
    if len(S) == 1:
        b = next(u for u in units(d) if axis(u) != axis(a))
        return [mk(Pz(o, a, b), (o, add(a, b)), (a, b), ([y],), 'Nint1')]
    b = neg(S[1])
    return [mk(Pz(o, a, b), (o, add(a, b)), (a, b), ([y, a],), 'Nint2')]


def ruleS(d, L, z, k, x):
    X = sub(x, z)
    ab = lambda v: add(z, v)
    if linf(z) == L:                                                    # S1: z_k = -L, only k saturated
        assert [i for i in range(d) if abs(z[i]) == L] == [k] and z[k] == -L
        ek = unit(d, k)
        if X == ek:
            g = next(u for u in units(d) if axis(u) != k)
            return [mk(Pz(z, ek, g), (z, ab(add(ek, g))), (x, ab(g)), ([x],), 'S1a')]
        assert ek in supp(X) and len(supp(X)) == 2
        h = next(u for u in supp(X) if u != ek)
        return [mk(Pz(z, ek, h), (ab(ek), ab(h)), (x, z), ([x],), 'S1b')]
    istar = min(i for i in range(d) if abs(z[i]) == L - 1)             # S2
    f = unit(d, istar, 1 if z[istar] > 0 else -1)
    S = supp(X)
    assert f not in S
    if len(S) == 1 and axis(S[0]) != istar:
        h = S[0]
        return [mk(Pz(z, h, f), (z, ab(add(h, f))), (x, ab(f)), ([x],), 'S2a')]
    if len(S) == 1:
        g = next(u for u in units(d) if axis(u) != istar)
        return [mk(Pz(z, f, g), (ab(f), ab(g)), (z, ab(add(f, g))), ([x, z],), 'S2b')]
    h = next(u for u in S if u != neg(f)) if neg(f) in S else S[0]
    return [mk(Pz(z, h, f), (z, ab(add(h, f))), (ab(h), ab(f)), ([x, ab(h)],), 'S2c')]


def expected(tag, d, arcs):
    if tag == 'G':
        return [(4, 4 * d - 4), (5, 6 * d - 7), (6, 8 * d - 10)][sum(1 for A in arcs if len(A) >= 2)]
    return {'F1': (6, 8 * d - 9), 'F2': (7, 10 * d - 12), 'F3': (6, 8 * d - 9), 'Nint1': (4, 4 * d - 4),
            'Nint2': (5, 6 * d - 7), 'S1a': (4, 4 * d - 4), 'S1b': (4, 4 * d - 4), 'S2a': (4, 4 * d - 4),
            'S2b': (5, 6 * d - 7), 'S2c': (5, 6 * d - 7)}[tag]


def entries(Lam):
    return sorted(v for v in Lam if any(u not in Lam for u in nbrs(v)))


def main():
    d, L = int(sys.argv[1]), int(sys.argv[2])
    bulk_clique = len(sys.argv) > 3 and sys.argv[3] == '1'
    storeS = not (len(sys.argv) > 4 and sys.argv[4] == '0')     # 0: check S designs without storing their keys
    o = tuple([0] * d)
    O0 = region(d)
    ent0 = entries(O0)
    assert ent0 == sorted(v for v in O0 if v != o)
    # p-interval used for the check of (4): the d=3 values of the theorem (any p <= 1/2 gives the same ordering)
    bracket = sys.argv[5] if len(sys.argv) > 5 else 'old'
    if bracket == 'old':
        p_lo, p_hi = (Fr(211, 1000), Fr(3473, 10000)) if d == 3 else (Fr(1, 2 * d), Fr(35, 100))
    elif bracket == 'R8':
        from gd_table import window
        p_lo, p_hi = window(d)[2:4]
    elif bracket == 'univ':
        p_lo, p_hi = Fr(1, 2 * d), Fr(7, 20)
    else:
        raise SystemExit(f'unknown bracket {bracket}')
    ps = (p_lo, p_hi)
    c = lambda a, b: max((1 / p) ** a * (1 / (1 - p)) ** b for p in ps)
    Kd = 2 * d * (max(c(7, 10 * d - 12) / (2 * (d - 2)), c(6, 8 * d - 9)) + c(6, 8 * d - 9))
    stats = defaultdict(set)
    ncfg = defaultdict(int)

    def mknode(typ, e, Lam, D, r):
        nO, nZ = len(r['O']), len(r['W']) - len(r['O'])
        assert (nO, nZ) == expected(D['tag'], d, r['arcs']), (D['tag'], nO, nZ)
        stats[(typ, D['tag'])].add((nO, nZ))
        pat = {f: (1 if f in r['O'] else 0) for f in r['W']}
        cs = tuple((1 / p) ** nO * (1 / (1 - p)) ** nZ for p in ps)
        return dict(pat=pat, variants={tuple(sorted(r['anchors']))}, Lam=frozenset(Lam), e=e, c=cs, W=Fr(0),
                    sub=typ, tag=D['tag'], tarc=frozenset(A[-1] for A in r['arcs']), grp=(r['P'], r['ip']))

    def add_config(store, lst):
        """lst: [(node, weight)] for one configuration; merge into store[grp][e][patkey]"""
        tot = defaultdict(Fr)
        for nd, w in lst:
            kk = (nd['grp'], nd['e'], frozenset(nd['pat'].items()))
            cur = store[nd['grp']][nd['e']].get(kk)
            if cur is None:
                store[nd['grp']][nd['e']][kk] = nd
                cur = nd
            else:
                assert cur['c'] == nd['c'] and cur['Lam'] == nd['Lam']
                cur['variants'] |= nd['variants']
            tot[kk] += w
        for kk, w in tot.items():
            nd = store[kk[0]][kk[1]][kk]
            nd['W'] = max(nd['W'], w)

    # ---- type B in the bulk frame (z = 0, source far away)
    KB = defaultdict(lambda: defaultdict(dict))
    for k in range(d):
        e = edge(o, unit(d, k))
        for x in ent0:
            for y in ent0:
                if x != y:
                    lst = []
                    for D in ruleB(d, o, x, y):
                        r = check_design('B', e, O0, D, xp=x, yp=y, src=None)
                        lst.append((mknode('B', e, O0, D, r), D['w']))
                    add_config(KB, lst)
                    ncfg['B'] += 1
    # ---- type N (absolute)
    KN = defaultdict(lambda: defaultdict(dict))
    for z in sorted({sub(o, v) for v in O0}):
        Lam = frozenset(add(z, v) for v in O0)
        for k in range(d):
            e = edge(z, add(z, unit(d, k)))
            for y in entries(Lam):
                if y != o:
                    lst = []
                    for D in ruleN(d, z, k, y):
                        typ = 'Nint' if o in D['intp'] else 'Nter'
                        r = check_design(typ, e, Lam, D, yp=y, src=o)
                        lst.append((mknode(typ, e, Lam, D, r), D['w']))
                    add_config(KN, lst)
                    ncfg['N'] += 1
    # ---- type S (absolute, radius L)
    KS = defaultdict(lambda: defaultdict(dict))
    for z in itertools.product(range(-L, L + 1), repeat=d):
        if linf(z) < L - 1:
            continue
        Lam = frozenset(add(z, v) for v in O0)
        if o in Lam or all(linf(v) <= L - 1 for v in Lam):
            continue
        for k in range(d):
            z2 = add(z, unit(d, k))
            if linf(z2) > L or (linf(z) > L - 1 and linf(z2) > L - 1):
                continue
            e = edge(z, z2)
            for x in sorted(Lam):
                if linf(x) > L - 1 or not any(u not in Lam and linf(u) <= L - 1 for u in nbrs(x)):
                    continue
                lst = []
                for D in ruleS(d, L, z, k, x):
                    r = check_design('S', e, Lam, D, L=L, xp=x, src=o)
                    assert [t for t in r['P'] - r['ip'] if linf(t) <= L - 1] == [D['arcs'][0][-1]]
                    lst.append((mknode('S', e, Lam, D, r), D['w']))
                if storeS:
                    add_config(KS, lst)
                ncfg['S'] += 1
    print(f'd={d}: every configuration has rule designs satisfying the Design Lemma: '
          f'B {ncfg["B"]}, N {ncfg["N"]}, S(L={L}) {ncfg["S"]}')
    for (typ, tag), v in sorted(stats.items()):
        print(f'  {typ:4s} {tag:5s}: (|O|, |W*\\O|) = {sorted(v)}')

    def tr_node(nd, t):
        trv = lambda v: add(v, t)
        return dict(nd, pat={edge(trv(f[0]), trv(f[1])): b for f, b in nd['pat'].items()},
                    variants={tuple(trv(u) for u in var) for var in nd['variants']},
                    Lam=frozenset(trv(v) for v in nd['Lam']), e=edge(trv(nd['e'][0]), trv(nd['e'][1])),
                    tarc=frozenset(trv(v) for v in nd['tarc']),
                    grp=(frozenset(trv(v) for v in nd['grp'][0]), frozenset(trv(v) for v in nd['grp'][1])))

    def keys_at(P, ip, z, k, where):
        """rule keys with labelled plaquette (P, ip) at the edge (z, z+e_k); where = 'bulk' / 'abs' (L fixed)"""
        e = edge(z, add(z, unit(d, k)))
        Lam = frozenset(add(z, v) for v in O0)
        if where == 'abs' and o in Lam:
            return list(KN.get((P, ip), {}).get(e, {}).values())
        if where == 'abs' and not all(linf(v) <= L - 1 for v in Lam):
            return list(KS.get((P, ip), {}).get(e, {}).values())
        g0 = (frozenset(sub(v, z) for v in P), frozenset(sub(v, z) for v in ip))
        e0 = edge(o, unit(d, k))
        return [tr_node(nd, z) for nd in KB.get(g0, {}).get(e0, {}).values()]

    # ---- (3), (4) per labelled plaquette: bulk frame (all B), and every labelled plaquette of an N or S key
    groups = []
    for (P0, ip0) in KB:
        m = min(P0)
        groups.append(('bulk', frozenset(sub(v, m) for v in P0), frozenset(sub(v, m) for v in ip0)))
    for (P0, ip0) in list(KN) + list(KS):
        groups.append(('abs', P0, ip0))
    groups = sorted(set(groups), key=repr)
    worst_ratio, npairs, nmix = Fr(0), 0, 0
    best_clique = Fr(0)
    for where, P0, ip0 in groups:
        tot = Fr(0)
        allnodes = []
        types = set()
        for z in P0:
            for k in range(d):
                lst = keys_at(P0, ip0, z, k, where)
                if not lst:
                    continue
                types |= {nd['sub'] for nd in lst}
                allnodes += lst
                for i, j in itertools.combinations(range(len(lst)), 2):
                    A, B = lst[i], lst[j]
                    npairs += 1
                    if compatible(A, B):
                        assert A['sub'] == B['sub'] == 'Nint' and A['tarc'] != B['tarc'], (A['tag'], B['tag'])
                mult = 2 if any(nd['sub'] == 'Nint' for nd in lst) else 1
                tot += mult * max(nd['W'] * max(nd['c']) for nd in lst)
        if len(types) > 1:
            nmix += 1
        worst_ratio = max(worst_ratio, tot / Kd)
        if bulk_clique and where == 'bulk':
            v, _, _, _ = max_clique_value(allnodes, ps)
            best_clique = max(best_clique, v)
    print(f'  {len(groups)} labelled plaquettes checked ({nmix} carry keys of two types); {npairs} same-edge key pairs: '
          f'all incompatible except N-int pairs with arcs at different terminals')
    print(f'  max over labelled P0 of sum_e max_(D at e) W c  /  K_d = {float(worst_ratio):.6f} (<= 1 required)')
    assert worst_ratio <= 1
    print(f'  hand constant: log2 K_d = {math.log2(Kd):.4f} at p in [{float(p_lo):.6f}, {float(p_hi):.6f}]')
    print(f'  [rev1] bracket "{bracket}": p_lo = {p_lo}, p_hi = {p_hi}; log2 K_d = {math.log2(Kd):.6f}')
    if bulk_clique:
        print(f'  exact bulk max-clique value: 2^{math.log2(best_clique):.4f}  (hand bound 2^{math.log2(Kd):.4f})')
        print(f'  [rev1] bracket "{bracket}": exact bulk max-clique value 2^{math.log2(best_clique):.6f} '
              f'<= K_d = 2^{math.log2(Kd):.6f}: {best_clique <= Kd}')
        assert best_clique <= Kd


if __name__ == '__main__':
    main()

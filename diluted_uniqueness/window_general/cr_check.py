"""Clean-room verification of the sparse-insulation surgery certificates (d = 3 unit cube, d = 2 box).

Independent of the programs that produced the rule tables: nothing from them is imported.  The only input is a
rule table (params/lp_<d>_<region>_rho*.json: randomized rules for types B and N).  Everything else is re-derived
from the definitions (regions, covering sets and types; the design lemma; key costs; compatibility; coins):

  * Design-Lemma hypotheses ((H1)-(H6) of Definition 5.5 and Lemma 5.6 of the diluted-model manuscript) for
    types B, N-ter, N-int, S;
  * exact per-key costs  c_kappa(p) = rho_hi^n(kappa) p^-|O| (1-p)^-|W*\\O|   (Definition 5.8, Lemma 5.9);
  * node weights W_kappa = max over configurations c of the total option weight of c with key kappa;
  * the compatibility relation (Definition 5.11, Corollary 5.12) and maximal-clique enumeration (own Bron-Kerbosch with pivoting);
  * K_B (bulk groups), K_BN (joint groups near the origin: type-N keys together with the translated type-B keys
    whose region avoids 0), K_S (an own picture-covariant type-S rule, verified at several L; bound = sum over
    ALL keys with the same labelled plaquette, which needs no compatibility argument);
  * q, N_Q (plaquette coins) and q_c, N_Qc for vertex coins (d = 3) and domino coins (d = 2), Section 5.5;
  * Delta = (1 - 2^-(q+1)) / ((q+1) (K + N_Q / (2 (1 - p_hi)))) as an exact rational, floored to 7 digits.

The output also reports K_B, K_N at each endpoint of the p-bracket (K is the maximum over both endpoints, Lemma
5.15), the type-B configurations with {x', y'} = e, the maximal single-design cost, and the composition of the
maximizing type-B clique.  main() returns the exact values for check_window_general.py.

Usage: python cr_check.py <rule.json> L_S [L_S ...]
"""
import sys
import json
import math
import itertools
from fractions import Fraction as Fr
from collections import defaultdict


# ------------------------------------------------------------------ lattice helpers
def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def unit(d, i, s=1):
    v = [0] * d
    v[i] = s
    return tuple(v)


def nbrs(v):
    for i in range(len(v)):
        for s in (1, -1):
            yield add(v, unit(len(v), i, s))


def edge(u, v):
    return (u, v) if u < v else (v, u)


def linf(v):
    return max(abs(x) for x in v)


def l1(v):
    return sum(abs(x) for x in v)


def region0(d, kind):
    if kind == 'cube':
        return frozenset(itertools.product((0, 1), repeat=d))
    if kind == 'box':
        return frozenset(itertools.product((-1, 0, 1), repeat=d))
    raise ValueError(kind)


def Fset(S):
    """F(S): edges with at least one endpoint in S"""
    out = set()
    for v in S:
        for u in nbrs(v):
            out.add(edge(u, v))
    return out


def plaq_frame(P):
    """(m, f, g) with P = {m, m+f, m+g, m+f+g}, f != g positive unit vectors; None if P is not a plaquette"""
    P = sorted(set(P))
    if len(P) != 4:
        return None
    m = P[0]
    D = [sub(v, m) for v in P[1:]]
    us = [w for w in D if l1(w) == 1 and min(w) >= 0]
    if len(us) != 2 or add(us[0], us[1]) not in D:
        return None
    return m, us[0], us[1]


def plaq_edges(P):
    m, f, g = plaq_frame(P)
    a, b = add(m, f), add(m, g)
    c = add(a, g)
    return {edge(m, a), edge(m, b), edge(a, c), edge(b, c)}


def diag_pairs(P):
    m, f, g = plaq_frame(P)
    return (frozenset({m, add(add(m, f), g)}), frozenset({add(m, f), add(m, g)}))


def plaquettes_in(S):
    S = set(S)
    d = len(next(iter(S)))
    out = []
    for m in sorted(S):
        for i, j in itertools.combinations(range(d), 2):
            a, b = add(m, unit(d, i)), add(m, unit(d, j))
            c = add(a, unit(d, j))
            if a in S and b in S and c in S:
                out.append(frozenset((m, a, b, c)))
    return out


def plaquettes_with_corner(v):
    d = len(v)
    out = set()
    for i, j in itertools.combinations(range(d), 2):
        for si in (1, -1):
            for sj in (1, -1):
                a, b = add(v, unit(d, i, si)), add(v, unit(d, j, sj))
                out.add(frozenset((v, a, b, add(a, unit(d, j, sj)))))
    return out


def lower_corner(P):
    return min(P)


# ------------------------------------------------------------------ Design-Lemma hypotheses
class Bad(Exception):
    pass


def check_design(typ, e, Lam, D, L=None, xp=None, yp=None, src=None):
    """Hypotheses (H1)-(H6) of Definition 5.5 (the design lemma is Lemma 5.6).  typ in {B, Nter, Nint, S}; L only for type S.
    src = position of the source vertex 0 of Pi in the current frame (None: type B in the bulk frame, where the
    source is outside Lambda by the type condition).  Returns dict with I, W* (as W), O, anchored attachments, n."""
    d = len(e[0])
    o = src
    P = frozenset(tuple(v) for v in D['P'])
    if plaq_frame(P) is None:
        raise Bad('not a plaquette')
    ip = frozenset(tuple(v) for v in D['intp'])
    tp = frozenset(tuple(v) for v in D['terp'])
    if ip not in diag_pairs(P) or tp not in diag_pairs(P) or ip == tp:
        raise Bad('labelling')
    if not P <= Lam:                                                     # (H1)
        raise Bad('P0 not in Lambda')
    arcs = [tuple(tuple(v) for v in A) for A in D['arcs']]
    used, ends = set(), set()
    for A in arcs:                                                       # (H2)
        if not A:
            raise Bad('empty arc')
        if any(l1(sub(u, v)) != 1 for u, v in zip(A, A[1:])):
            raise Bad('arc not a lattice path')
        if len(set(A)) != len(A):
            raise Bad('arc not self-avoiding')
        if not set(A) <= Lam:
            raise Bad('arc leaves Lambda')
        if A[-1] not in tp:
            raise Bad('arc does not end at a terminal')
        if set(A[:-1]) & P:
            raise Bad('arc meets V(P0) before its last vertex')
        if A[-1] in ends:
            raise Bad('two arcs end at the same terminal')
        if set(A) & used:
            raise Bad('arcs not vertex-disjoint')
        ends.add(A[-1])
        used |= set(A)
    att = [A[0] for A in arcs]
    I = set(ip)
    for A in arcs:
        I |= set(A[1:])
    if set(att) & I:                                                     # (H3)
        raise Bad('attachment insulated')
    if not I <= Lam:
        raise Bad('I not in Lambda')
    if typ == 'B':                                                       # (H4) type structure
        if att != [xp, yp] or (o is not None and o in Lam):
            raise Bad('B structure')
    elif typ == 'Nter':
        if att != [o, yp] or o in I or o not in Lam:
            raise Bad('N-ter structure')
    elif typ == 'Nint':
        if att != [yp] or o not in ip:
            raise Bad('N-int structure')
    elif typ == 'S':
        if att != [xp] or o in I:
            raise Bad('S structure')
        if any(linf(v) > L - 1 for v in arcs[0]):
            raise Bad('S: arc leaves B_{L-1}')
        if not any(linf(v) == L for v in ip):
            raise Bad('S: no internal vertex on S_L')
    else:
        raise Bad('unknown type')
    if typ in ('Nter', 'Nint') and o not in Lam:
        raise Bad('N: 0 not in Lambda')
    if not any(len(A) >= 2 or A[0] != o for A in arcs):                  # (H5) a terminal of degree >= 3
        raise Bad('no terminal of degree >= 3')
    O = set(plaq_edges(P))
    for A in arcs:
        O |= {edge(u, v) for u, v in zip(A, A[1:])}
    W = Fset(I)
    if typ != 'Nint':                                                    # (H6) W* = F(I) + e
        W.add(e)
    if not O <= W:
        raise Bad('O not inside W*')
    n = sum(1 for w in e if w in I and w != o and (L is None or linf(w) < L))
    return dict(P=P, ip=ip, I=frozenset(I), W=frozenset(W), O=frozenset(O),
                anchors=tuple(u for u in att if u != o), n=n, arcs=arcs)


def cost(nO, nZ, n, p, rho):
    return rho ** n * (1 / p) ** nO * (1 / (1 - p)) ** nZ


# ------------------------------------------------------------------ compatibility and cliques
def compatible(X, Y):
    """necessary condition for two keys to be realized by one output (Definition 5.11, Corollary 5.12).  A key may be produced by
    several designs that differ only in which vertices must carry a Pi-anchor (N-int with a trivial arc at
    either terminal); X['variants'] lists these anchor sets and the key is realizable if one of them is met."""
    px, py = X['pat'], Y['pat']
    if len(px) > len(py):
        px, py = py, px
    for f, v in px.items():
        w = py.get(f)
        if w is not None and w != v:
            return False
    for A, B in ((X, Y), (Y, X)):
        if not any(all(any(w not in A['Lam'] and B['pat'].get(edge(u, w)) != 0 and A['pat'].get(edge(u, w)) != 0
                           for w in nbrs(u)) for u in var) for var in A['variants']):
            return False
    return True


def maximal_cliques(n, adj):
    """Bron-Kerbosch with Tomita pivoting; adj: list of sets over range(n)"""
    out = []

    def bk(R, P, X):
        if not P and not X:
            out.append(R)
            return
        piv = max(P | X, key=lambda w: len(P & adj[w]))
        for v in list(P - adj[piv]):
            bk(R | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}
    bk(frozenset(), set(range(n)), set())
    return out


def max_clique_value(nodes, ps, by_endpoint=False):
    """max over maximal cliques and p in ps of sum W_kappa c_kappa(p); also returns #cliques and #edges.
    by_endpoint=True: additionally returns [(max over cliques at ps[i], argmax clique) for each i]."""
    n = len(nodes)
    adj = [set() for _ in range(n)]
    ne = 0
    for i, j in itertools.combinations(range(n), 2):
        if compatible(nodes[i], nodes[j]):
            adj[i].add(j)
            adj[j].add(i)
            ne += 1
    cls = maximal_cliques(n, adj)
    best, arg = Fr(0), None
    ep = [(Fr(0), None) for _ in ps]
    for cl in cls:
        for i, p in enumerate(ps):
            s = sum(nodes[j]['W'] * nodes[j]['c'][i] for j in cl)
            if s > best:
                best, arg = s, (sorted(cl), p)
            if s > ep[i][0]:
                ep[i] = (s, sorted(cl))
    if by_endpoint:
        return best, len(cls), ne, arg, ep
    return best, len(cls), ne, arg


# ------------------------------------------------------------------ type S: own picture-covariant rule
def simple_paths(allowed, s, t):
    if s == t:
        return [(s,)]
    out, path, seen = [], [s], {s}

    def rec(v):
        for u in nbrs(v):
            if u in allowed and u not in seen:
                path.append(u)
                seen.add(u)
                if u == t:
                    out.append(tuple(path))
                else:
                    rec(u)
                path.pop()
                seen.discard(u)
    rec(s)
    return out


def configs_S(d, Lam0, L):
    """(z, k, x'): e = {z, z+e_k} inside B_L with an endpoint in B_{L-1}; Lambda = z + Lam0 not inside B_{L-1};
    0 not in Lambda; x' in Lambda cap B_{L-1} with a neighbour in B_{L-1} minus Lambda"""
    o = tuple([0] * d)
    out = []
    for z in itertools.product(range(-L, L + 1), repeat=d):
        Lam = {add(z, v) for v in Lam0}
        if o in Lam or all(linf(v) <= L - 1 for v in Lam):
            continue
        for k in range(d):
            z2 = add(z, unit(d, k))
            if linf(z2) > L or (linf(z) > L - 1 and linf(z2) > L - 1):
                continue
            for xp in sorted(Lam):
                if linf(xp) <= L - 1 and any(u not in Lam and linf(u) <= L - 1 for u in nbrs(xp)):
                    out.append((z, k, xp))
    return out


def ruleS(d, Lam0, L, z, k, xp, p_lo, p_hi, rho):
    """min over candidates of max(c(p_lo), c(p_hi)); ties by the design written relative to z"""
    Lam = {add(z, v) for v in Lam0}
    e = edge(z, add(z, unit(d, k)))
    best = None
    for P in plaquettes_in(Lam):
        for ip, tp in (diag_pairs(P), diag_pairs(P)[::-1]):
            if xp in ip or not any(linf(v) == L for v in ip):
                continue
            us = [t for t in tp if linf(t) <= L - 1]
            if len(us) != 1:
                continue
            u = us[0]
            w = next(t for t in tp if t != u)
            allowed = {v for v in Lam if linf(v) <= L - 1} - set(ip) - {w}
            if xp not in allowed or u not in allowed:
                continue
            for A in simple_paths(allowed, xp, u):
                D = dict(P=sorted(P), intp=sorted(ip), terp=[u, w], arcs=[list(A)])
                r = check_design('S', e, Lam, D, L=L, xp=xp, src=tuple([0] * d))
                nO, nZ = len(r['O']), len(r['W']) - len(r['O'])
                val = max(cost(nO, nZ, r['n'], p_lo, rho), cost(nO, nZ, r['n'], p_hi, rho))
                rel = (tuple(sorted(sub(v, z) for v in P)), tuple(sorted(sub(v, z) for v in ip)),
                       tuple(sub(v, z) for v in A))
                if best is None or (val, rel) < (best[0], best[1]):
                    best = (val, rel, D, r)
    return best


# ------------------------------------------------------------------ main
def floor_sig(x, digits=7):
    """largest decimal with `digits` significant digits that is <= x (x > 0 rational)"""
    e10 = math.floor(math.log10(float(x)))
    scale = Fr(10) ** (digits - 1 - e10)
    m = (x * scale).numerator // (x * scale).denominator
    return f'{m / 10 ** (digits - 1):.{digits - 1}f}e{e10}', Fr(m) / scale


def window(q, NQ, K, p_hi):
    return (1 - Fr(1, 2 ** (q + 1))) / ((q + 1) * (K + Fr(NQ) / (2 * (1 - p_hi))))


def main():
    data = json.load(open(sys.argv[1]))
    Ls = [int(a) for a in sys.argv[2:]]
    d, kind = data['d'], data['region']
    p_lo, p_hi = Fr(data['p_lo']), Fr(data['p_hi'])
    ps = (p_lo, p_hi)
    rho = 1 - (1 - p_hi) ** (2 * d - 1)
    o = tuple([0] * d)
    Lam0 = region0(d, kind)
    ent0 = sorted(v for v in Lam0 if any(u not in Lam0 for u in nbrs(v)))
    print(f'clean-room check: d={d} region={kind} |Lambda_0|={len(Lam0)} entries={len(ent0)}; '
          f'p in [{p_lo}, {p_hi}]; rho_hi = 1-(1-p_hi)^{2*d-1}')

    # ---- q, N_Q, and vertex coins
    Lp0 = set(Lam0) | {u for v in Lam0 for u in nbrs(v)}
    Q0 = set()
    for v in Lp0:
        Q0 |= plaquettes_with_corner(v)
    q = len(Q0)
    planes0 = [frozenset((o, unit(d, i), unit(d, j), add(unit(d, i), unit(d, j))))
               for i, j in itertools.combinations(range(d), 2)]
    NQ = d * max(len({sub(v, w) for v in P for w in Lp0}) for P in planes0)
    print(f'  |Lambda^+|={len(Lp0)}  q = |Q(e)| = {q}   N_Q = max_P #{{e: P in Q(e)}} = {NQ}   (JSON: q={data["q"]}, N_Q={data["NQ"]})')
    coins = None
    if d in (2, 3):
        # coin partitions (Section 5.5): every coin is a set of plaquettes pairwise sharing an edge.
        #   d = 3: vertex coins, coin(P) = lower corner of P (the 3 plaquettes with a common lower corner);
        #   d = 2: domino coins, coin(P) = the vertical edge of P whose lower endpoint has even first coordinate
        #          (so a coin is a horizontal domino {P(x), P(x - e_1)}, x_1 even).
        if d == 3:
            coin_of = lambda P: min(P)
            name = 'vertex coins'
        else:
            def coin_of(P):
                m = min(P)
                x = m if m[0] % 2 == 0 else add(m, unit(2, 0))
                return (x, add(x, unit(2, 1)))
            name = 'domino coins'
        # members of a coin, and a check that members pairwise share an edge
        members = defaultdict(set)
        for x in itertools.product(range(-4, 5), repeat=d):
            for P in plaquettes_with_corner(x):
                members[coin_of(P)].add(P)
        for cn, Ps in members.items():
            if all(abs(t) <= 2 for t in min(min(P) for P in Ps)):
                for P1, P2 in itertools.combinations(Ps, 2):
                    assert len(set(plaq_edges(P1)) & set(plaq_edges(P2))) >= 1, 'coin members must share an edge'
        # q_c = max over the period classes of z(e) of #coins meeting Lambda^+(e); N_Qc = d max_coin #{z}
        qc = 0
        for z in itertools.product((0, 1), repeat=d):
            Qz = set()
            for v in Lp0:
                Qz |= plaquettes_with_corner(add(v, z))
            qc = max(qc, len({coin_of(P) for P in Qz}))
        NQc = 0
        for cn, Ps in members.items():
            if all(abs(t) <= 1 for t in min(min(P) for P in Ps)):
                V = set().union(*Ps)
                NQc = max(NQc, d * len({sub(v, w) for v in V for w in Lp0}))
        if d == 3:
            V3 = set().union(*planes0)
            assert {lower_corner(P) for P in Q0} == {sub(w, v) for w in Lp0 for v in V3}
        coins = (qc, NQc, name)
        print(f'  {name}: q_c = max_e #{{coins meeting Lambda^+(e)}} = {qc}   N_Qc = max_coin #{{e}} = {NQc}')

    # ---- types B and N from the rule table
    nodes = {'B': defaultdict(dict), 'N': defaultdict(dict)}
    seen = {'B': set(), 'N': set()}
    worst = {'B': Fr(0), 'N': Fr(0)}
    ncheck = 0
    opt_rows = {'B': [], 'N': []}          # (cfg, w, |I|, |O|, |W*\O|, n, max cost) per option
    for typ in ('B', 'N'):
        for entry in data['rule' + typ]:
            c = entry['config']
            if typ == 'B':
                k, xp, yp = c[0], tuple(c[1]), tuple(c[2])
                z = o
                Lam = set(Lam0)
                cfg = (k, xp, yp)
            else:
                z, k, yp = tuple(c[0]), c[1], tuple(c[2])
                xp = None
                Lam = {add(z, v) for v in Lam0}
                cfg = (z, k, yp)
            assert cfg not in seen[typ], ('duplicate configuration', cfg)
            seen[typ].add(cfg)
            e = edge(z, add(z, unit(d, k)))
            ws = [Fr(opt['w']) for opt in entry['options']]
            assert all(w > 0 for w in ws) and sum(ws) == 1, ('weights', cfg)
            wkey = defaultdict(Fr)
            for opt, w in zip(entry['options'], ws):
                D = opt['design']
                if typ == 'B':
                    sub_t = 'B'
                    r = check_design('B', e, Lam, D, xp=xp, yp=yp, src=None)
                else:
                    sub_t = 'Nint' if o in {tuple(v) for v in D['intp']} else 'Nter'
                    r = check_design(sub_t, e, Lam, D, yp=yp, src=o)
                ncheck += 1
                assert bool(D.get('closeE', sub_t != 'Nint')) == (sub_t != 'Nint') or e in Fset(r['I'])
                nO, nZ = len(r['O']), len(r['W']) - len(r['O'])
                cs = tuple(cost(nO, nZ, r['n'], p, rho) for p in ps)
                worst[typ] = max(worst[typ], *cs)
                if typ == 'B':
                    t = sub(o, min(r['P']))                          # frame: min corner of P0 at the origin
                else:
                    t = o
                tr = lambda v, t=t: add(v, t)
                grp = (frozenset(tr(v) for v in r['P']), frozenset(tr(v) for v in r['ip']))
                pat = {edge(tr(f[0]), tr(f[1])): (1 if f in r['O'] else 0) for f in r['W']}
                et = edge(tr(e[0]), tr(e[1]))
                key = (et, frozenset(pat.items()))
                wkey[(grp, key)] += w
                nd = nodes[typ][grp].get(key)
                var = tuple(sorted(tr(u) for u in r['anchors']))
                data_nd = dict(pat=pat, variants={var}, Lam=frozenset(tr(v) for v in Lam),
                               e=et, c=cs, sub=sub_t, nO=nO, nZ=nZ, n=r['n'], nI=len(r['I']))
                if nd is None:
                    data_nd['W'] = Fr(0)
                    nodes[typ][grp][key] = data_nd
                else:
                    for fld in ('Lam', 'c', 'sub'):
                        assert nd[fld] == data_nd[fld], ('inconsistent node', fld)
                    nd['variants'].add(var)
                opt_rows[typ].append((cfg, w, len(r['I']), nO, nZ, r['n'], max(cs)))
            for (grp, key), w in wkey.items():
                nd = nodes[typ][grp][key]
                nd['W'] = max(nd['W'], w)
                nd.setdefault('cfgw', {})[cfg] = w
    # completeness of the configuration lists
    allB = {(k, xp, yp) for k in range(d) for xp in ent0 for yp in ent0 if xp != yp}
    allN = set()
    for z in sorted({sub(o, v) for v in Lam0}):
        Lam = {add(z, v) for v in Lam0}
        ent = [v for v in Lam if any(u not in Lam for u in nbrs(v))]
        allN |= {(z, k, yp) for k in range(d) for yp in ent if yp != o}
    assert seen['B'] == allB, ('B configurations incomplete', len(seen['B']), len(allB))
    assert seen['N'] == allN, ('N configurations incomplete', len(seen['N']), len(allN))
    rN = max(linf(add(z, v)) for z in {sub(o, v) for v in Lam0} for v in Lam0)
    print(f'  type-N regions lie in B_{rN}: the N-type hypotheses I cap S_L = {{}} hold for L >= {rN + 1}')
    print(f'  rule table: {len(allB)} B and {len(allN)} N configurations (complete), {ncheck} options; '
          f'all satisfy the Design-Lemma hypotheses; weights positive, summing to 1')
    for typ in ('B', 'N'):
        print(f'  worst single {typ}-key cost 2^{math.log2(worst[typ]):.4f}')

    for typ in ('B', 'N'):
        mv = sum(1 for dct in nodes[typ].values() for nd in dct.values() if len(nd['variants']) > 1)
        print(f'  type {typ}: keys produced by designs with different anchor sets: {mv}')
    Ks = {}
    ep_arg = {}
    ep_vals = {}
    for typ in ('B', 'N'):
        best, ncl, nn, ne = Fr(0), 0, 0, 0
        ep_best = [(Fr(0), None, None) for _ in ps]
        for grp, dct in nodes[typ].items():
            lst = list(dct.values())
            nn += len(lst)
            v, n1, n2, arg, ep = max_clique_value(lst, ps, by_endpoint=True)
            ncl += n1
            ne += n2
            best = max(best, v)
            for i, (s, cl) in enumerate(ep):
                if s > ep_best[i][0]:
                    ep_best[i] = (s, grp, [lst[j] for j in cl])
        Ks[typ] = best
        print(f'  type {typ}: {len(nodes[typ])} groups, {nn} keys, {ne} compatible pairs, {ncl} maximal cliques; '
              f'K_{typ} = 2^{math.log2(best):.4f}')
        # Lemma 5.15: K = max(K(p_lo), K(p_hi)); print both endpoint values to 6 decimals
        assert best == max(s for s, _, _ in ep_best)
        iat = max(range(len(ps)), key=lambda i: ep_best[i][0])
        print(f'    K_{typ}(p_lo={ps[0]}) = 2^{math.log2(ep_best[0][0]):.6f}, '
              f'K_{typ}(p_hi={ps[1]}) = 2^{math.log2(ep_best[1][0]):.6f}; '
              f'K_{typ} = max = 2^{math.log2(best):.6f} = {float(best):.6f}, attained at p = {ps[iat]}')
        ep_arg[typ] = (iat, ep_best[iat])
        ep_vals[typ] = [v for v, _, _ in ep_best]

    # ---- the type-B configurations with {x', y'} = e, the maximal single-design cost, and the argmax clique
    rowsB = opt_rows['B']
    unitv = lambda k: unit(d, k)
    ecfg = sorted({cfg for cfg, *_ in rowsB if {cfg[1], cfg[2]} == {o, unitv(cfg[0])}})
    print(f'    type-B configurations with {{x\'-z, y\'-z}} = e: {len(ecfg)}')
    for cfg in ecfg:
        rr = [r for r in rowsB if r[0] == cfg]
        shapes = sorted({(r[2], r[3], r[4]) for r in rr})
        print(f'      cfg (k={cfg[0]}, x\'-z={cfg[1]}, y\'-z={cfg[2]}): {len(rr)} options, (|I|,|O|,|W*\\O|) in {shapes}, '
              f'max cost 2^{math.log2(max(r[6] for r in rr)):.4f}, max weight {max(r[1] for r in rr)}, '
              f'weights sum to {sum(r[1] for r in rr)}')
    cmax = max(r[6] for r in rowsB)
    at_max = sorted({r[0] for r in rowsB if r[6] == cmax})
    all_max = sorted(cfg for cfg in {r[0] for r in rowsB} if all(r[6] == cmax for r in rowsB if r[0] == cfg))
    print(f'    maximal single-design type-B cost 2^{math.log2(cmax):.4f}: attained by an option of '
          f'{len(at_max)} configurations; configurations all of whose options attain it: {len(all_max)} '
          f'(= the {{x\'-z, y\'-z}} = e set: {all_max == ecfg})')
    iat, (s, grp, cl) = ep_arg['B']
    print(f'    argmax B-clique at p = {ps[iat]}: {len(cl)} keys, value 2^{math.log2(s):.6f}; keys '
          f'(edge in the frame min(P0) = 0, (|I|,|O|,|W*\\O|), n, W, log2 W c, #producing configurations, '
          f'#of them with {{x\'-z, y\'-z}} = e):')
    for nd in sorted(cl, key=lambda nd: -nd['W'] * nd['c'][iat]):
        ncf = len(nd['cfgw'])
        ne_ = sum(1 for cfg in nd['cfgw'] if cfg in ecfg)
        print(f'      e={nd["e"]} ({nd["nI"]},{nd["nO"]},{nd["nZ"]}) n={nd["n"]} W={nd["W"]} '
              f'log2(Wc)={math.log2(nd["W"] * nd["c"][iat]):.4f} cfgs={ncf} e-cfgs={ne_}')

    # ---- joint B + N groups near the origin
    bshape = {}
    for grp, dct in nodes['B'].items():
        bshape[grp] = list(dct.values())
    KBN, nj = Fr(0), 0
    for grp, dct in nodes['N'].items():
        P, ip = grp
        m = min(P)
        shape = (frozenset(sub(v, m) for v in P), frozenset(sub(v, m) for v in ip))
        joint = list(dct.values())
        for nd in bshape.get(shape, []):
            Lam = frozenset(add(v, m) for v in nd['Lam'])
            if o in Lam:
                continue                                   # that edge is of type N, not B
            tr = lambda v: add(v, m)
            joint.append(dict(nd, pat={edge(tr(f[0]), tr(f[1])): b for f, b in nd['pat'].items()},
                              variants={tuple(tr(u) for u in var) for var in nd['variants']}, Lam=Lam,
                              e=edge(tr(nd['e'][0]), tr(nd['e'][1]))))
        nj += len(joint)
        v, _, _, _ = max_clique_value(joint, ps)
        KBN = max(KBN, v)
    Ks['BN'] = KBN
    print(f'  joint B+N groups at the origin ({len(nodes["N"])} labelled plaquettes, {nj} keys): K_BN = 2^{math.log2(KBN):.4f}')

    # ---- type S: own rule, verified at each L; crude multiplicity bound
    KS = Fr(0)
    worstS_all = Fr(0)
    pictures = {}
    for L in Ls:
        cf = configs_S(d, Lam0, L)
        keys = defaultdict(dict)
        worstS = Fr(0)
        picset = set()
        W0 = sorted(Lp0)
        for z, k, xp in cf:
            b = ruleS(d, Lam0, L, z, k, xp, p_lo, p_hi, rho)
            assert b is not None, ('no type-S design', z, k, xp)
            val, rel, D, r = b
            worstS = max(worstS, val)
            worstS_all = max(worstS_all, val)
            e = edge(z, add(z, unit(d, k)))
            grp = (r['P'], r['ip'])
            pat = frozenset((f, 1 if f in r['O'] else 0) for f in r['W'])
            keys[grp][(e, pat)] = val
            pic = tuple(min(2, max(0, linf(add(z, v)) - (L - 1))) for v in W0)
            picset.add((pic, k, sub(xp, z), rel))
        crude = max(sum(dct.values()) for dct in keys.values())
        maxkeys = max(len(dct) for dct in keys.values())
        pictures[L] = picset
        KS = max(KS, crude)
        print(f'  type S, L={L}: {len(cf)} configurations, each with a verified design; {len(picset)} distinct '
              f'(picture, k, x\'-z, design); worst key 2^{math.log2(worstS):.4f}; max #keys per labelled P0 = {maxkeys}; '
              f'K_S <= sum over keys = 2^{math.log2(crude):.4f}')
    if len(Ls) > 1:
        same = all(pictures[L] == pictures[Ls[0]] for L in Ls)
        print(f'  type S: identical (picture, k, x\'-z, design) sets for L in {Ls}: {same}')
        assert same
    Ks['S'] = KS

    K = max(Ks.values())
    print(f'  K = max(K_B, K_BN, K_S) = 2^{math.log2(K):.4f}   [K_B 2^{math.log2(Ks["B"]):.4f}, K_N 2^{math.log2(Ks["N"]):.4f}, '
          f'K_BN 2^{math.log2(Ks["BN"]):.4f}, K_S 2^{math.log2(Ks["S"]):.4f}]')
    Delta = window(q, NQ, K, p_hi)
    s, _ = floor_sig(Delta)
    print(f'  plaquette coins: Delta = (1-2^-{q+1}) / ({q+1} (K + {NQ}/(2(1-p_hi)))) >= {s};  g = 2 Delta >= {floor_sig(2*Delta)[0]}')
    if coins:
        qc, NQc, name = coins
        Dc = window(qc, NQc, K, p_hi)
        print(f'  {name}: Delta_c = (1-2^-{qc+1}) / ({qc+1} (K + {NQc}/(2(1-p_hi)))) >= {floor_sig(Dc)[0]};  '
              f'g = 2 Delta_c >= {floor_sig(2*Dc)[0]}')
    print('  exact K (numerator/denominator digits):', len(str(K.numerator)), len(str(K.denominator)))
    # For the displayed d=2 constant: the window if K were replaced by the value of the binding
    # family at one endpoint only -- for comparison only; the certified window uses K = max over both endpoints.
    for i, pe in enumerate(ps):
        Ki = max(ep_vals["B"][i], ep_vals["N"][i])
        if Ki < K:
            msg = (f'  comparison only: with K replaced by max(K_B, K_N) at p = {pe} alone (2^{math.log2(Ki):.6f} < K), '
                   f'Delta would be {floor_sig(window(q, NQ, Ki, p_hi))[0]}')
            if coins:
                msg += f', coin Delta {floor_sig(window(coins[0], coins[1], Ki, p_hi))[0]}'
            print(msg + ' (too large, not claimed)')
    return dict(d=d, p_lo=p_lo, p_hi=p_hi, q=q, NQ=NQ, K=K, Ks=Ks, ep_vals=ep_vals, worst_key=worst,
                worstS=worstS_all, Delta=Delta, coins=coins,
                Delta_coins=window(coins[0], coins[1], K, p_hi) if coins else None)


if __name__ == '__main__':
    main()

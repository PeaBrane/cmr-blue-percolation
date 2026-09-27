"""Extra independent checks for the clean-room verification (cr_check.py).  Sanity only; needs networkx.

  (1) maximal cliques: own Bron-Kerbosch (cr_check.maximal_cliques) versus networkx.find_cliques, all groups;
  (2) the anchor-variant issue: compatibility graphs built with ALL anchor variants of a key (correct, Definition 5.11)
      versus only the first variant (the one-anchor bookkeeping of the program that produced the rule tables);
      resulting K_N;
  (3) every type-N design of the table has 0 as a corner of P0 (used in Section 5.7);
  (4) negative controls: mutated designs must be rejected by cr_check.check_design.
Usage: python audit_extra.py <rule.json>
"""
import sys
import json
import math
import random
import itertools
from fractions import Fraction as Fr
from collections import defaultdict
import networkx as nx
import cr_check as C


def build(data, first_variant_only=False):
    d, kind = data['d'], data['region']
    p_lo, p_hi = Fr(data['p_lo']), Fr(data['p_hi'])
    rho = 1 - (1 - p_hi) ** (2 * d - 1)
    o = tuple([0] * d)
    Lam0 = C.region0(d, kind)
    nodes = {'B': defaultdict(dict), 'N': defaultdict(dict)}
    for typ in ('B', 'N'):
        for entry in data['rule' + typ]:
            c = entry['config']
            if typ == 'B':
                k, xp, yp = c[0], tuple(c[1]), tuple(c[2])
                z, Lam = o, set(Lam0)
            else:
                z, k, yp = tuple(c[0]), c[1], tuple(c[2])
                xp, Lam = None, {C.add(z, v) for v in Lam0}
            e = C.edge(z, C.add(z, C.unit(d, k)))
            wkey = defaultdict(Fr)
            for opt in entry['options']:
                D, w = opt['design'], Fr(opt['w'])
                if typ == 'B':
                    r = C.check_design('B', e, Lam, D, xp=xp, yp=yp, src=None)
                    t = C.sub(o, min(r['P']))
                else:
                    st = 'Nint' if o in {tuple(v) for v in D['intp']} else 'Nter'
                    r = C.check_design(st, e, Lam, D, yp=yp, src=o)
                    t = o
                tr = lambda v, t=t: C.add(v, t)
                grp = (frozenset(tr(v) for v in r['P']), frozenset(tr(v) for v in r['ip']))
                pat = {C.edge(tr(f[0]), tr(f[1])): (1 if f in r['O'] else 0) for f in r['W']}
                key = (C.edge(tr(e[0]), tr(e[1])), frozenset(pat.items()))
                nO, nZ = len(r['O']), len(r['W']) - len(r['O'])
                cs = tuple(C.cost(nO, nZ, r['n'], p, rho) for p in (p_lo, p_hi))
                var = tuple(sorted(tr(u) for u in r['anchors']))
                nd = nodes[typ][grp].get(key)
                if nd is None:
                    nodes[typ][grp][key] = dict(pat=pat, variants={var}, Lam=frozenset(tr(v) for v in Lam), c=cs, W=Fr(0))
                elif not first_variant_only:
                    nd['variants'].add(var)
                wkey[(grp, key)] += w
            for (grp, key), w in wkey.items():
                nodes[typ][grp][key]['W'] = max(nodes[typ][grp][key]['W'], w)
    return nodes, (p_lo, p_hi)


def graph(lst):
    n = len(lst)
    adj = [set() for _ in range(n)]
    for i, j in itertools.combinations(range(n), 2):
        if C.compatible(lst[i], lst[j]):
            adj[i].add(j)
            adj[j].add(i)
    return adj


def main():
    data = json.load(open(sys.argv[1]))
    d = data['d']
    o = tuple([0] * d)
    nodes, ps = build(data)
    # (1) cliques: own vs networkx
    ngr, ncl = 0, 0
    for typ in ('B', 'N'):
        for grp, dct in nodes[typ].items():
            lst = list(dct.values())
            adj = graph(lst)
            own = {frozenset(cl) for cl in C.maximal_cliques(len(lst), adj)}
            G = nx.Graph()
            G.add_nodes_from(range(len(lst)))
            G.add_edges_from((i, j) for i in range(len(lst)) for j in adj[i] if i < j)
            ref = {frozenset(cl) for cl in nx.find_cliques(G)}
            assert own == ref, ('clique mismatch', typ)
            ngr += 1
            ncl += len(own)
    print(f'(1) maximal cliques: own Bron-Kerbosch = networkx.find_cliques on all {ngr} groups ({ncl} cliques)')
    # (2) anchor variants
    nodes1, _ = build(data, first_variant_only=True)
    for label, nd_ in (('all anchor variants (correct)', nodes), ('first variant only (one-anchor bookkeeping)', nodes1)):
        K, npairs = Fr(0), 0
        for grp, dct in nd_['N'].items():
            lst = list(dct.values())
            v, _, ne, _ = C.max_clique_value(lst, ps)
            K = max(K, v)
            npairs += ne
        print(f'(2) type N, {label}: {npairs} compatible pairs, K_N = 2^{math.log2(K):.4f}')
    # (3) N designs contain 0
    tot = sum(len(ent['options']) for ent in data['ruleN'])
    with0 = sum(1 for ent in data['ruleN'] for opt in ent['options'] if o in {tuple(v) for v in opt['design']['P']})
    print(f'(3) type-N options with 0 in V(P0): {with0} of {tot}')
    # (4) negative controls
    random.seed(1)
    Lam0 = C.region0(d, data['region'])
    rejected, total = defaultdict(int), defaultdict(int)
    for entry in data['ruleB']:
        k, xp, yp = entry['config'][0], tuple(entry['config'][1]), tuple(entry['config'][2])
        e = C.edge(o, C.unit(d, k))
        for opt in entry['options']:
            D = opt['design']
            muts = {}
            muts['swap labels'] = dict(D, intp=D['terp'], terp=D['intp'])
            A = [list(map(tuple, a)) for a in D['arcs']]
            ip = [tuple(v) for v in D['intp']]
            muts['arc through internal vertex'] = dict(D, arcs=[a[:-1] + [ip[0]] + [a[-1]] if len(a) >= 1 else a for a in A[:1]] + A[1:])
            muts['arcs swapped attachments'] = dict(D, arcs=[A[1], A[0]])
            muts['arc ends at internal vertex'] = dict(D, arcs=[A[0][:-1] + [ip[0]]] + A[1:])
            muts['drop second arc'] = dict(D, arcs=A[:1])
            for name, M in muts.items():
                total[name] += 1
                try:
                    C.check_design('B', e, set(Lam0), M, xp=xp, yp=yp, src=None)
                except (C.Bad, Exception):
                    rejected[name] += 1
    for name in total:
        print(f'(4) mutation "{name}": rejected {rejected[name]} of {total[name]}')


if __name__ == '__main__':
    main()

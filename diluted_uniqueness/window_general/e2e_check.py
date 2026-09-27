"""End-to-end sanity test of the sparse pendant-arc surgery map and its counting; not part of the proof.

On random states of the diminished model (plaquette coins) in B_{L+2}, with a planted self-avoiding path from 0 to
S_L, every (+)-pivotal edge found among the candidates is processed literally:
  Step 1 (deactivate the activated plaquettes of Q(e) one by one), else Step 2: Pi = BFS-first shortest path in
  E(omega^{e,1}, alpha), type / attachments / configuration, and EVERY option of the rule (types B, N from the JSON
  table, type S from cr_check.ruleS).  For each option: hypotheses (cr_check.check_design), then the output is
  built and checked: P0 is an isolated diamond of omega' with the designed internal pair, (omega', alpha) in A_L,
  (omega', alpha^{P0,1}) not in A_L, the anchors are omega'-open, and the rho-event holds for omega.
Usage: python e2e_check.py <rule.json> L n_states seed p s
"""
import sys
import json
import random
import itertools
from collections import deque
import cr_check as C


def main():
    data = json.load(open(sys.argv[1]))
    L, nstates, seed, p, s = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
    d, kind = data['d'], data['region']
    from fractions import Fraction as Fr
    p_lo, p_hi = Fr(data['p_lo']), Fr(data['p_hi'])
    rho = 1 - (1 - p_hi) ** (2 * d - 1)
    o = tuple([0] * d)
    Lam0 = C.region0(d, kind)
    Lp0 = set(Lam0) | {u for v in Lam0 for u in C.nbrs(v)}
    ruleB = {(c['config'][0], tuple(c['config'][1]), tuple(c['config'][2])): c['options'] for c in data['ruleB']}
    ruleN = {(tuple(c['config'][0]), c['config'][1], tuple(c['config'][2])): c['options'] for c in data['ruleN']}
    R = L + 2
    V = [v for v in itertools.product(range(-R, R + 1), repeat=d)]
    Vset = set(V)
    edges = sorted({C.edge(u, w) for u in V for w in C.nbrs(u) if w in Vset})
    plaqs = [P for P in {frozenset(P) for u in V for P in C.plaquettes_with_corner(u)} if all(C.linf(v) <= L + 1 for v in P)]
    plaqs.sort(key=lambda P: sorted(P))
    pedges = {P: C.plaq_edges(P) for P in plaqs}
    by_edge = {}
    for P in plaqs:
        for f in pedges[P]:
            by_edge.setdefault(f, []).append(P)
    rng = random.Random(seed)
    stats = dict(states=0, cand=0, piv=0, step1=0, step2=0, B=0, N=0, S=0, options=0, fail=0)

    def Eopen(om, al):
        deg = {}
        for f in om:
            for v in f:
                deg[v] = deg.get(v, 0) + 1
        removed = set()
        seenP = set()
        for f in om:
            for P in by_edge.get(f, ()):
                if P in seenP or P not in al:
                    continue
                seenP.add(P)
                if not all(g in om for g in pedges[P]):
                    continue
                d1, d2 = C.diag_pairs(P)
                two1 = all(deg.get(v, 0) == 2 for v in d1)
                two2 = all(deg.get(v, 0) == 2 for v in d2)
                if two1 != two2:
                    removed |= pedges[P]
        return om - removed

    def is_diamond(om, P, ip):
        if not all(g in om for g in pedges.get(P, C.plaq_edges(P))):
            return False
        deg = lambda v: sum(1 for w in C.nbrs(v) if C.edge(v, w) in om)
        tp = frozenset(P) - ip
        return all(deg(v) == 2 for v in ip) and not all(deg(v) == 2 for v in tp)

    def bfs_path(E):
        adj = {}
        for u, w in E:
            adj.setdefault(u, []).append(w)
            adj.setdefault(w, []).append(u)
        par = {o: None}
        dq = deque([o])
        while dq:
            v = dq.popleft()
            if C.linf(v) == L:
                path = [v]
                while par[path[-1]] is not None:
                    path.append(par[path[-1]])
                return path[::-1]
            for w in sorted(adj.get(v, ())):
                if w not in par:
                    par[w] = v
                    dq.append(w)
        return None

    def AL(om, al):
        return bfs_path(Eopen(om, al)) is not None

    for st in range(nstates):
        om = {f for f in edges if rng.random() < p}
        # planted self-avoiding path from 0 to S_L (random walk avoiding itself, restart if stuck)
        while True:
            path, seen = [o], {o}
            while C.linf(path[-1]) < L:
                opts = [w for w in C.nbrs(path[-1]) if w not in seen and C.linf(w) <= L]
                if not opts:
                    break
                far = [w for w in opts if C.linf(w) >= C.linf(path[-1])]
                w = rng.choice(far if far and rng.random() < 0.6 else opts)
                path.append(w)
                seen.add(w)
            if C.linf(path[-1]) == L:
                break
        om |= {C.edge(a, b) for a, b in zip(path, path[1:])}
        al = {P for P in plaqs if rng.random() < s}
        stats['states'] += 1
        cands = {C.edge(a, b) for a, b in zip(path, path[1:])}
        comp = bfs_path(Eopen(om, al))
        cands |= set(rng.sample(edges, 40))
        for e in sorted(cands):
            if not (C.linf(e[0]) <= L and C.linf(e[1]) <= L):
                continue
            stats['cand'] += 1
            if not AL(om | {e}, al) or AL(om - {e}, al):
                continue
            stats['piv'] += 1
            z = e[0]
            k = next(i for i in range(d) if e[1][i] != e[0][i])
            Lam = {C.add(z, v) for v in Lam0}
            Lp = {C.add(z, v) for v in Lp0}
            QE = [P for P in plaqs if P & Lp]
            acts = [P for P in QE if P in al]
            alj = set(al)
            done = False
            for P in acts:
                alj.discard(P)
                if AL(om - {e}, alj):
                    stats['step1'] += 1
                    done = True
                    break
            if done:
                continue
            stats['step2'] += 1
            alk = set(al) - set(QE)
            Pi = bfs_path(Eopen(om | {e}, alk))
            assert Pi is not None
            ePi = {C.edge(a, b) for a, b in zip(Pi, Pi[1:])}
            assert e in ePi, 'Pi does not use e'
            idx = [i for i, v in enumerate(Pi) if v in Lam]
            i1, i2 = idx[0], idx[-1]
            if o in Lam:
                typ = 'N'
                yp = Pi[i2]
                opts = [(opt['design'], None) for opt in ruleN[(z, k, yp)]]
                anchors = {yp: Pi[i2 + 1]}
            elif all(C.linf(v) <= L - 1 for v in Lam):
                typ = 'B'
                xp, yp = Pi[i1], Pi[i2]
                raw = ruleB[(k, C.sub(xp, z), C.sub(yp, z))]
                tr = lambda A: [C.add(tuple(v), z) for v in A]
                opts = [(dict(P=tr(opt['design']['P']), intp=tr(opt['design']['intp']), terp=tr(opt['design']['terp']),
                              arcs=[tr(A) for A in opt['design']['arcs']]), None) for opt in raw]
                anchors = {xp: Pi[i1 - 1], yp: Pi[i2 + 1]}
            else:
                typ = 'S'
                xp = Pi[i1]
                b = C.ruleS(d, Lam0, L, z, k, xp, p_lo, p_hi, rho)
                opts = [(b[2], None)]
                anchors = {xp: Pi[i1 - 1]}
            stats[typ] += 1
            for D, _ in opts:
                stats['options'] += 1
                if typ == 'B':
                    r = C.check_design('B', e, Lam, D, xp=xp, yp=yp, src=o)
                elif typ == 'N':
                    st_ = 'Nint' if o in {tuple(v) for v in D['intp']} else 'Nter'
                    r = C.check_design(st_, e, Lam, D, yp=yp, src=o)
                else:
                    r = C.check_design('S', e, Lam, D, L=L, xp=xp, src=o)
                omp = (om - r['W']) | r['O']
                P0, ip = r['P'], r['ip']
                ok = is_diamond(omp, P0, ip)
                ok &= AL(omp, alk)
                ok &= not AL(omp, alk | {P0})
                ok &= all(C.edge(u, anchors[u]) in omp for u in r['anchors'])
                for w in e:
                    if w in r['I'] and w != o and C.linf(w) < L:
                        ok &= any(f in om for f in C.Fset([w]) if f != e)
                if not ok:
                    stats['fail'] += 1
                    print('FAIL', typ, e, D, flush=True)
        if (st + 1) % 10 == 0:
            print(st + 1, stats, flush=True)
    print('final', stats, flush=True)


if __name__ == '__main__':
    main()

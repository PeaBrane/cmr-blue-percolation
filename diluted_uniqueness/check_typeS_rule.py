"""Sanity check (not part of any proof) of the type-S surgery rule of Claim D, as written in the
manuscript (research-note source: theorem-A.md, Sec. 4.6, revision 2), in original coordinates.

For (d, L) = (2, 10), (2, 11), (2, 13) and (3, 10): every centre z with L - rho <= |z| <= L, every entry vertex x' in R = Lambda cap B_{L-1}
that has a neighbour x outside Lambda, and every such x. Build (i*, sigma*), j*, c, u, a, b, w, seg_1, O as in the
text and check the local facts used in the proof:
  (S1) c exists; u, a in R; b, w in Lambda cap S_L;
  (S2) seg_1 is a nearest-neighbour path, self-avoiding, contained in R, avoiding a, b, w;
  (S3) seg_1 meets V(P_0) only in u and no edge of seg_1 or xx' is a P_0-edge;
  (S4) O-degrees: a, b, w have degree 2, u has degree 3 (so (R1) holds with internal pair {a, b});
  (S5) O minus E(P_0) is exactly the edge set of the arc (x, x', ..., u), whose vertices other than x lie in
       Lambda minus S_L (the hypothesis of Lemma 4.8);
  (S6) O connects x to S_L (via u -> b), the local part of (R2).
"""
import itertools
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

RHO = 4


def ninf(v):
    return max(abs(t) for t in v)


def add(v, i, s):
    w = list(v)
    w[i] += s
    return tuple(w)


def run(d, L):
    count = 0
    rng = range(-L, L + 1)
    for z in itertools.product(rng, repeat=d):
        if not (L - RHO <= ninf(z) <= L):
            continue
        def in_lam(v):
            return all(abs(v[k] - z[k]) <= RHO for k in range(d))
        def in_R(v):
            return in_lam(v) and ninf(v) <= L - 1
        pairs = [(i, s) for i in range(d) for s in (1, -1)]
        istar, sg = next((i, s) for (i, s) in pairs if s * z[i] >= L - RHO)
        jstar = 0 if istar != 0 else 1
        I = [(max(z[k] - RHO, -(L - 1)), min(z[k] + RHO, L - 1)) for k in range(d)]
        for xp in itertools.product(*[range(lo, hi + 1) for (lo, hi) in I]):
            if not in_R(xp):
                continue
            outs = [add(xp, i, s) for i in range(d) for s in (1, -1) if not in_lam(add(xp, i, s))]
            if not outs:
                continue
            lo, hi = I[jstar]
            cands = [c for c in range(lo, hi) if c not in (xp[jstar] - 1, xp[jstar])]
            assert cands, ('no c', z, xp)
            c = cands[0]
            u = list(xp); u[istar] = sg * (L - 1); u[jstar] = c; u = tuple(u)
            a = add(u, jstar, 1)
            b = add(u, istar, sg)
            w = add(a, istar, sg)
            assert in_R(u) and in_R(a), 'S1 u,a'
            assert in_lam(b) and in_lam(w) and ninf(b) == L and ninf(w) == L, 'S1 b,w'
            # seg_1
            seg = [xp]
            cur = list(xp)
            target = sg * (L - 2)
            while cur[istar] != target:
                cur[istar] += 1 if target > cur[istar] else -1
                seg.append(tuple(cur))
            while cur[jstar] != c:
                cur[jstar] += 1 if c > cur[jstar] else -1
                seg.append(tuple(cur))
            cur[istar] = sg * (L - 1)
            seg.append(tuple(cur))
            assert seg[-1] == u
            assert all(sum(abs(p - q) for p, q in zip(seg[k], seg[k + 1])) == 1 for k in range(len(seg) - 1)), 'S2 nn'
            assert len(set(seg)) == len(seg), 'S2 self-avoiding'
            assert all(in_R(v) for v in seg), 'S2 in R'
            assert not ({a, b, w} & set(seg)), 'S2 avoid a,b,w'
            P0 = {u, a, w, b}
            P0_edges = {frozenset((u, a)), frozenset((a, w)), frozenset((w, b)), frozenset((b, u))}
            assert set(seg) & P0 == {u}, 'S3'
            seg_edges = {frozenset((seg[k], seg[k + 1])) for k in range(len(seg) - 1)}
            for x in outs:
                O = P0_edges | seg_edges | {frozenset((x, xp))}
                assert not ((seg_edges | {frozenset((x, xp))}) & P0_edges), 'S3 edges'
                deg = {}
                for e in O:
                    for v in e:
                        deg[v] = deg.get(v, 0) + 1
                assert deg[a] == 2 and deg[b] == 2 and deg[w] == 2 and deg[u] == 3, ('S4', deg[a], deg[b], deg[w], deg[u])
                arc = [x] + seg
                arc_edges = {frozenset((arc[k], arc[k + 1])) for k in range(len(arc) - 1)}
                assert O - P0_edges == arc_edges, 'S5 exact'
                assert len(set(arc)) == len(arc), 'S5 arc simple'
                assert all(in_lam(v) and ninf(v) <= L - 1 for v in arc[1:]) and not in_lam(x), 'S5 arc position'
                # S6: x -> ... -> u -> b in S_L inside O
                assert frozenset((u, b)) in O and ninf(b) == L
                count += 1
    return count


if __name__ == '__main__':
    # (d, L, number of configurations displayed in the manuscript)
    for d, L, shown in ((2, 10, 5820), (2, 11, 6660), (2, 13, 8340), (3, 10, 1780110)):
        n = run(d, L)
        assert n == shown, (d, L, n)
        print(f'd={d} L={L}: {n} (z, x\', x) configurations, all checks S1-S6 passed')
        sys.stdout.flush()

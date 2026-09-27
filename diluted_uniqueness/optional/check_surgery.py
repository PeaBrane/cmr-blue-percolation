"""Sanity check (NOT part of the proof) of the surgery map Phi (the map S of Proposition 4.7 of the diluted-model
manuscript; optional; needs networkx).

Diluted-model manuscript, Section 4 and Appendix A.3.

For random states (omega, alpha) on Z^d (d = 2, 3), L = 10, rho = 4, we enumerate ALL edges e that are
(+)-pivotal for A_L = {0 <-> S_L in E(omega, alpha)}, apply the map Phi exactly as specified in the
document (Step 1 deactivation in Q, then Step 2 surgery of type N / B / S), and verify:
  * the output plaquette P is s-pivotal in the output state;
  * omega' = omega off F(e), alpha' = alpha off Q(e), alpha' <= alpha, P in Q(e);
  * in Step 2: e is (+)-pivotal after deactivation, the chosen path uses e, the designed plaquette P_0 is a
    non-4-cycle isolated diamond of omega' with the designed internal pair, the two segments are disjoint.
Candidate edges: all edges with an endpoint within l_inf-distance 3 of the E-cluster of 0 (this contains every
(+)-pivotal edge, see the comment in plus_pivotal_candidates).
Usage: python check_surgery.py d nsamples seed [random|path] [L]
"""
import sys
import random
import itertools
from collections import deque
import networkx as nx

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

RHO = 4
L = 10


def add(v, w):
    return tuple(a + b for a, b in zip(v, w))


def unit(d, i, s=1):
    u = [0] * d
    u[i] = s
    return tuple(u)


def ninf(v):
    return max(abs(a) for a in v)


class State:
    def __init__(self, d, omega, alpha):
        self.d = d
        self.omega = omega      # dict edge -> 0/1, edge = (v, i) meaning {v, v+e_i}
        self.alpha = alpha      # dict plaquette -> 0/1, plaquette = (v, i, j), i<j

    def w(self, e, ov):
        if e in ov:
            return ov[e]
        return self.omega.get(e, 0)

    def a(self, P, av):
        if P in av:
            return av[P]
        return self.alpha.get(P, 0)


def edge_key(u, v):
    """canonical key of the edge {u,v}"""
    diff = [b - a for a, b in zip(u, v)]
    i = [k for k, x in enumerate(diff) if x != 0][0]
    if diff[i] == 1:
        return (u, i)
    return (v, i)


def endpoints(e, d):
    v, i = e
    return v, add(v, unit(d, i))


def edges_at(v, d):
    out = []
    for i in range(d):
        out.append((v, i))
        out.append((add(v, unit(d, i, -1)), i))
    return out


def corners(P, d):
    v, i, j = P
    ei, ej = unit(d, i), unit(d, j)
    return [v, add(v, ei), add(add(v, ei), ej), add(v, ej)]  # cyclic order


def plaq_edges(P, d):
    c = corners(P, d)
    return [edge_key(c[k], c[(k + 1) % 4]) for k in range(4)]


def plaqs_containing(e, d):
    v, i = e
    out = []
    for j in range(d):
        if j == i:
            continue
        a, b = min(i, j), max(i, j)
        out.append((v, a, b))
        out.append((add(v, unit(d, j, -1)), a, b))
    return out


def plaqs_at_vertex(v, d):
    out = set()
    for e in edges_at(v, d):
        for P in plaqs_containing(e, d):
            out.add(P)
    return out


def deg(st, v, ov):
    return sum(st.w(e, ov) for e in edges_at(v, st.d))


def diamond_info(st, P, ov):
    """None if P is not a non-4-cycle isolated diamond, else (internal pair, terminal pair)."""
    d = st.d
    if not all(st.w(e, ov) for e in plaq_edges(P, d)):
        return None
    c = corners(P, d)
    dg = [deg(st, x, ov) for x in c]
    q02 = dg[0] == 2 and dg[2] == 2
    q13 = dg[1] == 2 and dg[3] == 2
    if q02 == q13:
        return None
    if q02:
        return (c[0], c[2]), (c[1], c[3])
    return (c[1], c[3]), (c[0], c[2])


def E_open(st, e, ov, av):
    if not st.w(e, ov):
        return False
    for P in plaqs_containing(e, st.d):
        if st.a(P, av) and diamond_info(st, P, ov) is not None:
            return False
    return True


def nbrs(v, d):
    for i in range(d):
        for s in (1, -1):
            yield add(v, unit(d, i, s))


def A_L(st, ov=None, av=None, want_path=False):
    ov = ov or {}
    av = av or {}
    d = st.d
    o = tuple([0] * d)
    par = {o: None}
    dq = deque([o])
    while dq:
        v = dq.popleft()
        if ninf(v) == L:
            if not want_path:
                return True
            path = []
            while v is not None:
                path.append(v)
                v = par[v]
            return path[::-1]
        for u in nbrs(v, d):
            if ninf(u) > L or u in par:
                continue
            if E_open(st, edge_key(v, u), ov, av):
                par[u] = v
                dq.append(u)
    return False if not want_path else None


def cluster0(st):
    d = st.d
    o = tuple([0] * d)
    seen = {o}
    dq = deque([o])
    while dq:
        v = dq.popleft()
        for u in nbrs(v, d):
            if ninf(u) > L or u in seen:
                continue
            if E_open(st, edge_key(v, u), {}, {}):
                seen.add(u)
                dq.append(u)
    return seen


def plus_pivotal_candidates(st):
    # If e is (+)-pivotal, a path from 0 to S_L in E(omega^{e,1}, alpha) must use an edge whose E-status differs
    # between omega^{e,1} and omega^{e,0}; such edges lie on plaquettes with a corner at an endpoint of e.  The part
    # of the path before the first such edge lies in the E-cluster of 0 for whichever of omega^{e,0}, omega^{e,1}
    # equals omega; hence e has an endpoint within l_inf distance 2 of cluster0 (we use 3 for margin).
    d = st.d
    cl = cluster0(st)
    cand = set()
    for v in cl:
        for off in itertools.product(range(-3, 4), repeat=d):
            u = add(v, off)
            if ninf(u) <= L + 1:
                for e in edges_at(u, d):
                    cand.add(e)
    return sorted(cand)


def is_plus_pivotal(st, e, av=None):
    return bool(A_L(st, {e: 1}, av)) and not A_L(st, {e: 0}, av)


def is_s_pivotal(st, P, ov=None, av=None):
    ov = ov or {}
    av = dict(av or {})
    av0 = dict(av); av0[P] = 0
    av1 = dict(av); av1[P] = 1
    return bool(A_L(st, ov, av0)) and not A_L(st, ov, av1)


def box(c, r, d):
    return [add(c, off) for off in itertools.product(range(-r, r + 1), repeat=d)]


def lam_of(e, d):
    """surgery box (center, radius, type) determined by e alone"""
    z = endpoints(e, d)[0]          # lexicographically smaller endpoint
    o = tuple([0] * d)
    if ninf(z) <= RHO:
        return o, RHO + 1, 'N', z
    if ninf(z) >= L - RHO:
        return z, RHO, 'S', z
    return z, RHO, 'B', z


def F_of(c, r, d):
    F = set()
    for v in box(c, r, d):
        for e in edges_at(v, d):
            F.add(e)
    return F


def Q_of(c, r, d):
    Q = set()
    for v in box(c, r + 1, d):
        for P in plaqs_at_vertex(v, d):
            Q.add(P)
    return Q


def in_box(v, c, r):
    return max(abs(a - b) for a, b in zip(v, c)) <= r


def shell_graph(c, r, d):
    S = [v for v in box(c, r, d) if max(abs(a - b) for a, b in zip(v, c)) == r]
    Sset = set(S)
    G = nx.Graph()
    G.add_nodes_from(S)
    for v in S:
        for u in nbrs(v, d):
            if u in Sset:
                G.add_edge(v, u)
    return G


def cycle_pairing_2d(c, r, X, T):
    """explicit d=2 cycle-pairing lemma: returns {x: path from x to its target} (vertex-disjoint)"""
    # order the ring R_r(c) cyclically (counterclockwise from (c0+r, c1-r))
    cx, cy = c
    ring = []
    for k in range(-r, r):
        ring.append((cx + r, cy + k))
    for k in range(r, -r, -1):
        ring.append((cx + k, cy + r))
    for k in range(r, -r, -1):
        ring.append((cx - r, cy + k))
    for k in range(-r, r):
        ring.append((cx + k, cy - r))
    n = len(ring)
    assert n == 8 * r and len(set(ring)) == n
    pos = {v: k for k, v in enumerate(ring)}
    x1, x2 = X
    t1, t2 = T
    if set(X) == set(T):
        return {x1: [x1], x2: [x2]}
    for (xa, xb) in ((x1, x2), (x2, x1)):
        for (ta, tb) in ((t1, t2), (t2, t1)):
            if xa == ta:
                # path from xb to tb avoiding xa: go the way around that avoids xa
                for step in (1, -1):
                    path = [xb]
                    k = pos[xb]
                    ok = True
                    while ring[k] != tb:
                        k = (k + step) % n
                        if ring[k] == xa:
                            ok = False
                            break
                        path.append(ring[k])
                    if ok:
                        return {xa: [xa], xb: path}
    # X and T disjoint: try all four (pairing, direction) combinations, pick the first vertex-disjoint one
    for (ta, tb) in ((t1, t2), (t2, t1)):
        for s1 in (1, -1):
            for s2 in (1, -1):
                def walk(start, target, step, forbid):
                    path = [start]
                    k = pos[start]
                    while ring[k] != target:
                        k = (k + step) % n
                        if ring[k] in forbid:
                            return None
                        path.append(ring[k])
                    return path
                p1 = walk(x1, ta, s1, {x2, tb})
                p2 = walk(x2, tb, s2, {x1, ta})
                if p1 and p2 and not (set(p1) & set(p2)):
                    return {x1: p1, x2: p2}
    raise RuntimeError('cycle pairing failed')


def menger_pairing(G, X, T):
    H = G.copy()
    H.add_node('s'); H.add_node('t')
    for x in X:
        H.add_edge('s', x)
    for t in T:
        H.add_edge(t, 't')
    paths = list(nx.node_disjoint_paths(H, 's', 't'))
    assert len(paths) == 2
    out = {}
    for p in paths:
        q = p[1:-1]
        assert q[0] in X and q[-1] in T
        assert all(v not in X for v in q[1:]) and all(v not in T for v in q[:-1])
        out[q[0]] = q
    return out


def path_edges(path):
    return [edge_key(path[k], path[k + 1]) for k in range(len(path) - 1)]


def phi(st, e, stats):
    d = st.d
    c, r, typ, z = lam_of(e, d)
    F = F_of(c, r, d)
    Q = Q_of(c, r, d)
    assert e in F
    # ---- Step 1: deactivation (fixed order: sorted)
    acts = sorted(P for P in Q if st.alpha.get(P, 0))
    av = {}
    for P in acts:
        av[P] = 0
        if A_L(st, {e: 0}, av):
            stats['step1'] += 1
            return dict(omega_changes={e: 0}, alpha_changes=dict(av), P=P, F=F, Q=Q, typ='step1')
    # ---- Step 2: surgery; now no activation in Q
    assert is_plus_pivotal(st, e, av), 'e not (+)-pivotal after deactivation'
    Pi = A_L(st, {e: 1}, av, want_path=True)
    assert Pi is not None
    Pi_edges = set(path_edges(Pi))
    assert e in Pi_edges, 'Pi does not use e'
    inL = [k for k, v in enumerate(Pi) if in_box(v, c, r)]
    assert len(inL) >= 2
    i1, i2 = inL[0], inL[-1]
    newF = {}          # the designed configuration on F (all other F-edges closed)
    if typ == 'N':
        assert i1 == 0 and i2 < len(Pi) - 1
        yp, y = Pi[i2], Pi[i2 + 1]
        o = tuple([0] * d)
        e1, e2 = unit(d, 0), unit(d, 1)
        P0 = (o, 0, 1)
        internal = (e1, e2)
        wv = add(e1, e2)
        line = [add(add(e2, e1), unit(d, 0, k)) for k in range(0, r)]   # e1+e2 ... r e1 + e2
        T = line[-1]
        G = shell_graph(o, r, d)
        sp_ = nx.shortest_path(G, T, yp)
        seg2 = line + sp_[1:]
        assert len(set(seg2)) == len(seg2)
        open_edges = set(plaq_edges(P0, d)) | set(path_edges(seg2)) | {edge_key(yp, y)}
        seg1 = [o]
        stats['N'] += 1
    elif typ == 'B':
        assert i1 > 0 and i2 < len(Pi) - 1
        x, xp, yp, y = Pi[i1 - 1], Pi[i1], Pi[i2], Pi[i2 + 1]
        e1, e2 = unit(d, 0), unit(d, 1)
        P0 = (c, 0, 1)
        internal = (add(c, e1), add(c, e2))
        TL = add(c, unit(d, 0, -r))
        TR = add(add(c, unit(d, 0, r)), e2)
        lineL = [add(c, unit(d, 0, -k)) for k in range(r, -1, -1)]          # TL ... c
        lineR = [add(add(c, e2), unit(d, 0, k)) for k in range(r, 0, -1)]   # TR ... c+e1+e2
        G = shell_graph(c, r, d)
        assert xp in G and yp in G
        if d == 2:
            link = cycle_pairing_2d(c, r, (xp, yp), (TL, TR))
            link_m = menger_pairing(G, (xp, yp), (TL, TR))
            stats['menger_crosscheck'] += 1
        else:
            link = menger_pairing(G, (xp, yp), (TL, TR))
        segs = {}
        for s0 in (xp, yp):
            pth = link[s0]
            assert pth[0] == s0 and pth[-1] in (TL, TR)
            assert all(G.has_edge(pth[k], pth[k + 1]) for k in range(len(pth) - 1))
            line = lineL if pth[-1] == TL else lineR
            segs[s0] = pth + line[1:]
        seg1, seg2 = segs[xp], segs[yp]
        assert not (set(seg1) & set(seg2))
        assert not (set(seg1 + seg2) & set(internal))
        open_edges = (set(plaq_edges(P0, d)) | set(path_edges(seg1)) | set(path_edges(seg2))
                      | {edge_key(x, xp), edge_key(yp, y)})
        if xp in (TL, TR):
            stats['B_xp_is_T'] += 1
        if sum(1 for a, b in zip(xp, c) if abs(a - b) == r) >= 2:
            stats['B_xp_on_ridge_or_corner'] += 1
        stats['B'] += 1
    else:  # 'S'
        assert i1 > 0
        x, xp = Pi[i1 - 1], Pi[i1]
        assert ninf(xp) <= L - 1
        # face: first (i*, sigma*) with sigma* z_{i*} >= L - rho
        istar, sg = next((i, s) for i in range(d) for s in (1, -1) if s * z[i] >= L - RHO)
        jstar = 0 if istar != 0 else 1
        lo = max(z[jstar] - RHO, -(L - 1))
        hi = min(z[jstar] + RHO, L - 1)
        assert hi - lo + 1 >= 4
        cc = next(t for t in range(lo, hi) if t not in (xp[jstar] - 1, xp[jstar]))
        u = list(xp); u[jstar] = cc; u[istar] = sg * (L - 1); u = tuple(u)
        a_ = add(u, unit(d, jstar))
        b_ = add(u, unit(d, istar, sg))
        P0 = (tuple(min(p, q) for p, q in zip(u, add(a_, unit(d, istar, sg)))), min(istar, jstar), max(istar, jstar))
        internal = (a_, b_)
        # segment 1: move coordinate istar to sg*(L-2), then coordinate jstar to cc, then one step to u
        seg1 = [xp]
        cur = list(xp)
        target_h = sg * (L - 2)
        while cur[istar] != target_h:
            cur[istar] += 1 if target_h > cur[istar] else -1
            seg1.append(tuple(cur))
        while cur[jstar] != cc:
            cur[jstar] += 1 if cc > cur[jstar] else -1
            seg1.append(tuple(cur))
        cur[istar] = sg * (L - 1)
        seg1.append(tuple(cur))
        assert seg1[-1] == u and len(set(seg1)) == len(seg1)
        assert all(in_box(v, c, r) and ninf(v) <= L - 1 for v in seg1)
        assert not (set(seg1) & set(internal))
        assert set(corners(P0, d)) == {u, a_, b_, add(a_, unit(d, istar, sg))}
        open_edges = set(plaq_edges(P0, d)) | set(path_edges(seg1)) | {edge_key(x, xp)}
        seg2 = []
        stats['S'] += 1
        if i2 == len(Pi) - 1:
            stats['S_path_ends_in_box'] += 1
    assert open_edges <= F
    ov = {f: (1 if f in open_edges else 0) for f in F}
    info = diamond_info(st, P0, ov)
    assert info is not None and set(info[0]) == set(internal), ('P0 not the designed isolated diamond', info)
    return dict(omega_changes=ov, alpha_changes=dict(av), P=P0, F=F, Q=Q, typ=typ)


def verify(st, e, out):
    d = st.d
    ov = out['omega_changes']
    av = out['alpha_changes']
    P = out['P']
    assert set(ov) <= out['F']
    assert set(av) <= out['Q'] and all(v == 0 for v in av.values())
    assert all(st.alpha.get(Pp, 0) == 1 for Pp in av)   # only removals of existing activations
    assert P in out['Q']
    assert av.get(P, st.alpha.get(P, 0)) == 0
    assert is_s_pivotal(st, P, ov, av), 'output plaquette not s-pivotal'


def random_state(d, p, s, rng):
    omega, alpha = {}, {}
    R = L + 4
    for v in itertools.product(range(-R, R + 1), repeat=d):
        for i in range(d):
            omega[(v, i)] = 1 if rng.random() < p else 0
        for i in range(d):
            for j in range(i + 1, d):
                alpha[(v, i, j)] = 1 if rng.random() < s else 0
    return State(d, omega, alpha)


def path_state(d, q, s, rng, plant_rate=0.1):
    """targeted generator: a random geodesic (random edge weights) from 0 to S_L, noise edges with prob q,
    diamonds planted at turns of the path (activated with prob 1/2), random alpha with prob s elsewhere."""
    import heapq
    R = L + 4
    o = tuple([0] * d)
    wts = {}
    dist = {o: 0.0}
    par = {o: None}
    pq = [(0.0, o)]
    end = None
    while pq:
        dv, v = heapq.heappop(pq)
        if dv > dist[v]:
            continue
        if ninf(v) == L:
            end = v
            break
        for u in nbrs(v, d):
            if ninf(u) > L:
                continue
            ek = edge_key(v, u)
            if ek not in wts:
                wts[ek] = rng.random() ** 3
            nd = dv + wts[ek]
            if nd < dist.get(u, 1e18):
                dist[u] = nd
                par[u] = v
                heapq.heappush(pq, (nd, u))
    path = []
    v = end
    while v is not None:
        path.append(v)
        v = par[v]
    path = path[::-1]
    st = random_state(d, q, s, rng)
    onpath = set(path)
    for k in range(len(path) - 1):
        st.omega[edge_key(path[k], path[k + 1])] = 1
    # plant diamonds at some turns
    for k in range(1, len(path) - 1):
        a, b, c = path[k - 1], path[k], path[k + 1]
        d1 = tuple(y - x for x, y in zip(a, b)); d2 = tuple(y - x for x, y in zip(b, c))
        if d1 == d2 or rng.random() > plant_rate:
            continue
        qv = tuple(x + z - y for x, y, z in zip(a, b, c))
        if qv in onpath or ninf(qv) > L:
            continue
        st.omega[edge_key(a, qv)] = 1
        st.omega[edge_key(qv, c)] = 1
        # make internal corners b, qv of degree 2
        for vv in (b, qv):
            for e in edges_at(vv, d):
                if e not in (edge_key(a, b), edge_key(b, c), edge_key(a, qv), edge_key(qv, c)):
                    st.omega[e] = 0
        cs = [a, b, c, qv]
        lo = tuple(min(t) for t in zip(*cs))
        axes = [i for i in range(d) if len(set(x[i] for x in cs)) == 2]
        P = (lo, axes[0], axes[1])
        st.alpha[P] = 1 if rng.random() < 0.5 else 0
    return st


if __name__ == '__main__':
    d = int(sys.argv[1]); ns = int(sys.argv[2]); seed = int(sys.argv[3]); mode = sys.argv[4] if len(sys.argv) > 4 else 'random'
    if len(sys.argv) > 5:
        L = int(sys.argv[5])
    assert L >= 2 * RHO + 2
    rng = random.Random(seed)
    stats = dict(step1=0, N=0, B=0, S=0, B_xp_is_T=0, B_xp_on_ridge_or_corner=0, S_path_ends_in_box=0,
                 menger_crosscheck=0, pivotal_edges=0, samples=0)
    pc = 0.5 if d == 2 else 0.2488
    for t in range(ns):
        p = rng.choice([pc - 0.02, pc, pc + 0.02, pc + 0.05])
        s = rng.choice([0.25, 0.5, 1.0])
        if mode == 'path':
            st = path_state(d, rng.choice([0.05, 0.15, 0.3] if d == 2 else [0.02, 0.06, 0.12]), s, rng)
        else:
            st = random_state(d, p, s, rng)
        stats['samples'] += 1
        base = bool(A_L(st))
        for e in plus_pivotal_candidates(st):
            b = st.omega.get(e, 0)
            # (+)-pivotal needs A_L(omega^{e,1}) and not A_L(omega^{e,0}); one of them equals base
            if base and b == 0:
                continue
            if (not base) and b == 1:
                continue
            if base:
                if A_L(st, {e: 0}):
                    continue
            else:
                if not A_L(st, {e: 1}):
                    continue
            assert is_plus_pivotal(st, e)
            stats['pivotal_edges'] += 1
            out = phi(st, e, stats)
            verify(st, e, out)
        if (t + 1) % 10 == 0:
            print(t + 1, stats, flush=True)
    print('FINAL', stats)
    print('ALL CHECKS PASSED')

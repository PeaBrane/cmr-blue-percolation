"""Special Path-1 boxes near S_L (the classification lemma of the special boxes), exact.

(1) Intersection box centred at a non-corner point c of S_L, in local coordinates with S_L the column x = 0
    (c = (0,0), B_L on the side x <= 0). Relevant edges (both ends in B_L, not both in S_L):
        m-c (Gc), a-m, m-b (stubs), a-s-, b-s+ (stubs),   a=(-1,-1), m=(-1,0), b=(-1,1), s-=(0,-1), s+=(0,1).
    The vertices c, s-, s+ lie in S_L and are wired to the sink. Boundary conditions: ALL set partitions of {a,m,b}
    (non-crossing or not: a superset), a block X (joined to the origin), and any set of other blocks joined to the
    sink outside. For each beta, (+)-pivot counts over the other 4 relevant edges.
    Checks: no beta with a Gc pivot count > 0 and T count = 0; max N_Gc/N_T; and the per-record Path-1 inequality
        c_Gc * 120 * u4 * N_Gc + c_T * 120 * l4 * N_T <= 0,  u4 = (1+2eta)^4, l4 = (1-2eta)^4.
    The same enumeration is repeated with an arbitrary subset of the 5 relevant edges' parameter deviations:
    the exponent 4 uses that only the other 4 relevant edges matter.
(2) Corner intersection boxes (centre (L,L) and images): the 4 grid edges at the centre either leave B_L or have both
    ends in S_L, so no grid edge is ever pivotal; checked by listing.
(3) Side boxes whose segment lies in S_L: all 3 segment edges have both ends in S_L; checked by listing.
(4) Census for L in {9, 12, 15, 30}: every side / intersection box is classified as generic (inside B_L, origin
    not interior), origin (o = (1,0) interior), S_L-centred, corner, straddling (segment in S_L), or irrelevant
    (no edge with both ends in B_L and not both in S_L); every cell is inside B_L or irrelevant.
usage: python sl_box.py params/cert.json
"""
import json
import sys
from fractions import Fraction as F
from itertools import product

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from geom import box_edges, cls, int_boxes, norm, side_boxes  # noqa: E402

A, M, Bv, C, SM, SP = "a", "m", "b", "c", "s-", "s+"
EDGES = [(M, C, "Gc"), (A, M, "T"), (M, Bv, "T"), (A, SM, "T"), (Bv, SP, "T")]
SINK = {C, SM, SP}


def set_partitions(s):
    if not s:
        yield []
        return
    first, rest = s[0], s[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


def joined(z, blocks, xb, ysel):
    par = {}

    def f(u):
        par.setdefault(u, u)
        while par[u] != u:
            u = par[u]
        return u

    def un(u, v):
        ru, rv = f(u), f(v)
        if ru != rv:
            par[ru] = rv
    for bl in blocks:
        for v in bl[1:]:
            un(v, bl[0])
    for s in SINK:
        un(s, "SINK")
    for i in ysel:
        un(blocks[i][0], "SINK")
    for k, (u, v, _) in enumerate(EDGES):
        if z >> k & 1:
            un(u, v)
    return f(blocks[xb][0]) == f("SINK")


def sl_centred(cert):
    eta = F(cert["eta"])
    ci = [F(x) for x in cert["c_int"]]
    u4, l4 = (1 + 2 * eta) ** 4, (1 - 2 * eta) ** 4
    worst, arg, nbeta, ninf, nfail = F(0), None, 0, 0, 0
    for blocks in set_partitions([A, M, Bv]):
        for xb in range(len(blocks)):
            others = [i for i in range(len(blocks)) if i != xb]
            for mask in product([0, 1], repeat=len(others)):
                ysel = [others[i] for i in range(len(others)) if mask[i]]
                nbeta += 1
                NG = NT = 0
                for k, (_, _, c) in enumerate(EDGES):
                    for z in range(32):
                        if z >> k & 1:
                            continue
                        if joined(z | 1 << k, blocks, xb, ysel) and not joined(z, blocks, xb, ysel):
                            if c == "Gc":
                                NG += 1
                            else:
                                NT += 1
                if NG and not NT:
                    ninf += 1
                if NT and F(NG, NT) > worst:
                    worst, arg = F(NG, NT), (blocks, xb, ysel, NG, NT)
                lhs = ci[0] * 120 * (u4 if ci[0] > 0 else l4) * NG + ci[3] * 120 * (u4 if ci[3] > 0 else l4) * NT
                nfail += lhs > 0
    print(f"(1) S_L-centred intersection box: {nbeta} boundary conditions, G>0=T: {ninf}, "
          f"max N_Gc/N_T = {worst} at {arg}; Path-1 inequality failures: {nfail}")
    return ninf == 0 and nfail == 0 and worst == F(3, 4)


def relevant(e, L):
    (a, b) = e
    inB = max(abs(a[0]), abs(a[1])) <= L and max(abs(b[0]), abs(b[1])) <= L
    onS = max(abs(a[0]), abs(a[1])) == L and max(abs(b[0]), abs(b[1])) == L
    return inB and not onS


def census(L):
    o = (1, 0)
    R = L // 3 + 2
    counts = {}
    ok = True

    def inside(rect):
        x0, x1, y0, y1 = rect
        return -L <= x0 and x1 <= L and -L <= y0 and y1 <= L

    def interior(rect, v):
        x0, x1, y0, y1 = rect
        return x0 < v[0] < x1 and y0 < v[1] < y1

    for fam, boxes in (("side", side_boxes(R)), ("int", int_boxes(R))):
        for rect, seg in boxes:
            E = [norm(e) for e in box_edges(*rect)]
            rel = [e for e in E if relevant(e, L)]
            if not rel:
                kind = "irrelevant"
            elif inside(rect):
                kind = "origin" if interior(rect, o) else "generic"
                if kind == "generic" and interior(rect, o):
                    ok = False
            else:
                # not inside B_L but has relevant edges
                if fam == "side":
                    kind = "straddling"
                    # the counted grid edges (segment) must all have both ends in S_L
                    if any(relevant(e, L) for e in seg):
                        ok = False
                        print("  side box with a relevant segment edge outside B_L:", rect)
                    if interior(rect, o):
                        ok = False
                else:
                    cx, cy = (rect[0] + rect[1]) // 2, (rect[2] + rect[3]) // 2
                    if abs(cx) == L and abs(cy) == L:
                        kind = "corner"
                        if any(relevant(e, L) for e in seg):
                            ok = False
                            print("  corner box with a relevant grid edge:", rect)
                    elif max(abs(cx), abs(cy)) == L:
                        kind = "S_L-centred"
                        # relevant edges must be exactly the 5 of sl_centred (1 Gc + 4 T) in some orientation
                        if sorted(cls(e) for e in rel) != ["Gc", "T", "T", "T", "T"]:
                            ok = False
                            print("  S_L box with unexpected relevant edges:", rect, rel)
                    else:
                        kind = "other"
                        ok = False
                        print("  unclassified box", rect)
                if o in {v for e in E for v in e}:
                    ok = False
                    print("  special box contains the origin:", rect)
            counts[(fam, kind)] = counts.get((fam, kind), 0) + 1
    # origin box must be unique and equal to [0,3]x[-1,1]
    nor = counts.get(("side", "origin"), 0) + counts.get(("int", "origin"), 0)
    ok &= nor == 1
    # cells
    ccount = {"inside": 0, "irrelevant": 0}
    for kx in range(-R, R + 1):
        for ky in range(-R, R + 1):
            rect = (3 * kx, 3 * kx + 3, 3 * ky, 3 * ky + 3)
            E = [norm(e) for e in box_edges(*rect)]
            if inside(rect):
                ccount["inside"] += 1
                if interior(rect, o):
                    ok = False
            else:
                # star edges (with a candidate-corner endpoint) must be irrelevant
                if any(relevant(e, L) and cls(e) in ("P", "T") for e in E):
                    ok = False
                    print("  cell outside B_L with a relevant star edge:", rect)
                ccount["irrelevant"] += 1
    print(f"(4) census L={L}: boxes {dict(sorted(counts.items()))}; cells {ccount}: {'OK' if ok else 'FAIL'}")
    return ok


def main():
    cert = json.load(open(sys.argv[1]))
    ok = sl_centred(cert)
    for L in (9, 12, 15, 30):
        ok &= census(L)
    print("special boxes:", "VERIFIED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

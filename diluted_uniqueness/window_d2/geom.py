"""Geometry of the spacing-3 sparse diamond model in d = 2 (exact local-ratio certificate).

Written from the definitions of the sparse model and of the two box families.

Global objects on Z^2:
  grid lines       x = 0 mod 3 or y = 0 mod 3
  candidate corner  x != 0 and y != 0 mod 3; candidates P_k = (1,1) + 3k + [0,1]^2
  edge classes      P  (both ends candidate corners)  = edges of candidates
                    T  (exactly one end a candidate corner) = stubs
                    Gc (grid edge with an endpoint in 3Z^2), Gm (other grid edges)
  cells             C_k = 3k + [0,3]^2 (Path 2)
  side family       one box per grid segment (Path 1)
  intersection fam. one box c + [-1,1]^2 per c in 3Z^2 (Path 1)

This script
  (1) checks the structural facts used in the proof (every P-edge/stub/candidate lies in exactly one cell,
      stars are disjoint, classes are invariant under the pattern symmetries);
  (2) computes the family multiplicities and weights and checks the weighted coverage identity
      (every edge of a supported class receives total weight exactly 120 in each family);
  (3) writes the engine spec files cell.spec, side.spec, sideo.spec, int.spec.
usage: python geom.py OUTDIR
"""
import sys
from itertools import product

CLS = {"Gc": 0, "Gm": 1, "P": 2, "T": 3}
LW = 120


def m3(a):
    return a % 3


def is_cand(v):
    return m3(v[0]) != 0 and m3(v[1]) != 0


def is_int(v):
    return m3(v[0]) == 0 and m3(v[1]) == 0


def cls(e):
    a, b = e
    ca, cb = is_cand(a), is_cand(b)
    if ca and cb:
        return "P"
    if ca or cb:
        return "T"
    return "Gc" if (is_int(a) or is_int(b)) else "Gm"


def norm(e):
    a, b = e
    return (a, b) if a < b else (b, a)


def box_vertices(x0, x1, y0, y1):
    return [(x, y) for y in range(y0, y1 + 1) for x in range(x0, x1 + 1)]


def box_edges(x0, x1, y0, y1):
    E = []
    for y in range(y0, y1 + 1):
        for x in range(x0, x1 + 1):
            if x < x1:
                E.append(((x, y), (x + 1, y)))
            if y < y1:
                E.append(((x, y), (x, y + 1)))
    return E


def cyclic_boundary(x0, x1, y0, y1):
    """perimeter lattice points, counterclockwise from the lower-left corner"""
    B = [(x, y0) for x in range(x0, x1 + 1)]
    B += [(x1, y) for y in range(y0 + 1, y1 + 1)]
    B += [(x, y1) for x in range(x1 - 1, x0 - 1, -1)]
    B += [(x0, y) for y in range(y1 - 1, y0, -1)]
    return B


def transpose_edge(e):
    return norm(((e[0][1], e[0][0]), (e[1][1], e[1][0])))


# ---------------- families ----------------
def side_boxes(R):
    """(rect, segment edges) for all grid segments with index in [-R, R]."""
    out = []
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            rect = (3 * i, 3 * i + 3, 3 * j - 1, 3 * j + 1)          # horizontal segment [3i,3i+3] x {3j}
            seg = {norm(((3 * i + k, 3 * j), (3 * i + k + 1, 3 * j))) for k in range(3)}
            out.append((rect, seg))
            rect_t = (3 * j - 1, 3 * j + 1, 3 * i, 3 * i + 3)        # vertical segment {3j} x [3i,3i+3]
            seg_t = {transpose_edge(e) for e in seg}
            out.append((rect_t, seg_t))
    return out


def int_boxes(R):
    out = []
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            c = (3 * i, 3 * j)
            rect = (3 * i - 1, 3 * i + 1, 3 * j - 1, 3 * j + 1)
            seg = {norm(e) for e in box_edges(*rect) if c in e}
            out.append((rect, seg))
    return out


def multiplicity(boxes):
    m = {}
    for rect, _ in boxes:
        for e in box_edges(*rect):
            e = norm(e)
            m[e] = m.get(e, 0) + 1
    return m


def box_weights(rect, seg, mult, supported):
    """weight of every edge of the box: grid edges 120 iff on the box's own segment (and class supported);
    star edges 120 / multiplicity (if class supported)."""
    w = {}
    for e in box_edges(*rect):
        e = norm(e)
        c = cls(e)
        if c not in supported:
            w[e] = 0
        elif c in ("P", "T"):
            assert LW % mult[e] == 0, (e, mult[e])
            w[e] = LW // mult[e]
        else:
            w[e] = LW if e in seg else 0
    return w


def check_coverage(name, boxes, supported, R_test):
    mult = multiplicity(boxes)
    tot = {}
    for rect, seg in boxes:
        for e, we in box_weights(rect, seg, mult, supported).items():
            tot[e] = tot.get(e, 0) + we
    bad = 0
    n = 0
    for x in range(-R_test, R_test + 1):
        for y in range(-R_test, R_test + 1):
            for e in (((x, y), (x + 1, y)), ((x, y), (x, y + 1))):
                e = norm(e)
                want = LW if cls(e) in supported else 0
                n += 1
                if tot.get(e, 0) != want:
                    bad += 1
    print(f"coverage {name}: {n} edges tested in [-{R_test},{R_test}]^2, failures {bad}")
    assert bad == 0
    return mult


def check_structure(R_test=12):
    """every P-edge, stub and candidate lies in exactly one cell's star; stars are pairwise disjoint;
    grid edges lie on cell sides; classes invariant under 3Z^2 translations, x->-x and transpose."""
    star_owner = {}
    for kx in range(-R_test, R_test + 1):
        for ky in range(-R_test, R_test + 1):
            corners = [(3 * kx + 1, 3 * ky + 1), (3 * kx + 2, 3 * ky + 1), (3 * kx + 2, 3 * ky + 2), (3 * kx + 1, 3 * ky + 2)]
            for v in corners:
                for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    e = norm((v, (v[0] + d[0], v[1] + d[1])))
                    if e in star_owner and star_owner[e] != (kx, ky):
                        raise SystemExit(f"star overlap at {e}")
                    star_owner[e] = (kx, ky)
            cell = set(norm(e) for e in box_edges(3 * kx, 3 * kx + 3, 3 * ky, 3 * ky + 3))
            star = {e for e, k in star_owner.items() if k == (kx, ky)}
            assert len(star) == 12 and star <= cell
            assert sum(cls(e) == "P" for e in star) == 4 and sum(cls(e) == "T" for e in star) == 8
            assert all(cls(e) in ("Gc", "Gm") for e in cell - star) and len(cell - star) == 12
    n = 0
    for x in range(-3 * R_test + 3, 3 * R_test - 3):
        for y in range(-3 * R_test + 3, 3 * R_test - 3):
            for e in (((x, y), (x + 1, y)), ((x, y), (x, y + 1))):
                e = norm(e)
                c = cls(e)
                n += 1
                if c in ("P", "T"):
                    assert e in star_owner, e
                else:
                    assert e not in star_owner
                    # grid edge: both endpoints on a common grid line
                    (a, b) = e
                    assert (m3(a[0]) == 0 and a[0] == b[0]) or (m3(a[1]) == 0 and a[1] == b[1]), e
                for g in (lambda v: (v[0] + 3, v[1]), lambda v: (v[0], v[1] - 3), lambda v: (-v[0], v[1]),
                          lambda v: (v[1], v[0])):
                    assert cls(norm((g(e[0]), g(e[1])))) == c
    print(f"structure: stars disjoint, 12 edges each (4 P + 8 T) inside their cell; {n} edges classified; "
          f"classes invariant under 3Z^2, x->-x, transpose: OK")


def write_spec(path, name, verts, edges, weights, boundary, origin=None, cand=None):
    vid = {v: i for i, v in enumerate(verts)}
    with open(path, "w") as f:
        f.write(f"# {name}\n")
        f.write(f"NV {len(verts)}\n")
        for v in verts:
            f.write(f"V {vid[v]} {v[0]} {v[1]}\n")
        f.write(f"NE {len(edges)}\n")
        for e in edges:
            f.write(f"E {vid[e[0]]} {vid[e[1]]} {CLS[cls(e)]} {weights[norm(e)]}\n")
        f.write(f"NB {len(boundary)}\n")
        f.write("B " + " ".join(str(vid[v]) for v in boundary) + "\n")
        f.write(f"ORIGIN {vid[origin] if origin is not None else -1}\n")
        if cand is None:
            f.write("CAND -1\n")
        else:
            f.write("CAND " + " ".join(str(vid[v]) for v in cand) + "\n")
    print(f"wrote {path}: NV={len(verts)} NE={len(edges)} NB={len(boundary)} origin={origin} cand={cand}")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    check_structure()
    sb = side_boxes(8)
    ib = int_boxes(8)
    msb = check_coverage("side family (Gc,Gm,P,T)", sb, ("Gc", "Gm", "P", "T"), 12)
    mib = check_coverage("intersection family (Gc,T)", ib, ("Gc", "T"), 12)

    # cell (Path 2): all weights 1
    rect = (0, 3, 0, 3)
    V = box_vertices(*rect)
    E = [norm(e) for e in box_edges(*rect)]
    Bd = cyclic_boundary(*rect)
    assert len(E) == 24 and len(Bd) == 12
    write_spec(f"{out}/cell.spec", "cell 3k+[0,3]^2 at k=0; candidate {1,2}^2", V, E, {e: 1 for e in E}, Bd,
               cand=[(1, 1), (2, 1), (2, 2), (1, 2)])
    print("  cell classes:", {c: sum(cls(e) == c for e in E) for c in CLS})

    # side box (Path 1), reference = horizontal segment [0,3] x {0}
    rect = (0, 3, -1, 1)
    seg = {norm(((k, 0), (k + 1, 0))) for k in range(3)}
    V = box_vertices(*rect)
    E = [norm(e) for e in box_edges(*rect)]
    Bd = cyclic_boundary(*rect)
    w = box_weights(rect, seg, msb, ("Gc", "Gm", "P", "T"))
    print("  side box weights:", sorted(((e, cls(e), w[e]) for e in E), key=lambda t: t[0]))
    write_spec(f"{out}/side.spec", "side box [0,3]x[-1,1]", V, E, w, Bd)
    write_spec(f"{out}/sideo.spec", "side box [0,3]x[-1,1], origin (1,0) interior", V, E, w, Bd, origin=(1, 0))

    rect = (-1, 1, -1, 1)
    V = box_vertices(*rect)
    E = [norm(e) for e in box_edges(*rect)]
    Bd = cyclic_boundary(*rect)
    w = box_weights(rect, {e for e in E if (0, 0) in e}, mib, ("Gc", "T"))
    print("  int box weights:", sorted(((e, cls(e), w[e]) for e in E), key=lambda t: t[0]))
    write_spec(f"{out}/int.spec", "intersection box [-1,1]^2", V, E, w, Bd)


if __name__ == "__main__":
    main()

"""Shared standard-library helpers for the d = 2 exact local-ratio certificate.

Record conventions (those of the C engine engine/lrm.c):

* A box has NE <= 24 inside edges, NB boundary vertices in cyclic order and, for the origin box, one
  extra position (the origin, always its own outside block). A boundary condition is
  beta = (pi, X, Y): pi a non-crossing partition of the NB boundary positions, written as a restricted
  growth string (RGS), and two distinct blocks X < Y. In the origin box Y is the origin block.
* A record holds nine exact integers
      D, N[0][Gc], N[0][Gm], N[0][P], N[0][T], N[1][Gc], N[1][Gm], N[1][P], N[1][T]:
  D counts the inside configurations in which the candidate plaquette is s-pivotal, and N[a][c] is the
  class-c weighted number of pairs (f, z) with f (+)-pivotal for X <-> Y at candidate activation a.
  The Path-1 boxes have no candidate, so D = 0 and N[1] = 0 there.
* Raw engine dump: 80-byte records (uint32 pid, uint8 X, uint8 Y, uint16 0, uint64 counts[9]),
  little endian, in thread order, with a file DUMP.parts of NP-byte partition labels (NP = NB, plus 1
  for the origin box, whose last label is the origin block).
* Canonical form: the records sorted by (labels, X, Y), serialized column by column as
      labels (NP bytes each) || X (1 byte each) || Y (1 byte each) || counts (72 bytes each).
  Its SHA-256 is the order-independent digest recorded with the certificate. The committed
  certificates/*.canon.xz files are this stream, xz-compressed, so the SHA-256 of the decompressed
  file is the canonical digest.
"""

from fractions import Fraction as F
from itertools import product
from math import comb, gcd
import hashlib
import json
import lzma
import struct

RAW = struct.Struct("<IBBH9Q")
COUNTS = struct.Struct("<9Q")
GC, GM, PP, TT = 0, 1, 2, 3

# (file stem, number of boundary positions NB, positions NP, has origin, number of box edges E)
FAMILIES = {
    "cell_reps": (12, 12, False, 24),
    "side": (10, 10, False, 17),
    "sideo": (10, 11, True, 17),
    "int": (8, 8, False, 12),
}


# ------------------------------------------------------------------ canonical form

def canonical_stream(items, npos):
    """items: iterable of bytes objects labels(npos) + X + Y + counts(72). Returns the canonical stream."""
    rows = sorted(items)
    for a, b in zip(rows, rows[1:]):
        if a[:npos + 2] == b[:npos + 2]:
            raise ValueError("duplicate boundary condition")
    lab = b"".join(r[:npos] for r in rows)
    xs = bytes(r[npos] for r in rows)
    ys = bytes(r[npos + 1] for r in rows)
    cs = b"".join(r[npos + 2:] for r in rows)
    return lab + xs + ys + cs


def raw_items(dump_path, npos):
    """Yield labels + X + Y + counts for every record of a raw engine dump."""
    with open(dump_path + ".parts", "rb") as fh:
        parts = fh.read()
    assert len(parts) % npos == 0
    with open(dump_path, "rb") as fh:
        while True:
            chunk = fh.read(RAW.size * 65536)
            if not chunk:
                break
            assert len(chunk) % RAW.size == 0
            for off in range(0, len(chunk), RAW.size):
                pid, x, y, pad = struct.unpack_from("<IBBH", chunk, off)
                assert pad == 0 and (pid + 1) * npos <= len(parts)
                yield parts[pid * npos:(pid + 1) * npos] + bytes((x, y)) + chunk[off + 8:off + RAW.size]


class Canon:
    """A canonical record stream, held as four byte columns."""

    def __init__(self, data, npos):
        width = npos + 2 + COUNTS.size
        if len(data) % width:
            raise ValueError("stream length is not a multiple of the record width")
        n = len(data) // width
        self.n, self.npos = n, npos
        self.sha256 = hashlib.sha256(data).hexdigest()
        mv = memoryview(data)
        self.labels = mv[:n * npos]
        self.xs = mv[n * npos:n * npos + n]
        self.ys = mv[n * npos + n:n * npos + 2 * n]
        self.cs = mv[n * (npos + 2):]

    @classmethod
    def from_xz(cls, path, npos, nrec=None):
        """Decompress a committed stream. With nrec given, decompress into a buffer of the expected size (lower
        peak memory); a stream of any other length is rejected."""
        with lzma.open(path) as fh:
            if nrec is None:
                return cls(fh.read(), npos)
            buf = bytearray(nrec * (npos + 2 + COUNTS.size))
            view, filled = memoryview(buf), 0
            while filled < len(buf):
                got = fh.readinto(view[filled:])
                if not got:
                    raise ValueError("stream shorter than expected")
                filled += got
            if fh.read(1):
                raise ValueError("stream longer than expected")
            return cls(buf, npos)

    def label(self, i):
        return bytes(self.labels[i * self.npos:(i + 1) * self.npos])

    def key(self, i):
        return self.label(i) + bytes((self.xs[i], self.ys[i]))

    def counts(self, i):
        return COUNTS.unpack_from(self.cs, COUNTS.size * i)

    def rows(self):
        """Yield (labels, X, Y, counts) in canonical order."""
        for i, c in enumerate(COUNTS.iter_unpack(self.cs)):
            yield self.label(i), self.xs[i], self.ys[i], c

    def groups(self):
        """Yield (labels, [(X, Y), ...]) for each label, in canonical order (one group per distinct label if the
        stream is sorted)."""
        cur, pairs = None, []
        for i in range(self.n):
            lab = self.label(i)
            if lab != cur:
                if cur is not None:
                    yield cur, pairs
                cur, pairs = lab, []
            pairs.append((self.xs[i], self.ys[i]))
        if cur is not None:
            yield cur, pairs

    def find(self, key):
        """Index of the record with the given labels + X + Y (binary search on the sorted stream)."""
        lo, hi = 0, self.n
        while lo < hi:
            mid = (lo + hi) // 2
            if self.key(mid) < key:
                lo = mid + 1
            else:
                hi = mid
        if lo == self.n or self.key(lo) != key:
            raise KeyError(key)
        return lo

    def strictly_sorted(self):
        return all(self.key(i) < self.key(i + 1) for i in range(self.n - 1))


# ------------------------------------------------------------------ partitions

def rgs(lab):
    seen, out = {}, []
    for v in lab:
        out.append(seen.setdefault(v, len(seen)))
    return tuple(out)


def noncrossing(lab):
    """Four-index definition: no a < b < c < d with lab[a] = lab[c] != lab[b] = lab[d]."""
    n = len(lab)
    for a in range(n):
        for b in range(a + 1, n):
            if lab[b] == lab[a]:
                continue
            for c in range(b + 1, n):
                if lab[c] != lab[a]:
                    continue
                for d in range(c + 1, n):
                    if lab[d] == lab[b]:
                        return False
    return True


def all_rgs(n):
    out = [()]
    for _ in range(n):
        out = [r + (v,) for r in out for v in range(max(r, default=-1) + 2)]
    return out


def gen_noncrossing(n):
    """Pruned left-to-right RGS generation of the non-crossing partitions of n cyclically ordered points, as
    bytes objects (one label byte per point).

    Putting position i into an existing block B creates a crossing iff another block C has
    first(C) < last(B) < last(C)."""
    out = []
    lab = [0] * n

    def rec(i, first, last):
        if i == n:
            out.append(bytes(lab))
            return
        k = len(first)
        for b in range(k):
            lb = last[b]
            if any(c != b and first[c] < lb < last[c] for c in range(k)):
                continue
            lab[i] = b
            last[b] = i
            rec(i + 1, first, last)
            last[b] = lb
        lab[i] = k
        first.append(i)
        last.append(i)
        rec(i + 1, first, last)
        first.pop()
        last.pop()

    rec(0, [], [])
    return out


def validate_generator(nmax):
    """The pruned generator equals the four-index definition on all RGS of length <= nmax."""
    for n in range(1, nmax + 1):
        nc = {bytes(r) for r in all_rgs(n) if noncrossing(r)}
        g = gen_noncrossing(n)
        if len(g) != len(set(g)) or set(g) != nc:
            return False
    return True


def narayana(n, k):
    return comb(n, k) * comb(n, k - 1) // n


def cell_boundary():
    """Lattice points of the boundary of [0,3]^2, counterclockwise from (0,0)."""
    pts = [(x, 0) for x in range(0, 3)] + [(3, y) for y in range(0, 3)]
    pts += [(x, 3) for x in range(3, 0, -1)] + [(0, y) for y in range(3, 0, -1)]
    return pts


def d4_permutations():
    """The 8 symmetries of [0,3]^2 as permutations of the 12 cell-boundary positions."""
    pts = cell_boundary()
    idx = {p: k for k, p in enumerate(pts)}
    maps = [lambda x, y: (x, y), lambda x, y: (3 - y, x), lambda x, y: (3 - x, 3 - y), lambda x, y: (y, 3 - x),
            lambda x, y: (3 - x, y), lambda x, y: (x, 3 - y), lambda x, y: (y, x), lambda x, y: (3 - y, 3 - x)]
    perms = [tuple(idx[g(*p)] for p in pts) for g in maps]
    assert len(set(perms)) == 8
    return perms


def image(lab, perm):
    """RGS of g(pi), and the map from the blocks of pi to the blocks of g(pi)."""
    new = [0] * len(lab)
    for k, v in enumerate(lab):
        new[perm[k]] = v
    img = rgs(new)
    block_map = {}
    for k, v in enumerate(lab):
        block_map[v] = img[perm[k]]
    return img, block_map


# ------------------------------------------------------------------ certificate constants and inequalities

def load_cert(path):
    with open(path) as fh:
        c = json.load(fh)
    return {"aP": F(c["aP"]), "aT": F(c["aT"]), "g": F(c["g"]), "eta": F(c["eta"]),
            "c_side": [F(x) for x in c["c_side"]], "c_int": [F(x) for x in c["c_int"]]}


def check_c0(cert):
    """(C0): a_T = 9/5 a_P, eta = 1/25000 >= max(g, a_P, a_T), c_side + c_int = (g, g, g - a_P, g - a_T),
    c_int supported on {Gc, T}, and the signs used for the zero-grid special boxes."""
    aP, aT, g, eta, cs, ci = (cert[k] for k in ("aP", "aT", "g", "eta", "c_side", "c_int"))
    assert aT == F(9, 5) * aP
    assert eta == F(1, 25000) and eta >= max(g, aP, aT)
    assert ci[GM] == 0 and ci[PP] == 0
    assert [cs[k] + ci[k] for k in range(4)] == [g, g, g - aP, g - aT]
    assert cs[PP] < 0 and cs[TT] < 0 and ci[TT] < 0 and ci[GC] > 0 and cs[GM] > 0


def integer_form(coefs):
    """Scale rational coefficients to integers with a common positive factor."""
    den = 1
    for c in coefs:
        den = den * c.denominator // gcd(den, c.denominator)
    return [int(c * den) for c in coefs]


class RatioTracker:
    """Maximum of pos/neg over rows, compared exactly, with the rows that attain it."""

    def __init__(self):
        self.pos, self.neg, self.rows = 0, 1, []

    def add(self, pos, neg, row):
        if pos * self.neg > self.pos * neg:
            self.pos, self.neg, self.rows = pos, neg, [row]
        elif pos * self.neg == self.pos * neg and pos:
            self.rows.append(row)

    @property
    def value(self):
        return F(self.pos, self.neg)


def c2_forms(cert):
    """Integer forms of (C2) at u = 0 and u = 1/2, in the raw-count convention of the records:
         u = 0  : 4 (a_P n0_P + a_T n0_T) (1+2eta)^23                    <= D (1-2eta)^24
         u = 1/2: 2 (a_P (n0_P + n1_P) + a_T (n0_T + n1_T)) (1+2eta)^23  <= D (1-2eta)^24
    (with N = 2n this is 2 (a_P N^u_P + a_T N^u_T) (1+2eta)^23 <= D (1-2eta)^24)."""
    aP, aT, eta = cert["aP"], cert["aT"], cert["eta"]
    up23, lo24 = (1 + 2 * eta) ** 23, (1 - 2 * eta) ** 24
    return (integer_form([4 * aP * up23, 4 * aT * up23, lo24]),
            integer_form([2 * aP * up23, 2 * aT * up23, lo24]))


def check_c2(rows, cert):
    """(C2) on every cell record. rows yields (labels, X, Y, counts). Returns a summary dict."""
    (A0, B0, C0), (A1, B1, C1) = c2_forms(cert)
    n = viol = d0_star = 0
    worst = [RatioTracker(), RatioTracker()]
    for lab, x, y, c in rows:
        n += 1
        D, n0P, n0T, n1P, n1T = c[0], c[1 + PP], c[1 + TT], c[5 + PP], c[5 + TT]
        if D == 0 and n0P + n0T + n1P + n1T:
            d0_star += 1
        for j, (pos, neg) in enumerate(((A0 * n0P + B0 * n0T, C0 * D),
                                        (A1 * (n0P + n1P) + B1 * (n0T + n1T), C1 * D))):
            if pos > neg:
                viol += 1
            elif pos:
                worst[j].add(pos, neg, (bytes(lab), x, y, c))
    return {"records": n, "violations": viol, "D0_with_star_pivot": d0_star,
            "max_ratio_u0": worst[0].value, "max_ratio_u_half": worst[1].value,
            "argmax_u0": worst[0].rows, "argmax_u_half": worst[1].rows}


def c1_form(cert, family):
    """Integer form of (C1) for a Path-1 box: sum_k c_k rho_k v_k <= 0 with rho_k = (1 +- 2eta)^(E-1)
    according to the sign of c_k; returns (support, integer coefficients)."""
    eta = cert["eta"]
    E = FAMILIES[family][3]
    cc = cert["c_int"] if family == "int" else cert["c_side"]
    up, lo = (1 + 2 * eta) ** (E - 1), (1 - 2 * eta) ** (E - 1)
    support = [k for k in range(4) if cc[k] != 0]
    return support, integer_form([cc[k] * (up if cc[k] > 0 else lo) for k in support])


def check_c1(rows, cert, family):
    support, A = c1_form(cert, family)
    n = viol = extra = 0
    worst = RatioTracker()
    for lab, x, y, c in rows:
        n += 1
        if c[0] or any(c[5:]):
            extra += 1
        v = [c[1 + k] for k in support]
        pos = sum(a * t for a, t in zip(A, v) if a > 0)
        neg = -sum(a * t for a, t in zip(A, v) if a < 0)
        if pos > neg:
            viol += 1
        elif pos:
            worst.add(pos, neg, (bytes(lab), x, y, c))
    return {"records": n, "violations": viol, "nonzero_cell_slots": extra, "max_ratio": worst.value,
            "argmax": worst.rows}


def completeness(canon, nb, origin):
    """Every non-crossing partition of the nb boundary positions occurs, with exactly its admissible block
    pairs: all C(k,2) pairs, or (origin box) the k pairs (X, origin block). Uses sortedness of the stream.
    (Without an origin the one-block partition has no pair and hence no record.)
    Returns (complete, number of non-crossing partitions, number of records expected)."""
    expected = {}
    ncs = gen_noncrossing(nb)
    for lab in ncs:
        k = max(lab) + 1
        if origin:
            expected[bytes(lab) + bytes((k,))] = [(x, k) for x in range(k)]
        elif k > 1:
            expected[bytes(lab)] = [(x, y) for x in range(k) for y in range(x + 1, k)]
    ok, seen = True, 0
    for lab, pairs in canon.groups():
        seen += 1
        ok = ok and expected.get(lab) == pairs
    return ok and seen == len(expected), len(ncs), sum(len(v) for v in expected.values())


# ------------------------------------------------------------------ rigorous logarithms

SERIES_TERMS = 60


def atanh_bounds(z):
    """lo <= atanh(z) <= hi for rational 0 <= z < 1."""
    z = F(z)
    assert 0 <= z < 1
    lo = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(SERIES_TERMS))
    return lo, lo + z ** (2 * SERIES_TERMS + 1) / ((2 * SERIES_TERMS + 1) * (1 - z * z))


def log_bounds_int(n):
    """Enclosure of ln n for an integer n >= 1: e ln 2 + ln(n / 2^e) with 2^e <= n < 2^(e+1)."""
    if n == 1:
        return F(0), F(0)
    e = n.bit_length() - 1
    two_lo, two_hi = (2 * x for x in atanh_bounds(F(1, 3)))
    rest_lo, rest_hi = (2 * x for x in atanh_bounds(F(n - 2 ** e, n + 2 ** e)))
    return e * two_lo + rest_lo, e * two_hi + rest_hi


def log_bounds(x):
    """Enclosure of ln x for a positive rational x."""
    x = F(x)
    a_lo, a_hi = log_bounds_int(x.numerator)
    b_lo, b_hi = log_bounds_int(x.denominator)
    return a_lo - b_hi, a_hi - b_lo


def dec(x, digits, up):
    """x rounded down (up=False) or up (up=True) to `digits` decimals, exactly."""
    q = x.numerator * 10 ** digits
    q = -((-q) // x.denominator) if up else q // x.denominator
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def sci(x, sig, up):
    """Positive rational x in scientific notation with `sig` significant digits, rounded down or up."""
    e = len(str(x.numerator // x.denominator)) - 1 if x >= 1 else -len(str(x.denominator // x.numerator))
    while F(10) ** e > x:
        e -= 1
    while F(10) ** (e + 1) <= x:
        e += 1
    m = x / F(10) ** (e - sig + 1)
    q = -((-m.numerator) // m.denominator) if up else m.numerator // m.denominator
    if q == 10 ** sig:
        q, e = 10 ** (sig - 1), e + 1
    s = str(q)
    return f"{s[0]}.{s[1:]}e{e}"


# ------------------------------------------------------------------ engine spec files

def parse_spec(path):
    spec = {"V": {}, "E": [], "B": [], "origin": -1, "cand": None}
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t or t[0].startswith("#"):
                continue
            if t[0] == "V":
                spec["V"][int(t[1])] = (int(t[2]), int(t[3]))
            elif t[0] == "E":
                spec["E"].append(tuple(int(s) for s in t[1:5]))
            elif t[0] == "B":
                spec["B"] = [int(s) for s in t[1:]]
            elif t[0] == "ORIGIN":
                spec["origin"] = int(t[1])
            elif t[0] == "CAND" and int(t[1]) >= 0:
                spec["cand"] = [int(s) for s in t[1:5]]
            elif t[0] in ("NV", "NE", "NB"):
                spec[t[0]] = int(t[1])
    assert spec["NV"] == len(spec["V"]) and spec["NE"] == len(spec["E"]) and spec["NB"] == len(spec["B"])
    return spec


# ------------------------------------------------------------------ pure-Python engine (boxes without candidate)

class PathOneEngine:
    """Independent re-implementation, in Python, of the record definition for boxes without a candidate
    plaquette (the Path-1 side, origin and intersection boxes).

    Phase 1 enumerates all inside configurations z once. For every closed edge f it records the transition
    (partition of the positions induced by z, partition induced by z + f) with the class weight of f.
    Phase 2, for one outside partition pi, joins every inside partition with pi and adds the weight of each
    transition to the pairs (X, Y) that are joined after but not before. No single-merge shortcut is used."""

    def __init__(self, spec):
        assert spec["cand"] is None
        self.nv, self.edges = spec["NV"], spec["E"]
        self.nb = spec["NB"]
        self.pos = list(spec["B"]) + ([spec["origin"]] if spec["origin"] >= 0 else [])
        self.origin = spec["origin"] >= 0
        self.np = len(self.pos)
        self.parts, self.index = [], {}
        self.trans = {}
        self._phase1()

    def _pid(self, lab):
        i = self.index.get(lab)
        if i is None:
            i = self.index[lab] = len(self.parts)
            self.parts.append(lab)
        return i

    def _phase1(self):
        nv, pos, edges = self.nv, self.pos, self.edges
        ne = len(edges)
        for z in range(1 << ne):
            parent = list(range(nv))

            def find(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a

            for k, (u, v, _, _) in enumerate(edges):
                if z >> k & 1:
                    ru, rv = find(u), find(v)
                    if ru != rv:
                        parent[ru] = rv
            root = [find(a) for a in range(nv)]
            proots = [root[p] for p in pos]
            ib = self._pid(rgs(proots))
            for k, (u, v, c, w) in enumerate(edges):
                if z >> k & 1:
                    continue
                ru, rv = root[u], root[v]
                if ru == rv or ru not in proots or rv not in proots:
                    continue
                ia = self._pid(rgs([rv if r == ru else r for r in proots]))
                t = self.trans.get((ib, ia))
                if t is None:
                    t = self.trans[(ib, ia)] = [0, 0, 0, 0]
                t[c] += w

    def records(self, lab):
        """All records of the outside partition lab (RGS over the NB boundary positions), as a dict
        (X, Y) -> 9 counts, in the engine's convention."""
        pb = list(lab)
        nbt = max(lab) + 1
        if self.origin:
            pb.append(nbt)
            nbt += 1
        joins = []
        for ilab in self.parts:
            parent = list(range(nbt))

            def find(a):
                while parent[a] != a:
                    a = parent[a]
                return a

            first = {}
            for k, l in enumerate(ilab):
                if l in first:
                    a, b = find(pb[k]), find(first[l])
                    if a != b:
                        parent[a] = b
                else:
                    first[l] = pb[k]
            joins.append(tuple(find(b) for b in range(nbt)))
        agg = {}
        for (ib, ia), w in self.trans.items():
            key = (joins[ib], joins[ia])
            t = agg.get(key)
            if t is None:
                agg[key] = list(w)
            else:
                for c in range(4):
                    t[c] += w[c]
        out = {}
        for x in range(nbt):
            for y in range(x + 1, nbt):
                if self.origin and y != nbt - 1:
                    continue
                out[(x, y)] = [0, 0, 0, 0]
        for (jb, ja), w in agg.items():
            for (x, y), t in out.items():
                if ja[x] == ja[y] and jb[x] != jb[y]:
                    for c in range(4):
                        t[c] += w[c]
        return {k: (0, *t, 0, 0, 0, 0) for k, t in out.items()}

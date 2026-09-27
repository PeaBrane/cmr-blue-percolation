"""Macrostep route, code base R0: the macrostep family, states and orbits (independent of the L2 rig code).

A lateral word is a sequence of signed unit vectors in the m = 3 lateral coordinates, of length <= c, in which
no coordinate occurs with both signs ("monotone").  P(w) = y^|w| / Z.  A macrostep is (w, i), i uniform in the
f forward coordinates.  States of the pair chain are z = (A, D); the symmetry group is B_3 (signed permutations
of the lateral coordinates) x S_f (permutations of the forward coordinates).
"""
import itertools
from fractions import Fraction as Fr
import numpy as np
from ind_tables import fup

LETTERS = [(q, s) for q in range(3) for s in (1, -1)]


def words(c):
    out = []
    for l in range(c + 1):
        for seq in itertools.product(LETTERS, repeat=l):
            signs = {}
            ok = True
            for (q, s) in seq:
                if signs.setdefault(q, s) != s:
                    ok = False
                    break
            if ok:
                out.append(seq)
    return out


class Family:
    def __init__(self, f, c, y):
        self.f, self.c, self.y = f, c, Fr(y)
        W = words(c)
        self.W = W
        n = len(W)
        self.LEN = np.array([len(w) for w in W], np.int64)
        self.PTS = np.zeros((n, c + 1, 3), np.int64)
        for a, w in enumerate(W):
            cur = [0, 0, 0]
            for k, (q, s) in enumerate(w):
                cur[q] += s
                self.PTS[a, k + 1] = cur
            for k in range(len(w) + 1, c + 1):
                self.PTS[a, k] = cur          # padding (never read beyond LEN)
        self.END = np.array([self.PTS[a, len(W[a])] for a in range(n)], np.int64)
        self.Z = sum(self.y ** len(w) for w in W)
        self.P = [self.y ** len(w) / self.Z for w in W]
        # exact endpoint law
        el = {}
        for a in range(n):
            e = tuple(int(v) for v in self.END[a])
            el[e] = el.get(e, Fr(0)) + self.P[a]
        assert sum(el.values()) == 1
        self.elaw = el
        self.PAIRMASS = np.array([[fup(self.P[a] * self.P[b] / (f * f)) for b in range(n)] for a in range(n)])
        dAs = sorted({tuple(int(v) for v in (self.END[a] - self.END[b])) for a in range(n) for b in range(n)})
        self.dAs = dAs
        idx = {v: k for k, v in enumerate(dAs)}
        self.DAIDX = np.array([[idx[tuple(int(v) for v in (self.END[a] - self.END[b]))] for b in range(n)]
                               for a in range(n)], np.int64)


def akey(A):
    return tuple(sorted((abs(int(v)) for v in A), reverse=True))


def dkey(D):
    return tuple(sorted((int(v) for v in D), reverse=True))


def a_orbit(Ak):
    return sorted({tuple(s * v for s, v in zip(sg, p)) for p in itertools.permutations(Ak)
                   for sg in itertools.product((1, -1), repeat=3)})


def d_orbit(Dk):
    return sorted(set(itertools.permutations(Dk)))


def a_types(RA):
    return sorted({akey(A) for A in itertools.product(range(-RA, RA + 1), repeat=3) if sum(map(abs, A)) <= RA},
                  key=lambda k: (sum(k), k))


def d_types(f, Dmax1):
    """S_f-types of root-lattice vectors D (sum 0) with |D|_1 <= Dmax1 (even)."""
    R = Dmax1 // 2
    out = set()
    for D in itertools.product(range(-R, R + 1), repeat=f):
        if sum(D) == 0 and sum(map(abs, D)) <= Dmax1:
            out.add(dkey(D))
    return sorted(out, key=lambda k: (sum(map(abs, k)), k))

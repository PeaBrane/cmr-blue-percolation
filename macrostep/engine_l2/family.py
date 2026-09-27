"""The macrostep path family and the rigorous boost tables.

Family (f forward coordinates, m lateral coordinates, d = f + m):  a macrostep is a MONOTONE lateral word
w = (xi_1..xi_l), l <= c (each lateral coordinate is used with one sign only), followed by one forward step e_i,
i uniform in [f].  Word law P(w) = y^{|w|}/Z.  Points of a word: P_k(w) = xi_1 + ... + xi_k, |P_k|_1 = k.
"""
import itertools
import math
import numpy as np
from fractions import Fraction as Fr
from common import fup, add_up, mul_up, b_exact, INF


def monotone_words(m, c):
    out = []
    def rec(w, used):
        out.append(tuple(w))
        if len(w) >= c:
            return
        for j in range(m):
            if used[j] == 0:
                for s in (1, -1):
                    used[j] = s; w.append((j, s)); rec(w, used); w.pop(); used[j] = 0
            else:
                w.append((j, used[j])); rec(w, used); w.pop()
    rec([], [0] * m)
    return out


class Family:
    def __init__(self, d, f, c, y):
        self.d, self.f, self.m, self.c, self.y = d, f, d - f, c, Fr(y)
        m = self.m
        self.W = monotone_words(m, c)
        self.nW = len(self.W)
        self.LN = np.array([len(w) for w in self.W], np.int64)
        self.PT = np.zeros((self.nW, c + 1, m), np.int64)
        for k, w in enumerate(self.W):
            cur = [0] * m
            for i, (j, s) in enumerate(w):
                cur[j] += s
                self.PT[k, i + 1] = cur
        self.EN = np.array([self.PT[k, self.LN[k]] for k in range(self.nW)], np.int64)
        self.Nl = [int((self.LN == l).sum()) for l in range(c + 1)]
        self.Z = sum(self.Nl[l] * self.y ** l for l in range(c + 1))
        self.Pl = [self.y ** l / self.Z for l in range(c + 1)]          # probability of ONE word of length l
        # float upper bounds of pair masses P(w1)P(w2)/f^2, indexed by (l1, l2)
        self.mass_up = np.array([[fup(self.Pl[l1] * self.Pl[l2] / (f * f)) for l2 in range(c + 1)]
                                 for l1 in range(c + 1)])
        # lateral endpoint law (exact) and word-length law
        law = {}
        for k in range(self.nW):
            e = tuple(int(x) for x in self.EN[k])
            law[e] = law.get(e, Fr(0)) + self.Pl[self.LN[k]]
        self.elaw = law
        self.len_law = [self.Nl[l] * self.Pl[l] for l in range(c + 1)]
        assert sum(self.elaw.values()) == 1 and sum(self.len_law) == 1


class Tables:
    """Float upper bounds: bt[r] >= b(r) (Lemma 3.4(d)); Phi2[r, F] >= max_{rho >= r} sum_{k=0}^{c} b(F+|rho-k|);
    psi3[r0, D1] >= sum_{s>=3} Phi2[max(0, r0 - c(s-1)), max(s, D1 - s)]."""
    def __init__(self, d, t, qs, c, Rb=160, rmax=120, Fmax=420, Dmax=80, S=300):
        self.c = c
        be = [None] + [b_exact(d, t, mm) for mm in range(1, Rb + 1)]
        for mm in range(1, Rb):
            assert be[mm + 1] <= qs * be[mm]
        Rtot = Fmax + rmax + c + 10
        bt = np.zeros(Rtot + 1)
        bt[0] = 1e300  # never used: distance 0 pairs are shared vertices (excluded) -- sentinel
        for mm in range(1, Rtot + 1):
            if mm <= Rb:
                bt[mm] = fup(be[mm])
            else:
                bt[mm] = fup(be[Rb] * qs ** (mm - Rb))
        self.bt = bt
        self.qs = qs
        ph = np.zeros((rmax + 1, Fmax + 1))
        for F in range(Fmax + 1):
            vals = []
            for rho in range(0, rmax + c + 1):
                s = 0.0
                for k in range(c + 1):
                    arg = F + abs(rho - k)
                    if arg == 0:
                        s = INF; break
                    s = add_up(s, bt[arg])
                vals.append(s)
            for r in range(rmax + 1):
                ph[r, F] = max(vals[r:max(r, c) + 1])   # nonincreasing for rho >= c
        self.Phi2 = ph
        # tail: sum_{s>S} Phi2[0,s] <= b(S+1) (1 + 2q/(1-q)) / (1-q)
        q = qs
        bS1 = be[S + 1] if S + 1 <= Rb else be[Rb] * q ** (S + 1 - Rb)   # >= b(S+1)
        tail = fup(bS1 * (1 + 2 * q / (1 - q)) / (1 - q))
        ps = np.zeros((rmax + 1, Dmax + 1))
        for r0 in range(rmax + 1):
            for D1 in range(Dmax + 1):
                assert S >= D1 and c * (S - 1) >= r0
                s_ = 0.0
                for s in range(3, S + 1):
                    s_ = add_up(s_, ph[max(0, r0 - c * (s - 1)), max(s, D1 - s)])
                ps[r0, D1] = add_up(s_, tail)
        self.psi3 = ps

"""Independent checks of the Green-function ingredients (not part of the certificate):
 (1) forward u_n(z) from u_table vs brute-force convolution of the step law (e_i - e_j)/f^2 (exact Fractions);
 (2) lateral l_n(z) from lateral_table (upper bounds) vs an FFT float computation of the same laws."""
import sys, itertools
import numpy as np
from fractions import Fraction as Fr
from scipy.signal import fftconvolve
from family import Family
from green_rig import u_table, lateral_table, dkey
f = int(sys.argv[1]); d = int(sys.argv[2])
# (1) forward
N = 7
law = {tuple([0] * f): Fr(1)}
laws = [law]
for n in range(N):
    new = {}
    for z, p in law.items():
        for i in range(f):
            for j in range(f):
                zz = list(z); zz[i] += 1; zz[j] -= 1; zz = tuple(zz)
                new[zz] = new.get(zz, Fr(0)) + p / f ** 2
    law = new; laws.append(law)
zs = [tuple([0] * f), dkey([1] + [0] * (f - 2) + [-1]), dkey([2] + [0] * (f - 2) + [-2]), dkey([1, 1] + [0] * (f - 4) + [-1, -1]),
      dkey([2, -1, -1] + [0] * (f - 3)), dkey([3, 1] + [0] * (f - 4) + [-2, -2])]
U = u_table(f, N, zs)
bad = 0
for z in zs:
    for n in range(N + 1):
        bf = laws[n].get(tuple(z), Fr(0))
        if bf != U[z][n]:
            bad += 1; print("MISMATCH", z, n, bf, U[z][n])
print(f"(1) forward u_n(z): {len(zs)*(N+1)} values checked against brute force, mismatches = {bad}")
# (2) lateral
fam = Family(d, f, 3, Fr(1, 10))
N0 = 10
zreps = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 0), (3, 2, 1), (6, 0, 0), (5, 4, 3)]
L = lateral_table(fam, N0, zreps)
c = fam.c
PX = np.zeros((2 * c + 1,) * 3)
for e, p in fam.elaw.items():
    PX[tuple(x + c for x in e)] += float(p)
K = fftconvolve(PX, PX[::-1, ::-1, ::-1])
cur = np.ones((1, 1, 1))
worst_lo = 0.0; worst_hi = 0.0
for n in range(1, N0 + 1):
    cur = fftconvolve(cur, K)
    h = (cur.shape[0] - 1) // 2
    for z in zreps:
        v = cur[h + z[0], h + z[1], h + z[2]] if max(map(abs, z)) <= h else 0.0
        if v > 1e-12:
            r = L[z][n] / v - 1
            worst_lo = min(worst_lo, r); worst_hi = max(worst_hi, r)
print(f"(2) lateral l_n(z), n<={N0}: rigorous/FFT - 1 in [{worst_lo:.2e}, {worst_hi:.2e}] (FFT abs. error ~1e-16)")

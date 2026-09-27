"""Compute and cache the rigorous lateral table l_n(z) (upper bounds, upward rounding) for all B_3-orbit
representatives z with |z|_1 <= Zmax, n <= N0.  The lateral walk depends only on (m=3, c, y)."""
import sys, pickle, itertools, time
from fractions import Fraction as Fr
from family import Family
from green_rig import lateral_table
c, y, N0, Zmax, out = int(sys.argv[1]), Fr(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
fam = Family(9, 6, c, y)   # the lateral part does not depend on d, f
reps = sorted({tuple(sorted((abs(x) for x in z), reverse=True)) for z in itertools.product(range(-Zmax, Zmax + 1), repeat=3)
               if sum(map(abs, z)) <= Zmax}, key=lambda k: (sum(k), k))
t0 = time.time()
L = lateral_table(fam, N0, reps)
pickle.dump(dict(c=c, y=y, N0=N0, Zmax=Zmax, L=L), open(out, "wb"))
print(f"lateral cache: c={c} y={y} N0={N0} Zmax={Zmax} #reps={len(reps)} l_N0(0)<={L[(0,0,0)][N0]:.4e} [{time.time()-t0:.0f}s]")

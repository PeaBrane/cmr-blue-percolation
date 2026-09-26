"""Re-run the searches that chose (p, g_K, g_h) in ../params.json (optional; NumPy, SciPy).

Nothing here is needed for verification. Floating point only proposes
parameters; the exact checks here and in ../certify.py decide them.

Part-1 procedure (code base A; research notes scripts/certify_all.py, output
results_nearby.txt, local.md Sec. 6), at the nine research points of that
code base (eight certified rows and the uncertified d = 10, t = 3/25). For
each (d, t):
  1. bisect the largest r = p_B / p_A for which a Nelder-Mead point of the
     two-parameter tangent family makes every floating vertex sum <= 0, and
     take p = floor(p_A r * 10^4)/10^4 - 10^-4;
  2. compute N_k and D_k exactly, with tangent weights chosen by a bounded
     floating maximisation;
  3. maximise the floating mean-field density over lines 2h + 2K'S that stay
     0.001 below log(C^S N_k / D_k) in every environment, then round g_K down
     and lower g_h in steps of 10^-6 until the line passes the exact check.
The headline t = 3/25 row keeps this p but uses the exploration-phase line
K' = 0.013, h = -0.09 with e^(2K'), e^(2h) rounded down to five decimals.

New-point procedure for d = 10, t = 27/200 (code base B; research notes
assembly/certify_new_point.py): the same steps with code base B's search
functions (bisection over p directly, a grid for the weights and a
10^-5 grid for K').

The final lines compare the chosen values with ../params.json. The part-1
search takes several minutes.
"""

import sys
from fractions import Fraction as Fr
from math import atanh, ceil, exp, floor, log, tanh
from pathlib import Path
import json

import numpy as np
from scipy.optimize import minimize

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "crosscheck"))
sys.path.insert(0, str(HERE))

from certify_odds import choose_wk, denominator, make, numerator  # noqa: E402  (code base A)
from certify_pB import aligned_witness  # noqa: E402
from certify_pB import certify as certify_pB  # noqa: E402
from explore_pB_opt2 import build  # noqa: E402
from pB_family import Gfam  # noqa: E402
import freeze_params as B  # noqa: E402  (code base B searches and exact functions)

PARAMS = HERE.parent / "params.json"


# ---------------------------------------------------------------- part-1 procedure (certify_all.py)
def pB_params(d, tf):
    m = 2*d; beta = atanh(tf); pre = 1-exp(-4*beta)
    a, Us, PU, NU, SU = build(beta, m)
    def cert(r):
        lam = a*(1-r)/(r*(1-a))
        f = lambda x: np.max(Gfam(beta, lam, Us, PU, NU, SU, x[0], x[1], m))
        return min((minimize(f, [k0, l0], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 4000})
                    for k0 in [0.95, 1.0, 1.05] for l0 in [0.5, 0.7]), key=lambda r: r.fun)
    lo, hi = 0.5, 0.6
    for _ in range(40):
        r = (lo+hi)/2
        if cert(r).fun <= 0: lo = r
        else: hi = r
    b = cert(lo)
    return pre*lo, b.x[0], b.x[1]


def odds_rows(d, t):
    m, w, a, C, ch2 = make(d, t)
    rows = []
    for k in range(m+1):
        wk = Fr(round(choose_wk(d, float(t), k)*10**4), 10**4)
        N, _ = numerator(k, m, a, ch2, wk)
        D, dv = denominator(k, m, a, ch2)
        assert D == dv[(m-k)//2]
        rows.append((k, 2*k-m, N, D, C**(2*k-m)*N/D))
    return rows, m, C


def choose_line(rows, m, margin):
    Ls = np.array([log(r[4]) for r in rows]); Ss = np.array([r[1] for r in rows])
    best = None
    for Kp in np.linspace(0.0, 0.03, 30001):
        h = np.min((Ls-margin-2*Kp*Ss)/2)
        if Kp > abs(h) or m*Kp >= 1: continue
        x = 0.0
        for _ in range(3000): x = np.tanh(abs(h)+m*Kp*x)
        if best is None or (1-x)/2 > best[0]: best = ((1-x)/2, Kp, h)
    return best


def exact_line_check(rows, m, Gam, H):
    ok = all(r[4] >= H*Gam**r[1] for r in rows)
    slack = min(log(r[4])-log(H*Gam**r[1]) for r in rows)
    return ok, slack


def density_check(m, Gam, H, xbar):
    zx = m*xbar; p, q = zx.numerator, zx.denominator
    return (1/H)**q*Gam**p < ((1+xbar)/(1-xbar))**q


def exact_side_conditions(m, Gam, H):
    return (Gam*H <= 1) and (Gam**m < 7)


def part1(d, t):
    tf = float(t)
    pBf, kap, lx = pB_params(d, tf)
    tgt = Fr(floor(pBf*10**4)-1, 10**4)
    okB, _ = certify_pB(d, t, tgt, kap, lx, verbose=False)
    wit = aligned_witness(d, t)
    rows, m, C = odds_rows(d, t)
    rho_f, Kp, h = choose_line(rows, m, margin=0.001)
    Gam = Fr(round(exp(2*Kp)*10**6)-1, 10**6)
    H = Fr(floor(exp(2*h)*10**6)+5, 10**6)
    while not exact_line_check(rows, m, Gam, H)[0]:
        H -= Fr(1, 10**6)
    okL, slack = exact_line_check(rows, m, Gam, H)
    Kpe = log(Gam)/2; he = log(H)/2
    xs = 0.0
    for _ in range(5000): xs = np.tanh(abs(he)+m*Kpe*xs)
    xbar = Fr(ceil(xs*10**4), 10**4)
    okD = density_check(m, Gam, H, xbar) and exact_side_conditions(m, Gam, H)
    rho = (1-xbar)/2
    print(f"d={d} t={t}={tf}: p_B>={float(tgt):.4f} [{'OK' if okB else 'FAIL'}] (witness {float(wit):.5f}, float cert {pBf:.6f});"
          f" Holley Gamma={Gam} H={H} (K'={Kpe:.6f}, h={he:.6f}) [{'OK' if okL else 'FAIL'}, min log-slack {slack:.5f}];"
          f" rho>{float(rho):.4f} [{'OK' if okD else 'FAIL'}]; rho*p_B>{float(rho*tgt):.5f} vs 1/d={1/d:.5f}", flush=True)
    return tgt, Gam, H


# ---------------------------------------------------------------- new-point procedure (certify_new_point.py)
def pB_float(M):
    pre = float(M.pA)
    def cert(p):
        f = lambda x: B.float_max_delta(M, Fr(p), x[0], x[1])
        return min((minimize(f, [k0, l0], method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-14, 'maxiter': 2000})
                    for k0 in (0.98, 1.02) for l0 in (0.5, 0.6)), key=lambda r: r.fun)
    lo, hi = 0.45 * pre, 0.62 * pre
    for _ in range(22):
        mid = (lo + hi) / 2
        if cert(mid).fun <= 0: lo = mid
        else: hi = mid
    return lo


def new_point(d, t):
    M = B.Model(d, t)
    p = Fr(floor(pB_float(M) * 10 ** 4) - 1, 10 ** 4)
    kap, lx = B.search_family(M, p)
    B.exact.bond_floor(M, p, B.tangent_family(M, p, kap, lx))
    rows = []
    for k in range(M.m + 1):
        N = B.exact.numerator_lower(M, k, B.choose_wk(M, k))
        D = B.exact.denominator_upper(M, k)
        rows.append((2 * k - M.m, M.C ** (2 * k - M.m) * N / D))
    Ls = [(S, log(v)) for S, v in rows]
    best = None
    for i in range(3001):
        Kp = i * 1e-5
        h = min((L - 0.001 - 2 * Kp * S) / 2 for S, L in Ls)
        if h >= 0 or Kp > -h or M.m * Kp >= 1: continue
        x = 0.0
        for _ in range(3000): x = tanh(-h + M.m * Kp * x)
        if best is None or (1 - x) / 2 > best[0]: best = ((1 - x) / 2, Kp, h)
    _, Kp, h = best
    gK = Fr(round(exp(2 * Kp) * 10 ** 6) - 1, 10 ** 6)
    gh = Fr(floor(exp(2 * h) * 10 ** 6) + 5, 10 ** 6)
    while not all(v >= gh * gK ** S for S, v in rows):
        gh -= Fr(1, 10 ** 6)
    print(f"d={d} t={t}: p_B >= {p}; Holley line g_K={gK}, g_h={gh}", flush=True)
    return p, gK, gh


def main():
    chosen = {}
    for d, t in [(12, Fr(11, 100)), (12, Fr(23, 200)), (12, Fr(3, 25)), (12, Fr(1, 8)), (12, Fr(13, 100)),
                 (11, Fr(3, 25)), (11, Fr(13, 100)), (10, Fr(3, 25)), (10, Fr(13, 100))]:
        chosen[(d, t)] = part1(d, t)
    chosen[(10, Fr(27, 200))] = new_point(10, Fr(27, 200))
    mismatches = []
    params = json.loads(PARAMS.read_text())
    for row in params["rows"] + params["uncertified_rows"]:
        p, gK, gh = chosen[(row["d"], Fr(row["t"]))]
        frozen = (Fr(row["p"]), Fr(row["g_K"]), Fr(row["g_h"]))
        compare = frozen[:1] if row["label"] == "d12-t3/25" else frozen
        if (p, gK, gh)[:len(compare)] != compare:
            mismatches.append(row["label"])
    if mismatches:
        print("Different from params.json:", mismatches)
        sys.exit(1)
    print("PASS: the searches reproduce (p, g_K, g_h) of every row in params.json, including the "
          "uncertified row (p only for d12-t3/25, whose line is the exploration line)")


if __name__ == "__main__":
    main()

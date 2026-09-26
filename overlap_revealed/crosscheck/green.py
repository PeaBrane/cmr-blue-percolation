"""Rigorous bounds for the collision Green function of two uniform oriented walks.

u_n(d) = P(Z_n = 0), Z_n = S_n - S'_n, S, S' independent uniform oriented walks in Z^d.
u_n = (n!)^2 d^{-2n} [x^n] (sum_k x^k/(k!)^2)^d            (exact, rational)
Tail: u_n <= M_n <= M_{d*floor(n/d)} <= sqrt(d) (2 pi a)^{-(d-1)/2},  a = floor(n/d) >= 1
      (largest multinomial atom; M_n nonincreasing; Robbins' Stirling bounds).
Returns rational enclosures of G_d = sum u_n, S2_d = sum (n+1) u_n, and the u_n themselves.
"""
from fractions import Fraction as Fr
from math import factorial, comb
import mpmath as mp
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

mp.mp.prec = 200


def exact_u(d, N):
    """u_0..u_{N-1} as Fractions."""
    base = [Fr(1, factorial(k) ** 2) for k in range(N)]
    poly = [Fr(1)] + [Fr(0)] * (N - 1)
    for _ in range(d):
        new = [Fr(0)] * N
        for i, a in enumerate(poly):
            if a == 0:
                continue
            for j in range(N - i):
                new[i + j] += a * base[j]
        poly = new
    return [poly[n] * factorial(n) ** 2 / Fr(d) ** (2 * n) for n in range(N)]


iv = mp.iv
iv.prec = 200


def _fr(v):
    """Exact rational value of a raw mpf tuple (sign, man, exp, bc) or an mpf."""
    if isinstance(v, tuple):
        sign, man, exp, _bc = v
    else:
        sign, man, exp, _bc = mp.mpf(v)._mpf_
    val = Fr(man) * Fr(2) ** exp if exp >= 0 else Fr(man, 2 ** (-exp))
    return -val if sign else val


def hurwitz_upper(s, A):
    """Interval containing A^{-s} + A^{1-s}/(s-1) >= sum_{a>=A} a^{-s}  (s>1)."""
    s = iv.mpf(s)
    A = iv.mpf(A)
    return A ** (-s) + A ** (1 - s) / (s - 1)


def tail_bounds(d, A):
    """Rational upper bounds for sum_{n>=dA} u_n and sum_{n>=dA} (n+1) u_n (outward-rounded intervals)."""
    s = iv.mpf(d - 1) / 2
    pref = iv.sqrt(iv.mpf(d)) * (2 * iv.pi) ** (-s)
    t0 = d * pref * hurwitz_upper(s, A)
    t1 = pref * (d * d * hurwitz_upper(s - 1, A) + iv.mpf(d * (d + 1)) / 2 * hurwitz_upper(s, A))
    return _fr(t0._mpi_[1]), _fr(t1._mpi_[1])


def green_bounds(d, A=4):
    N = d * A
    u = exact_u(d, N)
    part = sum(u)
    part2 = sum((n + 1) * u[n] for n in range(N))
    t0q, t1q = tail_bounds(d, A)
    return {
        "u": u,
        "N": N,
        "G_lo": part,
        "G_hi": part + t0q,
        "S2_hi": part2 + t1q,
        "tail0": t0q,
        "tail1": t1q,
    }


if __name__ == "__main__":
    import sys
    ds = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [12]
    for d in ds:
        g = green_bounds(d)
        Fhi = 1 - 1 / g["G_hi"]
        print(f"d={d}: N={g['N']}  G in [{float(g['G_lo']):.10f}, {float(g['G_hi']):.10f}]  "
              f"tail0<={float(g['tail0']):.2e}  S2<={float(g['S2_hi']):.10f}  F_d<={float(Fhi):.10f}  "
              f"u1..u4={[float(x) for x in g['u'][1:5]]}")

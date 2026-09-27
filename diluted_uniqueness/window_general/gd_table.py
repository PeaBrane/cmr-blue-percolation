"""General-d window from the octahedral rule (closed form (E8.1)), exact rational arithmetic.

    q(d)   = 4 d^2 (d-1) (2d^2-3d+4) / 3          (Lemma R8.5, axis-type count; cross-checked by enumeration)
    N_Q(d) = 8 d^2 (2d^2-3d+4) / 3
    c(a,b) = max_{p in {p_-, p_+}} p^-a (1-p)^-b
    K_d    = 2d [ max( c(7,10d-12) / (2(d-2)), c(6,8d-9) ) + c(6,8d-9) ]
    Delta_d = (1 - 2^-(q+1)) / ((q+1) (K_d + N_Q / (2 (1-p_+)))),   g_d = 2 Delta_d
p_+ : p* (d) of Gomes-Pereira-Sanchis, Thm 2.1, bounded above by exact bisection (h is nonincreasing in p), + 1e-6;
p_- : 1/(2d) for d >= 5 (p_c >= 1/(2d-1), Grimmett (1.13)); 1/mu_d - 1e-6 with mu_3 <= 4.7387, mu_4 <= 6.8040
      (Poenitz-Tittmann 2000, Table 2, k = 14).
Also prints beta_*(d) = (1/2) ln(2 p_-(d) / g_d) (lower bound for beta_l(p) inside the window), rounded down.
After the table ([rev1] lines, added after the first audit) it prints the exact rational brackets (p_-, p_+) used
for every row, checks that p_+ has denominator dividing 10^6 (so the six-decimal value shown is exact), and
compares with the finer bracket p_+ = hi_60 + 10^-6 (hi_60 the 60-step bisection upper point; valid as well) and
with beta_* computed from p_c >= 1/(2d-1) for d >= 5.
The logarithm in beta_* is a rigorous rational lower bound (atanh series; the research version used 60-bit mpmath
interval arithmetic); the printed floors are the same.
Usage: python gd_table.py
"""
import itertools
import math
from fractions import Fraction as Fr

PT = {3: Fr(10000, 47387), 4: Fr(10000, 68040)}


def q_NQ(d):
    q = Fr(4 * d * d * (d - 1) * (2 * d * d - 3 * d + 4), 3)
    NQ = Fr(8 * d * d * (2 * d * d - 3 * d + 4), 3)
    assert q.denominator == 1 and NQ.denominator == 1
    return int(q), int(NQ)


def enumerate_q_NQ(d):
    """direct enumeration (cross-check only): Lambda^+ = {|v|_inf <= 2, |v|_1 <= 3}"""
    Lp = set()
    for s in range(4):
        for axes in itertools.combinations(range(d), s):
            for vals in itertools.product((-2, -1, 1, 2), repeat=s):
                if sum(abs(x) for x in vals) <= 3:
                    v = [0] * d
                    for a, x in zip(axes, vals):
                        v[a] = x
                    Lp.add(tuple(v))
    # plaquettes w + {0, e_i, e_j, e_i + e_j} meeting Lp: lower corner w in Lp - {0, e_i, e_j, e_i+e_j}
    cnt = 0
    for i, j in itertools.combinations(range(d), 2):
        corners = [(0, 0), (1, 0), (0, 1), (1, 1)]
        ws = set()
        for v in Lp:
            for a, b in corners:
                w = list(v)
                w[i] -= a
                w[j] -= b
                ws.add(tuple(w))
        cnt += len(ws)
    P = [(0,) * d, (1,) + (0,) * (d - 1), (0, 1) + (0,) * (d - 2), (1, 1) + (0,) * (d - 2)]
    zs = {tuple(a - b for a, b in zip(v, w)) for v in P for w in Lp}
    return len(Lp), cnt, d * len(zs)


def gps_upper(d, fine=False):
    ms = [(d + i) // 3 for i in range(3)]

    def h(p):
        xs = [(1 - p) ** m for m in ms]
        pr = Fr(1)
        for x in xs:
            pr *= 1 - x
        return pr - 2 + sum(xs)
    lo, hi = Fr(0), Fr(1)
    for _ in range(60):
        mid = (lo + hi) / 2
        if h(mid) < 0:
            hi = mid
        else:
            lo = mid
    ub = Fr(math.ceil(hi * 10 ** 6), 10 ** 6)
    assert h(ub) < 0          # h is nonincreasing in p, so p*(d) < ub
    if fine:
        assert h(hi) < 0      # rev1: the dyadic bisection point itself (finer bracket, same argument)
        return hi
    return ub


def floor_sig(x, digits=5):
    e10 = math.floor(math.log10(float(x)))
    scale = Fr(10) ** (digits - 1 - e10)
    m = (x * scale).numerator // (x * scale).denominator
    return f'{m / 10 ** (digits - 1):.{digits - 1}f}e{e10}'


def _atanh_bounds(z, terms=40):
    """lo <= atanh(z) <= hi for rational 0 <= z <= 1/3 (positive series, geometric tail)"""
    lo = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(terms))
    return lo, lo + z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))


def ln_lower(x):
    """rigorous lower bound of ln(x) for rational x >= 1, returned as a float rounded down.
    x = 2^e r with 1 <= r < 2; r is truncated downward to 64 bits, and ln x >= e ln 2 + 2 atanh((r'-1)/(r'+1))."""
    x = Fr(x)
    assert x >= 1
    e = x.numerator.bit_length() - x.denominator.bit_length()
    if Fr(2) ** e > x:
        e -= 1
    r = x / Fr(2) ** e
    assert 1 <= r < 2
    r_low = Fr(math.floor(r * 2 ** 64), 2 ** 64)
    lo = e * 2 * _atanh_bounds(Fr(1, 3))[0] + 2 * _atanh_bounds((r_low - 1) / (r_low + 1))[0]
    f = float(lo)
    return f if Fr(f) <= lo else math.nextafter(f, -math.inf)


def window(d, universal=False, fine=False):
    """universal=True: the closed-form choice p_- = 1/(2d), p_+ = 7/20 (valid for every d >= 3, since
    p_c(Z^d) <= p_c(Z^3) <= 2 sin(pi/18) < 0.3473 and p_c(Z^d) >= 1/(2d-1))"""
    q, NQ = q_NQ(d)
    if universal:
        p_lo, p_hi = Fr(1, 2 * d), Fr(7, 20)
    else:
        p_hi = gps_upper(d, fine) + Fr(1, 10 ** 6)
        p_lo = PT[d] - Fr(1, 10 ** 6) if d in PT else Fr(1, 2 * d)
    c = lambda a, b: max((1 / p) ** a * (1 / (1 - p)) ** b for p in (p_lo, p_hi))
    cF2, cF1 = c(7, 10 * d - 12), c(6, 8 * d - 9)
    K = 2 * d * (max(cF2 / (2 * (d - 2)), cF1) + cF1)
    tail = Fr(1, 2 ** (q + 1)) if q < 400 else Fr(1, 2 ** 400)     # 1 - 2^-(q+1) >= 1 - 2^-400
    Delta = (1 - tail) / ((q + 1) * (K + Fr(NQ) / (2 * (1 - p_hi))))
    return q, NQ, p_lo, p_hi, K, Delta


def main():
    print('cross-check of the hand counts (Lemma R8.5) by enumeration:')
    for d in range(2, 8):
        Lp, q, NQ = enumerate_q_NQ(d)
        qf, NQf = q_NQ(d)
        lpf = Fr((2 * d + 3) * (2 * d * d + 1), 3)
        print(f'  d={d}: |Lambda^+| = {Lp} (formula {lpf}), q = {q} (formula {qf}), N_Q = {NQ} (formula {NQf})')
        assert (Lp, q, NQ) == (lpf, qf, NQf)
    print()
    print('d | p_- | p_+ | q | N_Q | log2 K_d | Delta_d >= | g_d = 2 Delta_d >= | beta_*(d) >= | g_d (p_-=1/(2d), p_+=7/20) >=')
    for d in list(range(3, 13)) + [15, 20, 30, 50, 100]:
        q, NQ, p_lo, p_hi, K, D = window(d)
        g = 2 * D
        bstar = 0.5 * ln_lower(2 * p_lo / g)
        gu = 2 * window(d, universal=True)[5]
        print(f'{d} | {float(p_lo):.6f} | {float(p_hi):.6f} | {q} | {NQ} | {math.log2(K):.3f} | {floor_sig(D)} | '
              f'{floor_sig(g)} | {math.floor(bstar * 1000) / 1000:.3f} | {floor_sig(gu)}')
    # ---- rev1: exact brackets, a finer p_+, and beta_* from p_c >= 1/(2d-1)
    print()
    print('[rev1] exact brackets of the table rows (p_+ = ceil(10^6 hi_60)/10^6 + 10^-6, exactly the value shown):')
    print('d | p_- (exact) | p_+ (exact) | p_+ = 6-decimal value shown | fine p_+ = hi_60 + 10^-6: g_d >= | '
          'beta_*(d) with p_c^- = 1/(2d-1) (d >= 5) >=')
    for d in list(range(3, 13)) + [15, 20, 30, 50, 100]:
        q, NQ, p_lo, p_hi, K, D = window(d)
        assert (p_hi * 10 ** 6).denominator == 1
        shown = Fr(round(float(p_hi) * 10 ** 6), 10 ** 6)
        if d in PT:
            assert p_lo == PT[d] - Fr(1, 10 ** 6)
            lo_txt = {3: '10000/47387 - 1/10^6', 4: '10000/68040 - 1/10^6'}[d]
        else:
            assert p_lo == Fr(1, 2 * d)
            lo_txt = f'1/{2 * d}'
        gf = 2 * window(d, fine=True)[5]
        assert gf >= 2 * D
        pcm = PT[d] if d in PT else Fr(1, 2 * d - 1)
        b2 = 0.5 * ln_lower(2 * pcm / (2 * D))
        print(f'{d} | {lo_txt} | {p_hi.numerator}/{p_hi.denominator} | {shown == p_hi} | {floor_sig(gf)} | '
              f'{math.floor(b2 * 1000) / 1000:.3f}')


if __name__ == '__main__':
    main()

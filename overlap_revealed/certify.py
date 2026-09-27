"""Exact certificate for the overlap-revealed route (explicit dimensions 10, 11, 12).

For every certified row of Table 2 of the manuscript (Section 4.7), this
program reads the frozen inputs in params.json and recomputes, in exact
rational arithmetic, the conditions (C1)-(C6) listed in Section 4.7:

(C1) bond floor p <= p_B: the vertex sums Delta_j <= 0, j = 0, ..., 2d-1, of
     the finite floor certificate, with the frozen tangent points c_m;
(C2) Holley line: C^S Num_k / Den_k >= g_h g_K^S in all 2d+1 environments
     S = 2k - 2d, with the frozen tangent weights w_k of the numerator bound;
(C3) g_h < 1 <= g_K, g_K g_h <= 1 and g_K^(2d) < 2718/1000 < e (so 4dK' < 1);
     also g_K^d < 2718/1000 and d >= 6;
(C4) the mean-field test tanh(|h| + 2dK' mbar) < mbar, which gives the Ising
     plus-density bound rho_- = (1 - mbar)/2;
(C5) lambda_* = d tbar / (1 - 2d tbar) < 1 with tbar = tanh K' = (g_K-1)/(g_K+1);
(C6) eta' < 1 and Score < 1 (dimensions 11, 12), or rho_-^2 >= c_cov/4 and
     Score_c < 1 (dimension 10); and the root-probability bound theta_* > 0.

It then asserts every displayed number that params.json lists for the row
(the table entries and the Section 4 text displays: the Holley slack, the
Green-function bounds and the certified margins, for example), asserts that
the parameter kappa_2d of the manuscript's centered susceptibility bound
exceeds 1 at every row (so that bound does not apply there), and evaluates
the global criterion at the uncertified row d = 10, t = 3/25, where the
Score_c bound exceeds 1.

This is the second of the two implementations described in Section 4.7, with
the floating-point parameter search removed:
tangent points, tangent weights and mbar are read from params.json
(search/freeze_params.py re-derives them). Transcendental quantities are
enclosed by rationals: exp by a Taylor sum with a remainder bound, atanh by
its series with a geometric tail, square roots by integer square roots, pi by
Machin's formula. Upper bounds of the criterion are rounded outward to a
10^-50 grid. Assertions implement every decisive comparison; decimals in the
output are directed summaries of exact rationals. Only the standard library
is used.
"""

from fractions import Fraction as Fr
from math import comb, factorial, isqrt
from pathlib import Path
import hashlib
import json
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

PARAMS = Path(__file__).resolve().parent / "params.json"

R_CUTOFF = 40            # psibar(r) and |A_r| summed exactly for r <= 40, geometric tail beyond
GREEN_BLOCKS = 4         # u_k summed exactly for k < 4d, Stirling-type tail beyond
EXP_TERMS = 30           # degree of the Taylor polynomial bounding exp on [0, 1]
ATANH_TERMS = 25         # atanh series summed for k = 0..25, geometric tail beyond
SQRT_DIGITS = 30         # decimal digits of the integer square-root enclosures
GRID = 10 ** 50          # outward rounding grid for upper bounds
PI_LO = Fr(314159265358979, 10 ** 14)
E_LO = Fr(2718, 1000)


# ------------------------------------------------------------------ exact helpers
def up(x):
    """Smallest multiple of 10^-50 that is >= x."""
    x = Fr(x)
    return Fr(-((-x.numerator * GRID) // x.denominator), GRID)


def sqrt_lo(x):
    x = Fr(x)
    scale = 10 ** SQRT_DIGITS
    r = Fr(isqrt(x.numerator * scale * scale // x.denominator), scale)
    assert r * r <= x
    return r


def sqrt_hi(x):
    r = sqrt_lo(x) + Fr(1, 10 ** SQRT_DIGITS)
    assert r * r >= x
    return r


def exp_hi(x):
    """Upper bound for e^x, 0 <= x <= 1: Taylor sum plus 3 x^(K+1)/(K+1)!."""
    x = Fr(x)
    assert 0 <= x <= 1
    head = sum(x ** k / factorial(k) for k in range(EXP_TERMS + 1))
    return head + 3 * x ** (EXP_TERMS + 1) / factorial(EXP_TERMS + 1)


def atanh_bounds(z):
    """Rational lo <= atanh(z) <= hi for 0 <= z < 1 (positive series, geometric tail)."""
    z = Fr(z)
    assert 0 <= z < 1
    lo = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(ATANH_TERMS + 1))
    return lo, lo + z ** (2 * ATANH_TERMS + 3) / ((2 * ATANH_TERMS + 3) * (1 - z * z))


def atan_bounds(x, pairs=12):
    """Rational lo <= atan(x) <= hi for 0 < x < 1 (alternating series with decreasing terms)."""
    x = Fr(x)
    assert 0 < x < 1
    lo = sum((-1) ** k * x ** (2 * k + 1) / (2 * k + 1) for k in range(2 * pairs))
    return lo, lo + x ** (4 * pairs + 1) / (4 * pairs + 1)


def check_constants():
    """The rational stand-ins for e, pi and e^(-1/100) are valid bounds."""
    assert sum(Fr(1, factorial(k)) for k in range(12)) > E_LO           # Taylor partial sums of e are below e
    pi_lo = 16 * atan_bounds(Fr(1, 5))[0] - 4 * atan_bounds(Fr(1, 239))[1]   # Machin's formula
    assert PI_LO <= pi_lo
    assert exp_hi(Fr(1, 100)) <= Fr(100, 99)                            # e^(-1/100) >= 99/100


def dec_up(x, digits):
    x = Fr(x)
    return _decimal(-((-x.numerator * 10 ** digits) // x.denominator), digits)


def dec_dn(x, digits):
    x = Fr(x)
    return _decimal((x.numerator * 10 ** digits) // x.denominator, digits)


def _decimal(q, digits):
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def log_lower(x):
    """Rational lower bound for log(x), x >= 1, for display: x is first rounded down to 30 decimals,
    then log y = 2 atanh((y-1)/(y+1))."""
    x = Fr(x)
    assert x >= 1
    y = Fr((x.numerator * 10 ** 30) // x.denominator, 10 ** 30)
    return 2 * atanh_bounds((y - 1) / (y + 1))[0]


def binom_dist(n, p):
    """Exact law of the number of + among n independent signs with P(+) = p."""
    dist = [Fr(1)]
    for _ in range(n):
        new = [Fr(0)] * (len(dist) + 1)
        for i, x in enumerate(dist):
            new[i] += x * (1 - p)
            new[i + 1] += x * p
        dist = new
    return dist


def convolve(left, right):
    out = [Fr(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


class Model:
    def __init__(self, d, t):
        self.d, self.m, self.t = d, 2 * d, Fr(t)
        self.w = (1 + self.t) / (1 - self.t)        # e^{2 beta}
        self.a = self.w ** 2 / (1 + self.w ** 2)    # e^{2 beta} / (2 cosh 2 beta)
        self.C = (self.w + 1 / self.w) / 2          # cosh 2 beta
        self.pA = 1 - self.w ** -2                  # 1 - e^{-4 beta}

    def ch2(self, z):
        """2 cosh(beta z) for even z, exactly."""
        assert z % 2 == 0
        return self.w ** (z // 2) + self.w ** (-(z // 2))


# ------------------------------------------------------------------ (C1) bond floor
def tangent_points(M, spec):
    values = {int(U): Fr(v) for U, v in spec["values"].items()}
    if spec["form"] == "rho_times_omega":
        return {U: v * M.ch2(1 + U) for U, v in values.items()}
    assert spec["form"] == "direct", spec["form"]
    return values


def delta_j(M, j, lam, cU):
    """Delta_j: A, B = sums of j and n-j tilted signs, U = A+B (the manuscript's m), z = A-B."""
    n = M.m - 1
    PA, PB = binom_dist(j, M.a), binom_dist(n - j, M.a)
    total = Fr(0)
    for ia, pa in enumerate(PA):
        A = 2 * ia - j
        for ib, pb in enumerate(PB):
            B = 2 * ib - (n - j)
            U, z = A + B, A - B
            kap = M.w ** (-2 * U) - lam
            om = M.ch2(1 + z)                        # omega(z) = 2 cosh(beta (1+z))
            if kap >= 0:
                phi = 1 / om ** 2                    # Jensen branch
            else:
                c = cU[U]
                phi = 3 / c ** 2 - 2 * om / c ** 3  # tangent branch
            total += pa * pb * kap * phi
    return total


def frozen_success(M, f):
    """Opening probability at the aligned frozen cavity with f frozen successes at the target.

    f = 0 is the aligned witness, an upper bound for p_B."""
    n = M.m - 1 - f
    P = binom_dist(n, M.a)
    X = sum(x / M.ch2(1 + f + (2 * i - n)) ** 2 for i, x in enumerate(P))
    Y = sum(x / M.ch2(-1 + f + (2 * i - n)) ** 2 for i, x in enumerate(P))
    return M.pA * M.a * X / (M.a * X + (1 - M.a) * Y)


def susceptibility_kappa(M):
    """kappa_2d = 2d B A^(2d-1), A, B = (F(4 beta) +- F(2 beta))/2, of the centered susceptibility bound.

    For unit couplings F = cosh, so cosh 2beta = C and cosh 4beta = 2C^2 - 1 are rational."""
    cosh4 = 2 * M.C ** 2 - 1
    A, B = (cosh4 + M.C) / 2, (cosh4 - M.C) / 2
    return M.m * B * A ** (M.m - 1)


def bond_floor(M, p, cU):
    r = p / M.pA
    assert 0 < r < 1
    lam = M.w ** 2 * (1 - r) / r
    needed = {U for U in range(-(M.m - 1), M.m, 2) if M.w ** (-2 * U) < lam}
    assert set(cU) == needed, (sorted(cU), sorted(needed))
    assert all(c > 0 for c in cU.values())
    deltas = [delta_j(M, j, lam, cU) for j in range(M.m)]
    assert all(delta <= 0 for delta in deltas), max(deltas)
    return deltas


# ------------------------------------------------------------------ (C2) Holley line
def G_k(M, k):
    """Omega_k(z) = E 2cosh(beta (z + F)), F a sum of 2d-k fair signs; z = k mod 2, |z| <= k+2."""
    km = M.m - k
    P = binom_dist(km, Fr(1, 2))
    return {z: sum(x * M.ch2(z + 2 * l - km) for l, x in enumerate(P))
            for z in range(-k - 2, k + 3) if (z - k) % 2 == 0}


def numerator_lower(M, k, wk):
    """Num_k = min over frozen classes j of the tangent-line bound, tangent points mixed with weight w_k."""
    G = G_k(M, k)
    if k == 0:
        return 1 / G[0] ** 2
    c = [(1 - wk) * G[2 * i - k] + wk * (Fr(i, k) * G[2 * i - k - 2] + Fr(k - i, k) * G[2 * i - k + 2])
         for i in range(k + 1)]
    assert all(ci > 0 for ci in c)
    values = []
    for j in range(k + 1):
        PA, PB = binom_dist(j, M.a), binom_dist(k - j, M.a)
        total = Fr(0)
        for ia, pa in enumerate(PA):
            for ib, pb in enumerate(PB):
                A, B = 2 * ia - j, 2 * ib - (k - j)
                ci = c[ia + ib]
                total += pa * pb * (3 / ci ** 2 - 2 * G[A - B] / ci ** 3)
        values.append(total)
    return min(values)


def denominator_upper(M, k):
    """Den_k = max over frozen classes j of Y_k(j) = E[(2 cosh beta Z_j)^-2]; every class is computed.

    The maximum is asserted to sit at j = floor((2d-k)/2), the class named in Lemma 4.11
    of the manuscript ("Denominator")."""
    km = M.m - k
    fair = binom_dist(k, Fr(1, 2))
    values = []
    for j in range(km + 1):
        P = convolve(convolve(binom_dist(j, M.a), binom_dist(km - j, 1 - M.a)), fair)
        values.append(sum(x / M.ch2(2 * i - M.m) ** 2 for i, x in enumerate(P)))
    assert max(values) == values[km // 2], k
    return values[km // 2]


def holley_line(M, gK, gh, weights):
    assert len(weights) == M.m + 1 and all(0 <= wk <= 1 for wk in weights)
    rows = []
    for k in range(M.m + 1):
        N = numerator_lower(M, k, weights[k])
        D = denominator_upper(M, k)
        S = 2 * k - M.m
        lhs, rhs = M.C ** S * N / D, gh * gK ** S
        assert lhs >= rhs, (k, S)
        rows.append((S, N, D, lhs / rhs))
    return rows


# ------------------------------------------------------------------ (C3), (C4) hypotheses and density
def hypotheses(d, gK, gh):
    assert d >= 6
    assert gh < 1 <= gK and gK * gh <= 1          # h < 0 and 0 <= K' <= |h|
    assert gK ** (2 * d) < E_LO                    # e^{4dK'} < e, so 4dK' < 1
    assert gK ** d < E_LO                          # e^{2dK'} < 2718/1000
    K = atanh_bounds((gK - 1) / (gK + 1))          # K' = atanh((g_K-1)/(g_K+1))
    H = atanh_bounds((1 - gh) / (1 + gh))          # |h| = atanh((1-g_h)/(1+g_h))
    return K, H


def root_test(d, K_hi, H_hi, mbar):
    """tanh(|h| + 2dK' mbar) < mbar, i.e. |h| + 2dK' mbar < atanh(mbar), with enclosures."""
    assert 0 < mbar < 1
    assert H_hi + 2 * d * K_hi * mbar < atanh_bounds(mbar)[0]


# ------------------------------------------------------------------ (C5), (C6) oriented second moment
def partitions(k, max_parts, max_part=None):
    if max_part is None:
        max_part = k
    if k == 0:
        yield ()
        return
    if max_parts == 0:
        return
    for first in range(min(k, max_part), 0, -1):
        for rest in partitions(k - first, max_parts - 1, first):
            yield (first,) + rest


def u_exact(d, k):
    """u_k = P(Z_k = 0) = sum over compositions of (k!/(d^k prod k_i!))^2, via partitions."""
    fk, fd = factorial(k), factorial(d)
    total = 0
    for lam in partitions(k, d):
        mult = {}
        for v in lam:
            mult[v] = mult.get(v, 0) + 1
        arrangements = factorial(d - len(lam))
        for c in mult.values():
            arrangements *= factorial(c)
        arrangements = fd // arrangements
        denom = 1
        for v in lam:
            denom *= factorial(v)
        multinomial = fk // denom
        total += arrangements * multinomial * multinomial
    return Fr(total, d ** (2 * k))


def pow_neg_half(x, s2):
    """Upper bound for x^(-s2/2), x >= 1 rational, s2 a positive integer."""
    x = Fr(x)
    value = x ** (-(s2 // 2))
    if s2 % 2:
        value /= sqrt_lo(x)
    return value


def green(d):
    """Enclosures of G_d = sum u_k and G^(2)_d = sum (k+1) u_k (collision Green function)."""
    A = GREEN_BLOCKS
    N = d * A
    u = [u_exact(d, k) for k in range(N)]
    s2 = d - 1                                     # s = (d-1)/2
    prefactor = sqrt_hi(d) / sqrt_lo(2 * PI_LO) ** s2
    sum_s = pow_neg_half(A, s2) + pow_neg_half(A, s2 - 2) / (Fr(s2, 2) - 1)
    sum_s1 = pow_neg_half(A, s2 - 2) + pow_neg_half(A, s2 - 4) / (Fr(s2, 2) - 2)
    tail0 = d * prefactor * sum_s
    tail1 = prefactor * (d * d * sum_s1 + Fr(d * (d + 1), 2) * sum_s)
    G_lo = sum(u)
    return dict(u=u, G_lo=G_lo, G_hi=G_lo + tail0, tail0=tail0,
                S2_hi=sum((k + 1) * u[k] for k in range(N)) + tail1)


def multinomial_max(d, m):
    q, j = divmod(m, d)
    return Fr(factorial(m), factorial(q + 1) ** j * factorial(q) ** (d - j))


def gamma_r(d, t, m):
    return sum(comb(d, j) / (1 - 2 * t * (d - 2 * j)) ** (m + 1) for j in range(d + 1)) / 2 ** d


def n_radius(d, r):
    """|A_r|, the number of z with coordinate sum 0 and |z|_1 = 2r."""
    return sum(comb(d, k) * comb(d - k, l) * comb(r - 1, k - 1) * comb(r - 1, l - 1)
               for k in range(1, d + 1) for l in range(1, d - k + 1))


def criterion(d, p, gK, gh, rho, green_data):
    assert 0 < p <= 1 and 0 < rho < 1
    R = R_CUTOFF
    kappa = (1 - rho) / rho
    t = (gK - 1) / (gK + 1)                        # tbar = tanh K', exactly
    alpha = 2 * d * t
    assert alpha < 1
    qs = d * t / (1 - alpha)                       # lambda_*
    assert qs < 1
    top = 2 * R + 6
    b_exact = [t ** m * multinomial_max(d, m) * gamma_r(d, t, m) for m in range(top + 1)]
    assert all(b_exact[m + 1] <= qs * b_exact[m] for m in range(top))
    b = [up(x) for x in b_exact]

    def btail(S):
        return sum(b[S:]) + b[top] * qs / (1 - qs)

    psi = {r: up(sum(b[2 * r - abs(s)] for s in range(-r, r + 1)) + 2 * btail(r + 1))
           for r in range(1, R + 1)}
    x1 = up(kappa * psi[1])
    assert x1 <= 1
    E1 = up(exp_hi(x1))
    U1 = E1 - 1
    SU = up(sum(n_radius(d, r) * kappa * psi[r] * E1 for r in range(1, R + 1)))
    c = 2 * qs / (1 - qs)
    first = Fr(comb(R + d, d - 1) ** 2) * (2 * R + 3 + c) * b[R] * qs
    ratio = Fr(R + 1 + d, R + 2) ** 2 * (2 * R + 5 + c) / (2 * R + 3 + c) * qs
    assert ratio < 1
    SU = up(SU + kappa * E1 * first / (1 - ratio))
    g = green_data
    eta = up(sum(min(U1, uk * SU) for uk in g["u"]) + SU * g["tail0"])
    assert eta < 1
    Sh2 = up(g["S2_hi"] / g["G_lo"] ** 2 - 1)
    T1 = up(U1 * Sh2)
    F = up(1 - 1 / g["G_hi"])
    term1 = 1 / (d * rho * p)
    term2 = (F - Fr(1, d)) / rho
    term3 = T1 / ((1 - eta) * rho)
    # c_cov = 2 sinh 2K' / (e^{K'} cosh 2M + e^{-K'} cosh 2M')^2, M = |h| + (2d-1)K', M' = (2d-1)K';
    # a lower bound for c_cov may replace it in the shared-edge lemma.
    e2M = gK ** (2 * d - 1) / gh
    e2Mp = gK ** (2 * d - 1)
    cosh2M, cosh2Mp = (e2M + 1 / e2M) / 2, (e2Mp + 1 / e2Mp) / 2
    c_cov = (gK - 1 / gK) / (sqrt_hi(gK) * cosh2M + cosh2Mp / sqrt_lo(gK)) ** 2
    assert rho ** 2 >= c_cov / 4
    rho_c = rho + c_cov / (4 * rho)
    return dict(kappa=kappa, t=t, alpha=alpha, qs=qs, b1=b[1], psi1=psi[1], U1=U1, SU=SU, eta=eta,
                G_lo=g["G_lo"], G_hi=g["G_hi"], S2_hi=g["S2_hi"], F=F, Sh2=Sh2, T1=T1, term1=term1,
                term2=term2, term3=term3, score=term1 + term2 + term3, c_cov=c_cov, rho_c=rho_c,
                score_c=1 / (d * rho_c * p) + term2 + term3)


def theta_star(d, rho, eta, score):
    """e^{-1/100} rho_- (1-eta')(1-Score)/(1-1/d), with e^{-1/100} >= 99/100."""
    assert score < 1
    return Fr(99, 100) * rho * (1 - eta) * (1 - score) / (1 - Fr(1, d))


# ------------------------------------------------------------------ published displays
def check_published(published, values):
    """Each key is '<quantity>_<direction>'. The direction says how the display bounds the exact value:
    upper (value <= shown), lower (value >= shown), loglower (log value >= shown, via value >= e^shown)
    or exact (value == shown)."""
    for key, shown in published.items():
        name, direction = key.rsplit("_", 1)
        assert name in values, f"unknown published quantity {key}"
        value, bound = values[name], Fr(shown)
        if direction == "upper":
            assert value <= bound, (key, shown)
        elif direction == "lower":
            assert value >= bound, (key, shown)
        elif direction == "loglower":
            assert value >= exp_hi(bound), (key, shown)
        else:
            assert direction == "exact", key
            assert value == bound, (key, shown)


def certify_row(row, green_cache):
    d, t = row["d"], Fr(row["t"])
    p, gK, gh, mbar = (Fr(row[key]) for key in ("p", "g_K", "g_h", "mbar"))
    M = Model(d, t)
    cU = tangent_points(M, row["pB_tangent_points"])
    deltas = bond_floor(M, p, cU)
    witness = frozen_success(M, 0)
    assert p < witness
    weights = [Fr(w) for w in row["holley_weights"]]
    line = holley_line(M, gK, gh, weights)
    (K_lo, K_hi), (H_lo, H_hi) = hypotheses(d, gK, gh)
    root_test(d, K_hi, H_hi, mbar)
    rho = (1 - mbar) / 2
    if d not in green_cache:
        green_cache[d] = green(d)
    g = criterion(d, p, gK, gh, rho, green_cache[d])
    form = row["criterion"]
    assert form in ("basic", "refined_c"), form
    used = g["score"] if form == "basic" else g["score_c"]
    theta = theta_star(d, rho, g["eta"], used)
    assert theta > 0
    ratios = [ratio for _, _, _, ratio in line]
    gamma_min = gh * gK ** (-M.m)
    tanh_2beta = M.a * M.pA
    # The centered susceptibility bound needs kappa_2d < 1; the manuscript states kappa_2d > 1 at every row.
    sg_kappa = susceptibility_kappa(M)
    assert sg_kappa > 1
    values = dict(
        p=p, g_K=gK, g_h=gh, margin=1 - used, rho_c_gain=g["rho_c"] - rho,
        rho_minus=rho, eta_prime=g["eta"], score=g["score"], score_c=g["score_c"], theta_star=theta,
        kappa=g["kappa"], alpha=g["alpha"], abar=2 * d * K_hi, term1=g["term1"], term2=g["term2"],
        term3=g["term3"], pB_witness=witness, p_A=M.pA, tanh_2beta=tanh_2beta, closed_floor=M.pA / 2,
        p_over_tanh_2beta=p / tanh_2beta, pB_witness_over_tanh_2beta=witness / tanh_2beta, sg_kappa=sg_kappa,
        score_closed_floor=2 / (d * rho * M.pA) + g["term2"] + g["term3"],
        max_delta=max(deltas), argmax_delta=deltas.index(max(deltas)),
        **{f"frozen_success_{f}": frozen_success(M, f) for f in range(1, 5)},
        holley_ratio_S_max=ratios[-1], holley_ratio_S_min=ratios[0], holley_ratio_interior_min=min(ratios[1:-1]),
        gamma_min=gamma_min / (1 + gamma_min),
        K_prime_lo=K_lo, K_prime_hi=K_hi, abs_h_lo=H_lo, abs_h_hi=H_hi,
        G_d_lo=g["G_lo"], G_d_hi=g["G_hi"], G2_d=g["S2_hi"], F_d_minus_1_over_d=g["F"] - Fr(1, d))
    check_published(row["published"], values)
    slack_S, _, _, slack = min(line, key=lambda entry: entry[3])
    local_vector = deltas + [witness] + [x for _, N, D, _ in line for x in (N, D)]
    bound_vector = [rho, g["eta"], g["score"], g["score_c"], theta]
    report = {
        "label": row["label"], "d": d, "t": row["t"], "p": row["p"], "g_K": row["g_K"],
        "g_h": row["g_h"], "mbar": row["mbar"], "criterion": form,
        "pB_floor": {"vertex_sums": len(deltas), "max_delta_upper": dec_up(max(deltas), 10),
                     "aligned_witness_upper": dec_up(witness, 7),
                     "over_tanh_2beta": [dec_dn(p / tanh_2beta, 4), dec_up(witness / tanh_2beta, 4)]},
        "holley_line": {"environments": len(line), "min_log_slack_lower": dec_dn(log_lower(slack), 7),
                        "at_S": slack_S, "log_slack_at_S_min_lower": dec_dn(log_lower(ratios[0]), 7),
                        "log_slack_at_S_max_lower": dec_dn(log_lower(ratios[-1]), 7)},
        "hypotheses": {"gK_2d_upper": dec_up(gK ** (2 * d), 4), "abar_upper": dec_up(2 * d * K_hi, 6),
                       "alpha_upper": dec_up(g["alpha"], 6), "lambda_star_upper": dec_up(g["qs"], 6),
                       "kappa_upper": dec_up(g["kappa"], 6)},
        "rho_minus_lower": dec_dn(rho, 7),
        "green": {"G_d": [dec_dn(g["G_lo"], 10), dec_up(g["G_hi"], 10)], "F_d_upper": dec_up(g["F"], 10),
                  "G2_d_upper": dec_up(g["S2_hi"], 10), "G2_over_G2_minus_1_upper": dec_up(g["Sh2"], 8)},
        "boost": {"U1_upper": dec_up(g["U1"], 6), "sum_U_upper": dec_up(g["SU"], 4),
                  "eta_prime_upper": dec_up(g["eta"], 6), "T1_upper": dec_up(g["T1"], 8)},
        "terms_upper": [dec_up(g[key], 6) for key in ("term1", "term2", "term3")],
        "score_upper": dec_up(g["score"], 6),
        "c_cov_lower": dec_dn(g["c_cov"], 8), "rho_c_lower": dec_dn(g["rho_c"], 7),
        "score_c_upper": dec_up(g["score_c"], 6),
        "theta_star_lower": dec_dn(theta, 6),
        "susceptibility_kappa_2d": [dec_dn(sg_kappa, 4), dec_up(sg_kappa, 4)],
        "published_checked": sorted(row["published"]),
    }
    return report, local_vector, bound_vector


def evaluate_uncertified(row, green_cache):
    """Global criterion only, at a row that the manuscript reports as not certified."""
    d, p, gK, gh, mbar = row["d"], Fr(row["p"]), Fr(row["g_K"]), Fr(row["g_h"]), Fr(row["mbar"])
    (_, K_hi), (_, H_hi) = hypotheses(d, gK, gh)
    root_test(d, K_hi, H_hi, mbar)
    if d not in green_cache:
        green_cache[d] = green(d)
    g = criterion(d, p, gK, gh, (1 - mbar) / 2, green_cache[d])
    assert g["score"] > 1 and g["score_c"] > 1
    check_published(row["published"], dict(score=g["score"], score_c=g["score_c"]))
    return {"label": row["label"], "d": d, "t": row["t"], "score_upper": dec_up(g["score"], 6),
            "score_c_upper": dec_up(g["score_c"], 6), "published_checked": sorted(row["published"])}


def digest(values):
    return hashlib.sha256("\n".join(str(v) for v in values).encode()).hexdigest()


def main():
    check_constants()
    params = json.loads(PARAMS.read_text())
    green_cache = {}
    rows, local_all, bound_all = [], [], []
    for row in params["rows"]:
        start = time.monotonic()
        report, local_vector, bound_vector = certify_row(row, green_cache)
        rows.append(report)
        local_all.extend(local_vector)
        bound_all.extend(local_vector + bound_vector)
        used = report["score_upper"] if report["criterion"] == "basic" else report["score_c_upper"]
        print(f"PASS {report['label']}: p_B >= {row['p']}, Holley line in {report['holley_line']['environments']}"
              f" environments, rho_- >= {report['rho_minus_lower']}, {report['criterion']} score <= {used},"
              f" theta_* >= {report['theta_star_lower']}, {len(row['published'])} displays"
              f" ({time.monotonic() - start:.1f}s)", flush=True)
    uncertified = [evaluate_uncertified(row, green_cache) for row in params["uncertified_rows"]]
    for entry in uncertified:
        print(f"PASS {entry['label']} is not certified: Score <= {entry['score_upper']} and"
              f" Score_c <= {entry['score_c_upper']}, both above 1", flush=True)
    summary = {"certified_rows": len(rows), "rows": rows, "uncertified_rows": uncertified,
               "local_vector_sha256": digest(local_all), "bound_vector_sha256": digest(bound_all)}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

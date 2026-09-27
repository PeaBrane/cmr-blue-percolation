"""Exact deciders for the sharpened local inputs in the noisy-cavity class (Lemmas N, P and S).

Standard library only. All decisions are exact (fractions.Fraction and Python integers).

Model constants at t = tanh(beta):  w = e^{2 beta} = (1+t)/(1-t),  a = w^2/(1+w^2),  p_A = 1 - w^-2,
C = cosh(2 beta) = (w + 1/w)/2,  m = 2d,  n = m - 1.

Lemma N (noisy-cavity reduction) allows the cavity factor 2cosh(beta_c y) with any beta_c in
[atanh(t b), beta], b = tanh((2d-1) beta); the certificate stores e^{2 beta_c} as a rational wc' and this module
checks (1 + t b)/(1 - t b) <= wc' <= w exactly. Then f(y) = 2cosh(beta_c y) = wc'^{y/2} + wc'^{-y/2} for even y.

Lemma P (pair certificate for the bond floor p). With r = p/p_A and lam = w^2 (1-r)/r, for odd U in 1..n:
alpha_U = w^U - lam w^-U > 0, gamma_U = lam w^U - w^-U > 0 and phi_U(A, B) = alpha_U A^-2 - gamma_U B^-2.
p_x = (f(x-1), f(x+1)) for x in X = {-n, -n+2, ..., n}. The certificate is a table of rationals psi_U(x) with
  (V) theta psi_U(x) + (1-theta) psi_U(x') >= phi_U(theta p_x + (1-theta) p_x') for all x, x' and theta in [0, 1],
      decided as psi_U(x) > phi_U(p_x) at every point and, for every pair x < x', positivity on [0, 1] of the
      degree-5 polynomial L(theta) A(theta)^2 B(theta)^2 - alpha_U B(theta)^2 + gamma_U A(theta)^2 (exact
      Sturm root count on (0, 1) plus the signs at the end points);
  (C) Psi(j) = sum_{U odd} sum_h C(n-j, h) C(j, (n+U)/2 - h) psi_U(U - 2(2h - (n-j))) <= 0 for j = 0..n.

Lemma S (symmetrised coupled odds certificate, Holley line g_K, g_h). Environment k: N+ has k coordinates and
N- has m-k, S = 2k - m, Lambda'_k = g_h g_K^S C^-S. With i+- the plus counts of sigma in N+-:
  lambda+(i+) = a^{i+} (1-a)^{k-i+} 2^{-(m-k)},  lambda-(i-) = a^{i-} (1-a)^{m-k-i-} 2^{-k},
  kappa = lambda+ - Lambda' lambda-,  kappa_s(i+, i-) = (kappa(i+, i-) + kappa(k-i+, m-k-i-))/2,
  Gamma(l+, l-) = sum_sigma kappa_s * [3 c^-2 - 2 f(sigma.rho) c^-3 if kappa_s > 0; f(sigma.rho)^-2 if kappa_s < 0],
where rho has l+ plus coordinates in N+ and l- in N-, and c = c(i+, i-) > 0 are the stored tangent points. The
certificate requires Gamma(l+, l-) >= 0 for every class and every k. It also checks the stored Lambda'_k, the flip
symmetry Gamma(l+, l-) = Gamma(k-l+, m-k-l-), and g_h < 1 <= g_K, g_K g_h <= 1, g_K^d, g_K^{2d} < 2718/1000.

The counts #{sigma : i+, i-, sigma.rho = D} are products of binomial coefficients (hypergeometric block laws);
overlap_revealed/noisy_cavity/crosscheck/ recounts them by brute-force enumeration.
"""
from fractions import Fraction as Fr
from math import comb, lcm

E_LO = Fr(2718, 1000)


class Model:
    def __init__(self, d, t, wc):
        self.d, self.m, self.n = d, 2 * d, 2 * d - 1
        self.t, self.wc = Fr(t), Fr(wc)
        t = self.t
        self.w = (1 + t) / (1 - t)
        self.a = self.w ** 2 / (1 + self.w ** 2)
        self.pA = 1 - 1 / self.w ** 2
        self.C = (self.w + 1 / self.w) / 2
        n1, n2 = (1 + t) ** self.n, (1 - t) ** self.n
        b = (n1 - n2) / (n1 + n2)                          # tanh((2d-1) beta), exactly
        self.wc_min = (1 + t * b) / (1 - t * b)            # e^{2 atanh(t b)}
        self._f = {}

    def f(self, y):
        """2 cosh(beta_c y) for even y."""
        assert y % 2 == 0
        e = abs(y) // 2
        if e not in self._f:
            self._f[e] = self.wc ** e + 1 / self.wc ** e
        return self._f[e]

    def noisy_class_ok(self):
        return self.wc_min <= self.wc <= self.w


# ------------------------------------------------------------------ exact polynomial positivity (Sturm)
def _trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def _mul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                out[i + j] += x * y
    return out


def _rem(a, b):
    a, b = _trim(a), _trim(b)
    db = len(b) - 1
    while len(a) - 1 >= db and any(a):
        coef, sh = a[-1] / b[-1], len(a) - 1 - db
        for i in range(len(b)):
            a[i + sh] -= coef * b[i]
        a.pop()
        a = _trim(a) if a else [Fr(0)]
    return a


def _eval(p, x):
    s = Fr(0)
    for c in reversed(p):
        s = s * x + c
    return s


def roots_in_open_unit_interval(P):
    """Number of distinct real roots of P in (0, 1), by Sturm's theorem (needs P(0) != 0 != P(1))."""
    P = _trim(P)
    seq = [P, _trim([i * P[i] for i in range(1, len(P))] or [Fr(0)])]
    while any(seq[-1]):
        r = _rem(seq[-2], seq[-1])
        if not any(r):
            break
        seq.append([-c for c in r])

    def variations(x):
        vals = [v for v in (_eval(p, x) for p in seq) if v != 0]
        return sum(1 for u, v in zip(vals, vals[1:]) if (u > 0) != (v > 0))
    return variations(Fr(0)) - variations(Fr(1))


def sturm_positive(P):
    return _eval(P, Fr(0)) > 0 and _eval(P, Fr(1)) > 0 and roots_in_open_unit_interval(P) == 0


def _bernstein(P):
    """Bernstein coefficients on [0, 1] of the polynomial with power coefficients P."""
    n = len(P) - 1
    return [sum(Fr(comb(k, i), comb(n, i)) * P[i] for i in range(k + 1)) for k in range(n + 1)]


def _split(b):
    """de Casteljau subdivision at 1/2: Bernstein coefficients on the two halves."""
    left, right, cur = [b[0]], [b[-1]], list(b)
    while len(cur) > 1:
        cur = [(x + y) / 2 for x, y in zip(cur, cur[1:])]
        left.append(cur[0])
        right.append(cur[-1])
    return left, right[::-1]


def positive_on_unit_interval(P, depth=6):
    """P > 0 on [0, 1]. Sufficient test first: all Bernstein coefficients positive on the pieces of a dyadic
    subdivision (P is a convex combination of them). Where that is inconclusive, Sturm's theorem decides."""
    stack = [(_bernstein(_trim(P)), 0)]
    while stack:
        b, lev = stack.pop()
        if min(b) > 0:
            continue
        if b[0] <= 0 or b[-1] <= 0:               # the value of P at an end point of this piece
            return False
        if lev >= depth:
            return sturm_positive(P)
        left, right = _split(b)
        stack += [(left, lev + 1), (right, lev + 1)]
    return True


# ------------------------------------------------------------------ Lemma P
def pair_certificate(d, entry):
    """Decide Lemma P for one temperature. entry: t, p, wc_prime, psi[U][x] (rational strings).
    Returns a dict of exact results; asserts every condition."""
    M = Model(d, entry["t"], entry["wc_prime"])
    n, w = M.n, M.w
    p = Fr(entry["p"])
    psi = {int(U): {int(x): Fr(v) for x, v in row.items()} for U, row in entry["psi"].items()}
    Us, X = list(range(1, n + 1, 2)), list(range(-n, n + 1, 2))
    assert set(psi) == set(Us) and all(set(psi[U]) == set(X) for U in Us), "psi table shape"
    assert M.noisy_class_ok(), "wc' outside [e^{2 atanh(t b)}, e^{2 beta}]"
    tanh2b = M.pA * M.a
    assert tanh2b == 2 * M.t / (1 + M.t ** 2)
    assert M.pA / 2 < p < tanh2b
    r = p / M.pA
    lam = w * w * (1 - r) / r
    al = {U: w ** U - lam / w ** U for U in Us}
    ga = {U: lam * w ** U - 1 / w ** U for U in Us}
    assert all(al[U] > 0 and ga[U] > 0 for U in Us)
    f = M.f
    points = chords = 0
    for U in Us:
        for x in X:
            points += 1
            assert psi[U][x] > al[U] / f(x - 1) ** 2 - ga[U] / f(x + 1) ** 2, ("point", U, x)
        for i in range(len(X)):
            for j in range(i + 1, len(X)):
                x, y = X[i], X[j]
                A = [f(y - 1), f(x - 1) - f(y - 1)]            # theta p_x + (1 - theta) p_y, first coordinate
                B = [f(y + 1), f(x + 1) - f(y + 1)]
                L = [psi[U][y], psi[U][x] - psi[U][y]]
                A2, B2 = _mul(A, A), _mul(B, B)
                P = _mul(_mul(L, A2), B2)
                for k, c in enumerate(B2):
                    P[k] -= al[U] * c
                for k, c in enumerate(A2):
                    P[k] += ga[U] * c
                chords += 1
                assert positive_on_unit_interval(P), ("chord", U, x, y)
    Psi = []
    for j in range(n + 1):
        tot = Fr(0)
        for U in Us:
            P_plus = (n + U) // 2
            for h in range(0, n - j + 1):
                g = P_plus - h
                if 0 <= g <= j:
                    tot += comb(n - j, h) * comb(j, g) * psi[U][U - 2 * (2 * h - (n - j))]
        Psi.append(tot)
    assert max(Psi) <= 0, "class sums"
    return dict(model=M, p=p, points=points, chords=chords, Psi=Psi, noisy_witness=aligned_witness(M, noisy=True),
                noisy_witness_beta_c=aligned_witness(M, noisy=True, wc=M.wc_min),
                plain_witness=aligned_witness(M, noisy=False))


def aligned_witness(M, noisy=True, wc=None):
    """p_A r at the aligned frozen cavity with the cavity factor 2cosh(beta' y): beta' from the certificate's wc'
    (noisy=True, or e^{2 beta'} = wc if given) or beta' = beta (noisy=False, the plain class). It is an exact upper
    bound for every class-uniform floor of that class."""
    n, a = M.n, M.a
    base = M.w if not noisy else (M.wc if wc is None else Fr(wc))
    g = lambda y: base ** (abs(y) // 2) + 1 / base ** (abs(y) // 2)
    X = sum(comb(n, k) * a ** k * (1 - a) ** (n - k) / g(2 * k - n + 1) ** 2 for k in range(n + 1))
    Y = sum(comb(n, k) * a ** k * (1 - a) ** (n - k) / g(2 * k - n - 1) ** 2 for k in range(n + 1))
    return M.pA * a * X / (a * X + (1 - a) * Y)


# ------------------------------------------------------------------ Lemma S
def _block_counts(nb, l):
    """cnt[i][e] = #{s in {+-1}^nb : #plus(s) = i, s.rho = 2e - nb} for rho = (+^l, -^(nb-l))."""
    cnt = [[0] * (nb + 1) for _ in range(nb + 1)]
    for i in range(nb + 1):
        for a in range(max(0, i - (nb - l)), min(i, l) + 1):   # a plus coordinates of s inside rho's plus block
            dot = 4 * a - 2 * l - 2 * i + nb
            cnt[i][(dot + nb) // 2] += comb(l, a) * comb(nb - l, i - a)
    return cnt


def holley_certificate(d, entry):
    """Decide Lemma S in all 2d+1 environments for one temperature. entry: t, wc_prime and
    holley = {gK, gh, envs: [{k, Lam, c: {"i+,i-": value}}]}. Asserts every condition."""
    M = Model(d, entry["t"], entry["wc_prime"])
    m, a, C = M.m, M.a, M.C
    H = entry["holley"]
    gK, gh = Fr(H["gK"]), Fr(H["gh"])
    assert M.noisy_class_ok(), "wc' outside [e^{2 atanh(t b)}, e^{2 beta}]"
    assert gh < 1 <= gK and gK * gh <= 1 and gK ** d < E_LO and gK ** (2 * d) < E_LO
    envs = {e["k"]: e for e in H["envs"]}
    assert sorted(envs) == list(range(m + 1))
    # f(D) = (P^(2e) + Q^(2e)) / (P Q)^e with wc' = P/Q and e = |D|/2: integer numerators over the common
    # denominator (PQ)^(m/2) for f, and over prod_e (P^(2e) + Q^(2e))^2 for f^-2, so that the class sums are
    # integer arithmetic until the last step.
    P_, Q_ = M.wc.numerator, M.wc.denominator
    half = m // 2
    num_e = [P_ ** (2 * e) + Q_ ** (2 * e) for e in range(half + 1)]
    den_f = (P_ * Q_) ** half
    Fint = [num_e[abs(2 * j - m) // 2] * (P_ * Q_) ** (half - abs(2 * j - m) // 2) for j in range(m + 1)]
    den_g = 1
    for x in num_e:
        den_g *= x * x
    Gint = [den_g // num_e[abs(2 * j - m) // 2] ** 2 * (P_ * Q_) ** abs(2 * j - m) for j in range(m + 1)]
    assert all(Fr(Fint[j], den_f) == M.f(2 * j - m) for j in range(m + 1))
    assert all(Fr(Gint[j], den_g) == 1 / M.f(2 * j - m) ** 2 for j in range(m + 1))
    rows = []
    for k in range(m + 1):
        kp, km = k, m - k
        S = kp - km
        Lam = gh * gK ** S / C ** S
        e = envs[k]
        assert Fr(e["Lam"]) == Lam, ("Lambda'", k)
        c = {tuple(map(int, key.split(","))): Fr(v) for key, v in e["c"].items()}
        lp = [a ** i * (1 - a) ** (kp - i) / 2 ** km for i in range(kp + 1)]
        lm = [a ** i * (1 - a) ** (km - i) / 2 ** kp for i in range(km + 1)]
        # summand of Gamma: ks (3 c^-2 - 2 f c^-3) = A_i - B_i f  (ks > 0), or ks f^-2 (ks < 0)
        const, pos, neg = Fr(0), [], []
        for ip in range(kp + 1):
            for im in range(km + 1):
                ks = (lp[ip] - Lam * lm[im] + lp[kp - ip] - Lam * lm[km - im]) / 2
                if ks > 0:
                    cc = c[(ip, im)]
                    assert cc > 0
                    const += 3 * ks / cc ** 2 * comb(kp, ip) * comb(km, im)     # sum over sigma of A_i
                    pos.append((ip, im, 2 * ks / cc ** 3 / den_f))
                elif ks < 0:
                    neg.append((ip, im, ks / den_g))
        # one common denominator for all coefficients of this environment: integer arithmetic below
        den_k = lcm(*[coef.denominator for _, _, coef in pos + neg])
        pos = [(ip, im, -(coef * den_k).numerator) for ip, im, coef in pos]
        neg = [(ip, im, (coef * den_k).numerator) for ip, im, coef in neg]
        G = {}
        for lminus in range(km + 1):
            cM = _block_counts(km, lminus)
            # VF[im][e+] = sum_{e-} cM[im][e-] Fint[e+ + e-], likewise VG
            VF = [[sum(cM[im][em] * Fint[ep + em] for em in range(km + 1) if cM[im][em]) for ep in range(kp + 1)]
                  for im in range(km + 1)]
            VG = [[sum(cM[im][em] * Gint[ep + em] for em in range(km + 1) if cM[im][em]) for ep in range(kp + 1)]
                  for im in range(km + 1)]
            for lplus in range(kp + 1):
                cP = _block_counts(kp, lplus)
                tot = 0
                for ip, im, coef in pos:
                    tot += coef * sum(x * v for x, v in zip(cP[ip], VF[im]) if x)
                for ip, im, coef in neg:
                    tot += coef * sum(x * v for x, v in zip(cP[ip], VG[im]) if x)
                G[(lplus, lminus)] = const + Fr(tot, den_k)
        assert all(G[(x, y)] == G[(kp - x, km - y)] for (x, y) in G), ("flip symmetry", k)
        gmin = min(G.values())
        assert gmin >= 0, ("Gamma", k)
        rows.append(dict(k=k, S=S, classes=len(G), min_gamma=gmin, argmin=min(G, key=G.get)))
    return dict(gK=gK, gh=gh, rows=rows, min_gamma=min(r["min_gamma"] for r in rows))

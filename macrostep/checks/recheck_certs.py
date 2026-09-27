"""Macrostep route: fresh exact re-check of engine certificate pickles.
[It refuses optimized mode (-O) and exits with status 1 when a certificate fails. It reads the engine pickles written by
../engine/ind_certify.py and ../engine/ind_certify_t.py; ../reproduce_engine.py runs it on every freshly computed
pickle. The arithmetic is exact integer arithmetic; NumPy is needed only to unpickle the engine output.]

For each pickle (written by ind_certify.py or ind_certify_t.py) and for BOTH stored certificates
(the minimal-lambda one: w, lam, H; the best-theta one: w_th, lam_th, H_th) this script decides, in exact integer
arithmetic on the dyadic rationals represented by the stored binary64 numbers:
    (i)  (Mbar w)(z) + H etaF etabar(z) <= lam w(z)   for every z in C',
    (ii) max_z (Gnear w)(z) + H etaF <= lam H,         lam < 1,  w > 0,  etaF < lam,
and computes theta_* = (99/100) / (C_fin (1 + Lambda lam H / (1 - lam))) exactly, Lambda = max(max etabar/w, 1/H),
with C_fin = rho_-^{-1} exp(kappa (b(1) + 2 psi3(0,0))) recomputed here from Lemma 4.26(iv) of the first manuscript (exact rationals,
Taylor upper bound for exp).  It shares no code with the certificate producers.

Usage: python recheck_certs.py file1.pkl [file2.pkl ...]
"""
import sys, pickle, time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")
from fractions import Fraction as Fr
from math import comb, factorial


def dy(x):
    """exact (numerator, exponent) with x = n / 2^k for a finite binary64 x."""
    n, den = float(x).as_integer_ratio()
    k = den.bit_length() - 1
    assert den == 1 << k
    return n, k


def scaled(vals):
    """list of integers N_i and a common K with vals[i] = N_i / 2^K exactly."""
    pairs = [dy(v) for v in vals]
    K = max(k for _, k in pairs)
    return [n << (K - k) for n, k in pairs], K


# ------------------------------------------------------------------------ C_fin from Lemma 4.26(iv)
def largest_multinomial(d, m):
    q, j = divmod(m, d)
    return Fr(factorial(m), factorial(q + 1) ** j * factorial(q) ** (d - j))


def b_exact(d, t, m):
    Gam = sum(Fr(comb(d, j)) / (1 - 2 * t * (d - 2 * j)) ** (m + 1) for j in range(d + 1)) / 2 ** d
    return t ** m * largest_multinomial(d, m) * Gam


def cfin_upper(d, t, c, rho, kappa, S=40):
    t = Fr(t)
    q = d * t / (1 - 2 * d * t)
    assert 0 < 2 * d * t < 1 and q < 1
    b = {m: b_exact(d, t, m) for m in range(1, S + 2 * c + 3)}
    for m in range(1, S + 2 * c + 2):
        assert b[m + 1] <= q * b[m]
    g = 1 + 2 * q / (1 - q)
    phi2_0 = lambda F: max(sum(b[F + abs(r - k)] for k in range(c + 1)) for r in range(c + 1))
    psi = sum(phi2_0(s) for s in range(3, S + 1)) + g * b[S + 1] / (1 - q)
    x = Fr(kappa) * (b[1] + 2 * psi)
    assert 0 <= x <= 1
    ex = sum(x ** k / factorial(k) for k in range(30)) + 3 * x ** 30 / factorial(30)
    return ex / Fr(rho), float(b[1]), float(psi)


def dec(n, k):
    """the decimal string of n / 10^k."""
    if k <= 0:
        return str(n * 10 ** (-k))
    t = str(n).rjust(k + 1, "0")
    return t[:-k] + "." + t[-k:]


def down(x, sig=6):
    """decimal string <= the rational x (x > 0), sig significant digits (safe-direction display of a lower bound)."""
    x = Fr(x)
    assert x > 0
    e = 0
    while x * Fr(10) ** e < Fr(10) ** (sig - 1):
        e += 1
    while x * Fr(10) ** e >= Fr(10) ** sig:
        e -= 1
    n = (x * Fr(10) ** e).numerator // (x * Fr(10) ** e).denominator
    return dec(n, e)


def up(x, sig=7):
    """decimal string >= the rational x (x > 0)."""
    x = Fr(x)
    e = 0
    while x * Fr(10) ** e < Fr(10) ** (sig - 1):
        e += 1
    while x * Fr(10) ** e >= Fr(10) ** sig:
        e -= 1
    y = x * Fr(10) ** e
    n = -((-y.numerator) // y.denominator)
    return dec(n, e)


def check(M, eta, Gnear, etaF, w, lam, H):
    n = len(w)
    lamF, HF, eF = Fr(float(lam)), Fr(float(H)), Fr(float(etaF))
    if not (lamF < 1 and eF < lamF and all(float(x) > 0 for x in w)):
        return False, "basic"
    W, Kw = scaled(list(w))
    ok_i = True
    worst_i = None
    for s in range(n):
        R, Kr = scaled(list(M[s]))
        lhs = sum(a * b for a, b in zip(R, W))                      # scale 2^(Kr+Kw)
        lhsF = Fr(lhs, 1 << (Kr + Kw)) + HF * eF * Fr(float(eta[s]))
        rhs = lamF * Fr(float(w[s]))
        slack = float((rhs - lhsF) / rhs)
        worst_i = slack if worst_i is None else min(worst_i, slack)
        if not lhsF <= rhs:
            ok_i = False
    gam = None
    for s in range(n):
        R, Kr = scaled(list(Gnear[s]))
        v = Fr(sum(a * b for a, b in zip(R, W)), 1 << (Kr + Kw))
        gam = v if gam is None else max(gam, v)
    ok_ii = gam + HF * eF <= lamF * HF
    return ok_i and ok_ii, dict(ok_i=ok_i, ok_ii=ok_ii, min_rel_slack_i=worst_i,
                                slack_ii=float((lamF * HF - gam - HF * eF) / (lamF * HF)), Gamma_w=float(gam))


def theta(eta, w, lam, H, Cfin):
    lamF, HF = Fr(float(lam)), Fr(float(H))
    Lam = max(max(Fr(float(e)) / Fr(float(x)) for e, x in zip(eta, w)), 1 / HF)
    bound = Cfin * (1 + Lam * lamF * HF / (1 - lamF))
    return Fr(99, 100) / bound, bound, Lam


def main():
    failed = 0
    for path in sys.argv[1:]:
        t0 = time.time()
        o = pickle.load(open(path, "rb"))
        cert = o['cert']
        if not cert.get('ok'):
            print(f"{path}: no certificate stored (rho(Mbar)={cert.get('rhoM')})", flush=True)
            continue
        d, c = o['d'], o['c']
        Cfin, b1, psi = cfin_upper(d, o['t'], c, o['rho'], o['kappa'])
        etaF = o['far']['etaF']
        out = []
        for tag, (w, lam, H) in (("min-lam", (cert['w'], cert['lam'], cert['H'])),
                                 ("best-theta", (cert['w_th'], cert['lam_th'], cert['H_th']))):
            ok, info = check(o['M'], o['eta'], o['Gnear'], etaF, w, lam, H)
            failed += not ok
            th, bound, Lam = theta(o['eta'], w, lam, H, Cfin)
            out.append(f"  {tag:10s}: exact (i),(ii) {ok}; lam <= {up(Fr(float(lam)))}; Gamma_w <= {up(Fr(info['Gamma_w']))}; "
                       f"H = {float(H):.4f}; min rel slack (i) = {info['min_rel_slack_i']:.2e}, (ii) = "
                       f"{info['slack_ii']:.2e}; Lambda = {float(Lam):.4f}; C_fin(1+Lam lam H/(1-lam)) <= {up(bound)}; "
                       f"theta_* >= {down(th, 4)}")
        print(f"{path}: d={d} |C'|={len(o['eta'])} etaF={etaF:.6e} C_fin<={up(Cfin)} (b(1)={b1:.7f}, "
              f"psi3(0,0)={psi:.4e}) [{time.time()-t0:.0f}s]", flush=True)
        for line in out:
            print(line, flush=True)
    if failed:
        print(f"FAIL: {failed} certificate(s) do not pass the exact check", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

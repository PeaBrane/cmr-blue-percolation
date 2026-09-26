"""Theorem 4.5 evaluated at part 1's certified comparison models (local.md Thm C4 / Cor D2 / Sec. 6;
results_headline.txt, results_nearby.txt).  Inputs: d, p_B (rational lower bound), g_K = e^{2K'}, g_h = e^{2h}
(exact rationals).  K' = log(g_K)/2, H = |h| = -log(g_h)/2 are enclosed by rational intervals (mpmath.iv, outward).
Upper endpoints are used for the mean-field root (m* is increasing in K', H), for tbar >= tanh K', and in the
denominator of c_cov; the lower endpoint of K' in the numerator of c_cov."""
from fractions import Fraction as Fr
from criterion import evaluate, iv, ivq, q_up, q_lo
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

def enclose(g):
    """Rational [lo, hi] containing log(g)/2, rounded outward to 1e-12."""
    from math import floor, ceil
    x = iv.log(ivq(Fr(g))) / 2
    lo, hi = q_lo(x), q_up(x)
    S = 10 ** 12
    return Fr(floor(lo * S), S), Fr(ceil(hi * S), S)

def run(label, d, pB, gK, gh):
    Klo, Khi = enclose(gK)
    hlo, hhi = enclose(gh)            # h in [hlo, hhi] (negative)
    Hlo, Hhi = -hhi, -hlo             # |h| in [Hlo, Hhi]
    assert Khi <= Hlo and 2 * d * Khi < 1
    o = evaluate(d, Khi, Hhi, Fr(pB), verbose=False)
    rho = o['rho']
    ch = lambda x: (iv.exp(x) + iv.exp(-x)) / 2
    sh = lambda x: (iv.exp(x) - iv.exp(-x)) / 2
    den = iv.exp(ivq(Khi)) * ch(2 * ivq(Hhi + (2 * d - 1) * Khi)) + iv.exp(-ivq(Klo)) * ch(2 * ivq((2 * d - 1) * Khi))
    c = q_lo(2 * sh(2 * ivq(Klo)) / den ** 2)
    assert rho ** 2 >= c / 4
    rho_c = rho + c / (4 * rho)
    t1c = 1 / (d * rho_c * Fr(pB))
    sc = t1c + o['term2'] + o['term3']
    thc = rho * (1 - o['eta']) * (1 - sc) / (1 - Fr(1, d)) if sc < 1 else None
    f = lambda x, k=6: "None" if x is None else f"{float(x):.{k}f}"
    print(f"{label}: d={d} pB>={pB} K' in [{f(Klo,7)},{f(Khi,7)}] |h| in [{f(Hlo,7)},{f(Hhi,7)}] rho_->={f(rho,7)} "
          f"eta'<={f(o['eta'])} | Score<={f(o['score'])} (theta>={f(o['theta'])}) | rho_c>={f(rho_c)} Score_c<={f(sc)} (theta_c>={f(thc)})")
    return o, sc

if __name__ == "__main__" and "sharp" not in sys.argv:
    run("headline C4", 12, "2161/10000", "102634/100000", "83527/100000")
    rows = [(12, "2023/10000", "255577/250000", "873143/1000000", "t=0.110"),
            (12, "2094/10000", "1024283/1000000", "855173/1000000", "t=0.115"),
            (12, "2161/10000", "513169/500000", "208997/250000", "t=0.120 opt"),
            (12, "2224/10000", "1028469/1000000", "407843/500000", "t=0.125"),
            (12, "2283/10000", "515339/500000", "794291/1000000", "t=0.130"),
            (11, "2180/10000", "513251/500000", "427717/500000", "t=0.120"),
            (11, "2310/10000", "6443/6250", "25553/31250", "t=0.130"),
            (10, "2199/10000", "1026679/1000000", "437233/500000", "t=0.120"),
            (10, "2336/10000", "1031101/1000000", "840843/1000000", "t=0.130")]
    for d, pB, gK, gh, lab in rows:
        run(lab, d, pB, gK, gh)


def run_sharp(label, d, pB, gK, gh):
    from covariance import cov_lower_sharp
    from math import floor
    Klo, Khi = enclose(gK)
    hlo, hhi = enclose(gh)
    Hlo, Hhi = -hhi, -hlo
    o = evaluate(d, Khi, Hhi, Fr(pB), verbose=False)
    rho, mbar = o['rho'], o['mbar']
    c, aux = cov_lower_sharp(d, Klo, Khi, Hhi, mbar, Khi)
    assert rho ** 2 >= c / 4
    rho_c = rho + c / (4 * rho)
    sc = 1 / (d * rho_c * Fr(pB)) + o['term2'] + o['term3']
    th = rho * (1 - o['eta']) * (1 - sc) / (1 - Fr(1, d)) if sc < 1 else None
    f = lambda x, k=6: "None" if x is None else f"{float(x):.{k}f}"
    print(f"{label}: d={d} c_cov'>={f(c)} (E(a+b)^2<={f(aux['Eplus'])}, E(a-b)^2<={f(aux['Eminus'])}) rho_c'>={f(rho_c)} "
          f"Score_c'<={f(sc)} theta'>={f(th)}")

if __name__ == "__main__" and "sharp" in sys.argv:
    run_sharp("headline C4", 12, "2161/10000", "102634/100000", "83527/100000")
    run_sharp("t=0.125", 12, "2224/10000", "1028469/1000000", "407843/500000")
    run_sharp("t=0.130", 11, "2310/10000", "6443/6250", "25553/31250")
    run_sharp("t=0.130", 10, "2336/10000", "1031101/1000000", "840843/1000000")
    run_sharp("t=0.120", 10, "2199/10000", "1026679/1000000", "437233/500000")

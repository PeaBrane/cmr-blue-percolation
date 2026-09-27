"""Exact check of the macrostep interaction-matrix certificates (Theorems Q9, Q8, Q8-sharp, Q7 and the Q7 addendum).

Standard library only. For every instance (d, t) of params.json this program recomputes, in exact rational
arithmetic:

(L)  the local inputs:
     - instance 1 (d = 9 and 8, the local methods of the overlap-revealed route): the bond floor p <= p_B from the
       vertex sums Delta_j <= 0 with the frozen tangent points, and the Holley line in all 2d+1 environments with
       the frozen tangent weights (functions of ../overlap_revealed/certify.py);
     - instances 2 and 3 (d = 8 sharpened and d = 7, the noisy-cavity class): Lemma P (chord conditions and class
       sums) and Lemma S (all classes Gamma >= 0) with the frozen data of local_sharpened.json (functions of
       ../overlap_revealed/noisy_cavity/local_certificates.py);
(H)  the hypotheses of Theorem M: g_h < 1 <= g_K, g_K g_h <= 1, g_K^(2d) < 2718/1000, the mean-field root test
     at the frozen xbar, alpha = 2d tanh K' < 1, q_* = d tanh K'/(1 - alpha) < 1, rho_-^2 >= c_cov/4 and
     rho_c p < rho_- (code base B, engine/indep_global.py, which also gives rho_-, rho_c, kappa and t);
(C)  C_fin >= rho_-^-1 exp(kappa (b(1) + 2 psi_3(0,0))) with b(m) = t^m M_m Gamma_m(t) and a Taylor bound for exp,
     and the torus sizes L_0(n) = 2cn + 1 + ceil(log(100 kappa ((c+1)n+1)^2/(1-alpha))/log(1/alpha)) at
     n = 1, 10, 100, with rational enclosures of the logarithms.

For every stored engine certificate (certificates/<label>.json and .bin.xz) it then
(E)  checks the SHA-256 digest of the payload, that the certificate's inputs and the local constants used by the
     engine equal the ones recomputed in (H), and that its state list is exactly the set of orbit types
     {|D|_1 <= 2 R_D, |A|_1 <= R_A} with R_A >= 2c (so C' contains every state that carries intersections);
(G)  re-decides the criterion of Theorem G for both stored certificates (the minimal-lambda one and the
     theta_*-maximising one), reading every binary64 number as the exact dyadic rational it is:
       (i)  (Mbar w)(z) + H eta_F etabar(z) <= lambda w(z) for every z in C',
       (ii) max_z (Gnear w)(z) + H eta_F <= lambda H,
     with w > 0, H > 0 and 0 < eta_F < lambda < 1; and computes exactly
       Lambda = max(max etabar/w, 1/H),   theta_* = (99/100) / (C_fin (1 + Lambda lambda H/(1 - lambda))).

Finally it asserts every displayed number listed in params.json ("displays") and prints a JSON summary equal to
../expected/macrostep.json. The entries Mbar, etabar, Gnear and eta_F of a certificate are outputs of the engine
(engine/ind_certify.py, engine/ind_certify_t.py); reproduce_engine.py recomputes them from scratch.

    python certify_macrostep.py                         # everything (default run)
    python certify_macrostep.py --runs R0-d9-246,...    # a subset of the stored certificates
    python certify_macrostep.py --cert-dir DIR --certificates-only   # re-decide freshly exported certificates
"""
from array import array
from fractions import Fraction as Fr
from itertools import permutations, product
from math import comb, factorial
from pathlib import Path
import argparse
import hashlib
import json
import lzma
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE / "engine"))
sys.path.insert(0, str(ROOT / "overlap_revealed"))
sys.path.insert(0, str(ROOT / "overlap_revealed" / "noisy_cavity"))
import indep_global as IG          # noqa: E402  code base B: local constants (standard library)
import certify as OR               # noqa: E402  overlap-revealed certificate: instance-1 local inputs
import local_certificates as LC    # noqa: E402  Lemmas P and S in the noisy-cavity class

PARAMS = json.loads((HERE / "params.json").read_text())
C_WORD = PARAMS["engine"]["c"]
PSI_CUT = 40                       # psi_3(0,0): exact terms s <= 40, geometric tail beyond


# ------------------------------------------------------------------ exact helpers
def dec_up(x, digits):
    return OR.dec_up(x, digits)


def dec_dn(x, digits):
    return OR.dec_dn(x, digits)


def sig_dn(x, sig=4):
    """Decimal string <= x > 0 with `sig` significant digits (lower bounds)."""
    x = Fr(x)
    assert x > 0
    e = 0
    while x * Fr(10) ** e < 10 ** (sig - 1):
        e += 1
    while x * Fr(10) ** e >= 10 ** sig:
        e -= 1
    y = x * Fr(10) ** e
    n = y.numerator // y.denominator
    return _sci(n, e)


def sig_up(x, sig=7):
    x = Fr(x)
    assert x > 0
    e = 0
    while x * Fr(10) ** e < 10 ** (sig - 1):
        e += 1
    while x * Fr(10) ** e >= 10 ** sig:
        e -= 1
    y = x * Fr(10) ** e
    n = -((-y.numerator) // y.denominator)
    return _sci(n, e)


def _sci(n, e):
    if e <= 0:
        return str(n * 10 ** (-e))
    s = str(n).rjust(e + 1, "0")
    return s[:-e] + "." + s[-e:]


def log_bounds(x, terms=80):
    """Rational lo <= log(x) <= hi for rational x > 0: x = 2^k y with 1 <= y < 2, log y = 2 atanh((y-1)/(y+1))."""
    x = Fr(x)
    assert x > 0
    k = x.numerator.bit_length() - x.denominator.bit_length()
    y = x / Fr(2) ** k
    while y >= 2:
        y, k = y / 2, k + 1
    while y < 1:
        y, k = y * 2, k - 1

    def atanh(z):
        lo = sum(z ** (2 * i + 1) / (2 * i + 1) for i in range(terms))
        return lo, lo + z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    l2 = atanh(Fr(1, 3))
    ly = atanh((y - 1) / (y + 1))
    lo = 2 * ly[0] + (2 * k * l2[0] if k >= 0 else 2 * k * l2[1])
    hi = 2 * ly[1] + (2 * k * l2[1] if k >= 0 else 2 * k * l2[0])
    return lo, hi


def ceil_ratio_of_logs(a, b):
    """ceil(log a / log b) for rationals a > 1, b > 1, asserted to be determined by the enclosures."""
    la, lb = log_bounds(a), log_bounds(b)
    assert la[0] > 0 and lb[0] > 0
    lo, hi = la[0] / lb[1], la[1] / lb[0]
    c_lo, c_hi = -((-lo.numerator) // lo.denominator), -((-hi.numerator) // hi.denominator)
    assert c_lo == c_hi and lo.denominator != 1, "enclosure straddles an integer"
    return c_lo


# ------------------------------------------------------------------ boost constants and C_fin
def boost_constants(d, t, c):
    """b(1) and an upper bound for psi_3(0,0) = sum_{s>=3} Phi_2(0, s), Phi_2(0, F) = max_{rho in [0,c]}
    sum_{k=0}^{c} b(F + |rho - k|), with b(m) = t^m M_m Gamma_m(t) (global.md Lemma 3.4) and the tail
    sum_{s > S} Phi_2(0, s) <= g b(S+1)/(1-q), g = 1 + 2q/(1-q), q = d t/(1 - 2 d t) (macrostep write-up, Lemma G.4(d))."""
    t = Fr(t)
    q = d * t / (1 - 2 * d * t)
    assert 0 < 2 * d * t < 1 and q < 1
    top = PSI_CUT + 2 * c + 2
    b = {m: t ** m * IG.Mmax_multinomial(d, m) * IG.Gam(d, t, m) for m in range(1, top + 1)}
    assert all(b[m + 1] <= q * b[m] for m in range(1, top))
    g = 1 + 2 * q / (1 - q)
    phi2_0 = lambda F: max(sum(b[F + abs(r - k)] for k in range(c + 1)) for r in range(c + 1))
    psi = sum(phi2_0(s) for s in range(3, PSI_CUT + 1)) + g * b[PSI_CUT + 1] / (1 - q)
    return b[1], psi


def cfin_upper(rho, kappa, b1, psi):
    x = Fr(kappa) * (b1 + 2 * psi)
    assert 0 <= x <= 1
    ex = sum(x ** k / factorial(k) for k in range(30)) + 3 * x ** 30 / factorial(30)
    return ex / Fr(rho)


# ------------------------------------------------------------------ instances
def check_instance(key, inst, sharpened, full=True):
    d, t = inst["d"], Fr(inst["t"])
    p, gK, gh, xbar = (Fr(inst[k]) for k in ("p", "g_K", "g_h", "xbar"))
    o = IG.evaluate(d, p, gK, gh, xbar)        # asserts g_K >= 1 > g_h, g_K g_h <= 1, g_K^(2d) < 2718/1000,
    #                                             d >= 6, the root test, alpha < 1, q_* < 1, rho_-^2 >= c_cov/4
    assert gh < 1 <= gK and gK * gh <= 1 and gK ** (2 * d) < Fr(2718, 1000) and IG.mf_ok(d, gK, gh, xbar)
    rho, rho_c, kappa, tb, alpha, qs = o["rho"], o["rho_c"], o["kappa"], o["t"], o["alpha"], o["qs"]
    assert tb == (gK - 1) / (gK + 1) and alpha == 2 * d * tb and alpha < 1 and qs < 1
    assert rho ** 2 >= o["c_cov"] / 4 and rho_c * p < rho      # rho_e p < rho_- with rho_e = rho_c (Lemma 3.5)
    b1, psi = boost_constants(d, tb, C_WORD)
    cfin = cfin_upper(rho, kappa, b1, psi)
    L0 = {n: 2 * C_WORD * n + 1 + ceil_ratio_of_logs(100 * kappa * ((C_WORD + 1) * n + 1) ** 2 / (1 - alpha), 1 / alpha)
          for n in (1, 10, 100)}
    v = dict(rho_minus=rho, rho_c=rho_c, kappa=kappa, t=tb, alpha=alpha, qs=qs, rho_over_rhoc_p=rho / (rho_c * p),
             rhoc_p=rho_c * p, b1=b1, psi3_00=psi, C_fin=cfin, c_cov=o["c_cov"],
             L0_1=Fr(L0[1]), L0_10=Fr(L0[10]), L0_100=Fr(L0[100]))
    report = {"d": d, "t": inst["t"], "p": inst["p"], "g_K": inst["g_K"], "g_h": inst["g_h"], "xbar": inst["xbar"],
              "local_inputs": inst["local_inputs"], "rho_minus_lower": dec_dn(rho, 7), "rho_c_lower": dec_dn(rho_c, 7),
              "kappa_upper": dec_up(kappa, 6), "tanh_Kprime_upper": dec_up(tb, 7), "alpha_upper": dec_up(alpha, 5),
              "q_star_upper": dec_up(qs, 5), "b1_upper": dec_up(b1, 7), "psi3_00_upper": sig_up(psi, 5),
              "C_fin_upper": dec_up(cfin, 6), "L0": [L0[1], L0[10], L0[100]]}
    if not full:
        return v, report, []
    local_vector = []
    if inst["local_inputs"] == "tangent":
        M = OR.Model(d, t)
        OR.hypotheses(d, gK, gh)
        cU = OR.tangent_points(M, inst["pB_tangent_points"])
        deltas = OR.bond_floor(M, p, cU)                         # asserts Delta_j <= 0 for all j
        witness = OR.frozen_success(M, 0)
        assert p < witness
        line = OR.holley_line(M, gK, gh, [Fr(w) for w in inst["holley_weights"]])   # asserts all 2d+1 lines
        slack = min(r[3] for r in line)
        v.update(max_delta=max(deltas), aligned_witness=witness, min_log_slack=OR.log_lower(slack))
        report["local"] = {"bond_floor": f"p_B >= {inst['p']}", "vertex_sums": len(deltas),
                           "max_delta_upper": "-" + sig_dn(-max(deltas), 4),
                           "aligned_witness_upper": dec_up(witness, 6), "holley_environments": len(line),
                           "min_log_slack_lower": dec_dn(OR.log_lower(slack), 6)}
        local_vector = deltas + [witness] + [x for _, N, D, _ in line for x in (N, D)]
    else:
        E = sharpened[f"d{d}-t{inst['t']}"]
        assert E["d"] == d and Fr(E["t"]) == t and Fr(E["p"]) == p
        assert Fr(E["holley"]["gK"]) == gK and Fr(E["holley"]["gh"]) == gh
        P = LC.pair_certificate(d, E)
        H = LC.holley_certificate(d, E)
        v.update(max_Psi=max(P["Psi"]), chords=Fr(P["chords"]), noisy_witness=P["noisy_witness"],
                 plain_witness=P["plain_witness"], p_minus_plain_witness=p - P["plain_witness"],
                 min_gamma=H["min_gamma"])
        report["local"] = {"pair_certificate": f"p_B >= {inst['p']} (Lemmas N, P)", "wc_prime": E["wc_prime"],
                           "chord_polynomials": P["chords"], "point_conditions": P["points"],
                           "max_Psi_upper": "-" + dec_dn(-max(P["Psi"]), 2),
                           "noisy_witness_upper": dec_up(P["noisy_witness"], 6),
                           "plain_witness_upper": dec_up(P["plain_witness"], 6),
                           "holley_environments": len(H["rows"]), "min_Gamma_lower": sig_dn(H["min_gamma"], 3)}
        local_vector = P["Psi"] + [P["noisy_witness"]] + [r["min_gamma"] for r in H["rows"]]
    return v, report, local_vector


# ------------------------------------------------------------------ engine certificates
def orbit_types(f, RA, RD):
    """The engine's state list: [(A-type, D-type) for D-types for A-types] (engine/ind_family.py)."""
    akey = lambda A: tuple(sorted((abs(x) for x in A), reverse=True))
    dkey = lambda D: tuple(sorted(D, reverse=True))
    Ats = sorted({akey(A) for A in product(range(-RA, RA + 1), repeat=3) if sum(map(abs, A)) <= RA},
                 key=lambda k: (sum(k), k))
    Dts = sorted({dkey(D) for D in product(range(-RD, RD + 1), repeat=f) if sum(D) == 0 and sum(map(abs, D)) <= 2 * RD},
                 key=lambda k: (sum(map(abs, k)), k))
    return [(A, D) for D in Dts for A in Ats]


def load_payload(cert_dir, header):
    raw = lzma.decompress((cert_dir / header["payload"]["file"]).read_bytes())
    assert hashlib.sha256(raw).hexdigest() == header["payload"]["sha256_uncompressed"], "payload digest"
    arr = array("d")
    arr.frombytes(raw)
    if sys.byteorder == "big":
        arr.byteswap()
    out, pos = {}, 0
    for item in header["payload"]["layout"]:
        out[item["name"]] = arr[pos:pos + item["count"]]
        pos += item["count"]
    assert pos == len(arr)
    return out


def scaled(values):
    """Integers N_i and K with values[i] = N_i / 2^K exactly (binary64 numbers are dyadic rationals)."""
    pairs = [x.as_integer_ratio() for x in values]
    K = max(den.bit_length() - 1 for _, den in pairs)
    return [n << (K - (den.bit_length() - 1)) for n, den in pairs], K


def decide(M, Gn, eta, etaF, w, lam, H, n):
    """Exact decision of (i), (ii) and the side conditions; returns Gamma_w, Lambda and the minimal relative slack."""
    assert all(x > 0 for x in w) and 0 < lam < 1 and H > 0 and 0 < etaF < lam
    lamF, HF, eF = Fr(lam), Fr(H), Fr(etaF)
    W, Kw = scaled(w)
    HE = HF * eF
    worst = None
    for s in range(n):
        R, Kr = scaled(M[s * n:(s + 1) * n])
        lhs = Fr(sum(a * b for a, b in zip(R, W)), 1 << (Kr + Kw)) + HE * Fr(eta[s])
        rhs = lamF * Fr(w[s])
        assert lhs <= rhs, ("condition (i)", s)
        slack = (rhs - lhs) / rhs
        worst = slack if worst is None else min(worst, slack)
    gam = None
    for s in range(n):
        R, Kr = scaled(Gn[s * n:(s + 1) * n])
        v = Fr(sum(a * b for a, b in zip(R, W)), 1 << (Kr + Kw))
        gam = v if gam is None else max(gam, v)
    assert gam + HE <= lamF * HF, "condition (ii)"
    Lam = max(max(Fr(e) / Fr(x) for e, x in zip(eta, w)), 1 / HF)
    return gam, Lam, worst


def check_run(run, header, cert_dir, inst, iv):
    d = inst["d"]
    fam = PARAMS["engine"]["family"][str(d)]
    I = header["inputs"]
    assert header["label"] == run["label"] and header["code_base"] == run["code_base"]
    assert (I["d"], I["f"], I["c"], I["N0"], I["R_A"], I["R_D"]) == (d, fam["f"], C_WORD, PARAMS["engine"]["N0"],
                                                                      run["R_A"], run["R_D"])
    assert all(Fr(I[k]) == Fr(inst[k]) for k in ("p", "g_K", "g_h", "xbar")) and Fr(I["y"]) == Fr(fam["y"])
    lc = header["engine_local_constants"]
    assert (Fr(lc["rho_minus"]), Fr(lc["rho_c"]), Fr(lc["kappa"]), Fr(lc["t"])) == \
        (iv["rho_minus"], iv["rho_c"], iv["kappa"], iv["t"]), "engine local constants"
    assert run["R_A"] >= 2 * C_WORD
    states = [(tuple(A), tuple(D)) for A, D in header["states"]]
    assert states == orbit_types(fam["f"], run["R_A"], run["R_D"]), "state list"
    n = header["n_states"]
    assert n == len(states) == run["n_states"]
    arrs = load_payload(cert_dir, header)
    etaF = float.fromhex(header["etaF"])
    out = {"etaF": Fr(etaF)}
    for tag, key in (("min", "min_lambda"), ("best", "best_theta")):
        c = header["certificates"][key]
        lam, H = float.fromhex(c["lambda"]), float.fromhex(c["H"])
        gam, Lam, slack = decide(arrs["M"], arrs["Gnear"], arrs["eta"], etaF, arrs[c["w"]], lam, H, n)
        bound = iv["C_fin"] * (1 + Lam * Fr(lam) * Fr(H) / (1 - Fr(lam)))
        theta = Fr(99, 100) / bound
        out.update({f"lambda_{tag}": Fr(lam), f"H_{tag}": Fr(H), f"Gamma_w_{tag}": gam, f"Lambda_{tag}": Lam,
                    f"bound_{tag}": bound, f"theta_{tag}": theta, f"slack_i_{tag}": slack})
    out["margin"] = 1 - out["lambda_min"]
    report = {"label": run["label"], "code_base": run["code_base"], "instance": run["instance"], "near_types": n,
              "theorems": run["theorems"], "etaF_upper": sig_up(Fr(etaF), 5),
              "min_lambda": {"lambda_upper": dec_up(out["lambda_min"], 6), "Gamma_w_upper": dec_up(out["Gamma_w_min"], 3),
                             "H": dec_dn(out["H_min"], 3), "theta_star_lower": sig_dn(out["theta_min"], 4)},
              "best_theta": {"lambda_upper": dec_up(out["lambda_best"], 6),
                             "bound_upper": dec_up(out["bound_best"], 2),
                             "theta_star_lower": sig_dn(out["theta_best"], 4)},
              "near_data_sha256": header["near_data_sha256"]}
    vec = [out[k] for k in ("etaF", "lambda_min", "H_min", "Gamma_w_min", "theta_min", "lambda_best", "H_best",
                            "Gamma_w_best", "theta_best")]
    return out, report, vec


# ------------------------------------------------------------------ displays
def check_displays(displays, values):
    """Each display: scope (instance key or run label), quantity, direction (upper/lower/exact), value."""
    for item in displays:
        v = values[item["scope"]][item["quantity"]]
        shown = Fr(item["value"])
        if item["direction"] == "upper":
            assert v <= shown, item
        elif item["direction"] == "lower":
            assert v >= shown, item
        else:
            assert item["direction"] == "exact", item
            assert v == shown, item
    return len(displays)


def digest(values):
    """SHA-256 of the exact rationals, each written as hexadecimal numerator/denominator."""
    return hashlib.sha256("\n".join(f"{Fr(v).numerator:x}/{Fr(v).denominator:x}" for v in values).encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cert-dir", default=str(HERE / "certificates"))
    ap.add_argument("--runs", default=None, help="comma-separated labels (default: every run of params.json)")
    ap.add_argument("--certificates-only", action="store_true",
                    help="skip the local certificates and the displays (for freshly recomputed certificates)")
    a = ap.parse_args()
    OR.check_constants()
    cert_dir = Path(a.cert_dir).resolve()
    full = not a.certificates_only
    sharpened = json.loads((HERE / "local_sharpened.json").read_text())
    runs = PARAMS["runs"] if a.runs is None else [r for r in PARAMS["runs"] if r["label"] in a.runs.split(",")]
    needed = {r["instance"] for r in runs} if a.runs else set(PARAMS["instances"])
    values, inst_reports, local_all, cert_all = {}, {}, [], []
    for key, inst in PARAMS["instances"].items():
        if key not in needed:
            continue
        t0 = time.monotonic()
        v, rep, lv = check_instance(key, inst, sharpened, full)
        values[key], inst_reports[key] = v, rep
        local_all += lv
        extra = (f"; local inputs ({inst['local_inputs']}) certified" if full else "")
        print(f"PASS {key}: hypotheses of Theorem M, rho_- >= {rep['rho_minus_lower']}, rho_c >= {rep['rho_c_lower']},"
              f" C_fin <= {rep['C_fin_upper']}, L_0(1, 10, 100) = {rep['L0']}{extra} ({time.monotonic() - t0:.1f}s)",
              flush=True)
    run_reports = []
    for run in runs:
        t0 = time.monotonic()
        header = json.loads((cert_dir / f"{run['label']}.json").read_text())
        v, rep, vec = check_run(run, header, cert_dir, PARAMS["instances"][run["instance"]], values[run["instance"]])
        values[run["label"]] = v
        run_reports.append(rep)
        cert_all += vec
        print(f"PASS {run['label']} ({rep['near_types']} types): (i), (ii) exact for both certificates; "
              f"lambda <= {rep['min_lambda']['lambda_upper']} (theta_* >= {rep['min_lambda']['theta_star_lower']}), "
              f"best theta_* >= {rep['best_theta']['theta_star_lower']} at lambda <= {rep['best_theta']['lambda_upper']}"
              f" ({time.monotonic() - t0:.1f}s)", flush=True)
    summary = {"instances": inst_reports, "certificates": run_reports}
    if full:
        displays = [x for x in PARAMS["displays"] if x["scope"] in values]
        summary["displays_checked"] = check_displays(displays, values)
        print(f"PASS {summary['displays_checked']} displayed numbers", flush=True)
        summary["local_vector_sha256"] = digest(local_all)
    summary["certificate_vector_sha256"] = digest(cert_all)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

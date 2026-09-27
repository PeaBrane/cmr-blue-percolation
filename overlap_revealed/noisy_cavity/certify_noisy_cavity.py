"""Exact certificate for the second d = 9 proof: the oriented engine with noisy-cavity local inputs (Theorem Q9').

Standard library only. For each of the five temperatures t = 3/20 (the headline row), 29/200, 31/200, 4/25 and
7/50, with the frozen inputs of params.json and data/CERT_DATA_d9.json, this program decides in exact rational
arithmetic:

(N)  the noisy-cavity class of Lemma N: e^{2 atanh(t b)} <= wc' <= e^{2 beta}, b = tanh((2d-1) beta);
(P)  Lemma P (bond floor p <= every conditional bond probability): the 162 point conditions, the 1377 chord
     polynomials of degree 5 (Bernstein subdivision, then Sturm) and the 18 class sums Psi(j) <= 0;
(S)  Lemma S (Holley line g_h g_K^S in all 19 environments): every two-block class Gamma(l+, l-) >= 0, with the
     stored Lambda'_k recomputed and the flip symmetry checked;
(W)  plain aligned-frozen witness < p <= noisy aligned-frozen witness (so the noisy class is strictly smaller than
     the plain class at the certified floor);
(H)  g_h < 1 <= g_K, g_K g_h <= 1, g_K^d < g_K^(2d) < 2718/1000, d >= 6, and the mean-field root test
     tanh(|h| + 2dK' xbar) < xbar at the frozen xbar;
(G)  the oriented second-moment criterion (G Thm 4.5) with code base B (functions of ../certify.py): rho_-, eta',
     Score and Score_c with rho_-^2 >= c_cov/4 (G Lemma 3.5);
(G') the sharper covariance constant of G Lemma 3.6 (lemma36.py, the second implementation): c'_cov, rho'_c and
     Score'_c, with 2d Kb <= 35/100 where Lemma 3.6 is used, and the least L_* of Lemma 3.6' otherwise;
(T)  theta_* = (99/100) rho_- (1 - eta') (1 - Score_e)/(1 - 1/d) for the score used by the theorem (Score_c, or
     Score'_c at t = 7/50) and theta' with Score'_c.

It then asserts every displayed number listed in params.json ("displays"), each the worse of the research code
bases and rounded in the safe direction, and prints a JSON summary equal to ../../expected/noisy_cavity.json.
"""
from fractions import Fraction as Fr
from pathlib import Path
import hashlib
import json
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import certify as OR                 # noqa: E402  code base B of the overlap-revealed route
import lemma36 as L36                # noqa: E402
import local_certificates as LC      # noqa: E402

PARAMS = json.loads((HERE / "params.json").read_text())
DATA = json.loads((HERE / "data" / "CERT_DATA_d9.json").read_text())


def sig_dn(x, sig=4):
    x = Fr(x)
    e = 0
    while x * Fr(10) ** e < 10 ** (sig - 1):
        e += 1
    while x * Fr(10) ** e >= 10 ** sig:
        e -= 1
    y = x * Fr(10) ** e
    digits = str(y.numerator // y.denominator)
    return f"{digits[0]}.{digits[1:]}e{len(digits) - 1 - e}"


def certify_row(row, green_cache):
    d, t = row["d"], Fr(row["t"])
    p, gK, gh, xbar = (Fr(row[k]) for k in ("p", "g_K", "g_h", "xbar"))
    E = DATA[row["label"]]
    assert E["d"] == d and Fr(E["t"]) == t and Fr(E["p"]) == p
    assert Fr(E["holley"]["gK"]) == gK and Fr(E["holley"]["gh"]) == gh
    # (N), (P), (S), (W)
    P = LC.pair_certificate(d, E)
    H = LC.holley_certificate(d, E)
    assert P["plain_witness"] < p <= P["noisy_witness_beta_c"] and p <= P["noisy_witness"]
    # (H)
    (K_lo, K_hi), (H_lo, H_hi) = OR.hypotheses(d, gK, gh)
    OR.root_test(d, K_hi, H_hi, xbar)
    rho = (1 - xbar) / 2
    # (G)
    if d not in green_cache:
        green_cache[d] = OR.green(d)
    g = OR.criterion(d, p, gK, gh, rho, green_cache[d])
    # (G')
    ccp, aux = L36.lemma36(d, gK, gh, xbar)
    assert rho ** 2 >= ccp / 4
    rho_cp = rho + ccp / (4 * rho)
    score_cp = 1 / (d * rho_cp * p) + g["term2"] + g["term3"]
    lemma = "3.6" if aux["abar"] <= Fr(35, 100) else "3.6'"
    Ls, Ls_alpha = L36.least_Lstar(aux["abar"]), L36.least_Lstar(aux["alpha"])
    assert Ls <= 40                  # Lemma 3.6' changes no torus size, since L_0(n) >= 40
    # (T)
    used = g["score_c"] if row["score_used"] == "score_c" else score_cp
    theta = OR.theta_star(d, rho, g["eta"], used)
    theta_p = OR.theta_star(d, rho, g["eta"], score_cp)
    if row["score_used"] == "score_c_prime":
        assert lemma == "3.6"        # t = 7/50: G Lemma 3.6 itself applies
    M = P["model"]
    values = dict(
        p=p, noisy_witness_beta_c=P["noisy_witness_beta_c"], noisy_witness=P["noisy_witness"],
        p_over_noisy_witness=p / P["noisy_witness_beta_c"], plain_witness=P["plain_witness"], max_Psi=max(P["Psi"]),
        min_gamma=H["min_gamma"], rho_minus=rho, kappa=g["kappa"], alpha=g["alpha"], two_d_Kprime=aux["abar"],
        eta_prime=g["eta"], score=g["score"], score_c=g["score_c"], score_c_prime=score_cp, theta_star=theta,
        theta_prime=theta_p, c_cov=g["c_cov"], rho_c=g["rho_c"], c_cov_prime=ccp, rho_c_prime=rho_cp,
        least_Lstar=Fr(Ls), least_Lstar_alpha=Fr(Ls_alpha), Kprime_lo=K_lo, Kprime_hi=K_hi, abs_h_lo=H_lo,
        abs_h_hi=H_hi, margin_c=1 - g["score_c"], margin_c_prime=1 - score_cp, term1=g["term1"],
        term1_c=1 / (d * g["rho_c"] * p), term1_c_prime=1 / (d * rho_cp * p), term2=g["term2"], term3=g["term3"],
        G_lo=g["G_lo"], G_hi=g["G_hi"], F=g["F"], Sh2=g["Sh2"], U1=g["U1"], SU=g["SU"],
        nine_rho_c_p=d * g["rho_c"] * p, nine_rho_c_prime_p=d * rho_cp * p, prefactor_base=M.w / (1 + M.w ** 2),
        **{f"Psi_{j}": v for j, v in enumerate(P["Psi"])},
        **{f"min_gamma_k{r['k']}": r["min_gamma"] for r in H["rows"]})
    worst = min(H["rows"], key=lambda r: r["min_gamma"])
    report = {
        "label": row["label"], "theorem": row["theorem"], "t": row["t"], "p": row["p"], "g_K": row["g_K"],
        "g_h": row["g_h"], "xbar": row["xbar"], "wc_prime": E["wc_prime"],
        "pair_certificate": {"point_conditions": P["points"], "chord_polynomials": P["chords"],
                             "max_Psi_upper": "-" + OR.dec_dn(-max(P["Psi"]), 2),
                             "noisy_witness_upper": OR.dec_up(P["noisy_witness_beta_c"], 9),
                             "plain_witness_upper": OR.dec_up(P["plain_witness"], 9),
                             "p_over_noisy_witness_lower": OR.dec_dn(p / P["noisy_witness_beta_c"], 6)},
        "holley_line": {"environments": len(H["rows"]), "classes": sum(r["classes"] for r in H["rows"]),
                        "min_Gamma_lower": sig_dn(H["min_gamma"], 3), "at_S": worst["S"],
                        "at_class": list(worst["argmin"])},
        "rho_minus_lower": OR.dec_dn(rho, 7), "kappa_upper": OR.dec_up(g["kappa"], 6),
        "alpha_upper": OR.dec_up(g["alpha"], 6), "two_d_Kprime_upper": OR.dec_up(aux["abar"], 6),
        "eta_prime_upper": OR.dec_up(g["eta"], 7), "score_upper": OR.dec_up(g["score"], 7),
        "score_c_upper": OR.dec_up(g["score_c"], 7), "c_cov_prime_lower": OR.dec_dn(ccp, 8),
        "rho_c_prime_lower": OR.dec_dn(rho_cp, 8), "score_c_prime_upper": OR.dec_up(score_cp, 7),
        "covariance_lemma": lemma, "least_Lstar": Ls, "score_used": row["score_used"],
        "theta_star_lower": OR.dec_dn(theta, 6), "theta_prime_lower": OR.dec_dn(theta_p, 6)}
    vector = P["Psi"] + [P["noisy_witness_beta_c"], P["plain_witness"]] + [r["min_gamma"] for r in H["rows"]] + \
        [rho, g["eta"], g["score"], g["score_c"], ccp, score_cp, theta, theta_p]
    return values, report, vector


def check_displays(displays, values):
    for item in displays:
        v, shown = values[item["scope"]][item["quantity"]], Fr(item["value"])
        if item["direction"] == "upper":
            assert v <= shown, item
        elif item["direction"] == "lower":
            assert v >= shown, item
        else:
            assert item["direction"] == "exact", item
            assert v == shown, item
    return len(displays)


def main():
    OR.check_constants()
    green_cache, values, rows, vector = {}, {}, [], []
    for row in PARAMS["rows"]:
        start = time.monotonic()
        v, rep, vec = certify_row(row, green_cache)
        values[row["label"]] = v
        rows.append(rep)
        vector += vec
        used = rep["score_c_upper"] if row["score_used"] == "score_c" else rep["score_c_prime_upper"]
        print(f"PASS d=9 {row['label']} ({row['theorem']}): p_B >= {row['p']} (Lemma P), Holley line in "
              f"{rep['holley_line']['environments']} environments (Lemma S, min Gamma >= "
              f"{rep['holley_line']['min_Gamma_lower']}), {row['score_used']} <= {used}, theta_* >= "
              f"{rep['theta_star_lower']} ({time.monotonic() - start:.1f}s)", flush=True)
    n = check_displays(PARAMS["displays"], values)
    print(f"PASS {n} displayed numbers", flush=True)
    digest = hashlib.sha256("\n".join(f"{Fr(x).numerator:x}/{Fr(x).denominator:x}" for x in vector).encode())
    print(json.dumps({"certified_rows": len(rows), "rows": rows, "displays_checked": n,
                      "value_vector_sha256": digest.hexdigest()}, indent=2))


if __name__ == "__main__":
    main()

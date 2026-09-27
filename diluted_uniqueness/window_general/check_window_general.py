"""Default light check of the sparse-insulation windows (d = 2 box, d = 3 cube, octahedral rule for d >= 3).

Standard library only; about 15 s, at most about 150 MB per process. Decides, in exact rational arithmetic:

 1. The rule tables params/lp_3_cube_rho_256.json (T_3) and params/lp_2_box_rho.json (T_2) have the
    recorded SHA-256 hashes.
 2. cr_check.py on both tables (type S at L = 5, 6, 7): every design satisfies the design-lemma
    hypotheses, the configuration lists are complete, the weights are positive and sum to 1; the exact
    clique constants K_B, K_N, K_BN, the type-S bounds, q, N_Q and the coin counts, and the windows
        d = 3 cube:  Delta_3 >= 3.232324e-9 (plaquette coins), >= 7.405526e-9 (vertex coins),
        d = 2 box:   Delta_2 >= 4.552624e-8 (plaquette coins), >= 7.153991e-8 (domino coins),
    with g_d = 2 Delta_d, and the classical-window thresholds beta_*(d) = (1/2) ln(p_c^- / Delta_d).
 3. oct_check.py for d = 3 at L = 5 (brackets "old" and "R8"): every configuration of the
    octahedral hand rule gets designs satisfying the hypotheses with the costs of its case, the same-edge
    key pairs are incompatible except N-int pairs at different terminals, the realizable-key sum per
    labelled plaquette is at most K_d, and the exact bulk clique value is at most K_d.
 4. The octahedral closed form for Delta_d: the counts q(d), N_Q(d) (formula = enumeration for d = 2..7), and every
    displayed entry of the octahedral table (d = 3..8, 10, 20, 100) and of the axis-star-coin table
    (d = 3..8, 10, 12), recomputed here from gd_table.window and gd_coins.counts: exact brackets, q, N_Q,
    q_c, N_c, and each displayed window or threshold is the floor of the exact value (displayed value
    <= exact value < displayed value + one unit of its last digit); log2 K_d to the displayed 3 decimals.
The printed outputs of cr_check.py, oct_check.py, gd_table.py and gd_coins.py are also compared line by
line with expected_outputs/ (the outputs of the audited runs). That comparison is a regression check;
the decisions above are assertions on recomputed exact values. The heavier d = 4, 5 octahedral checks,
d = 3 at L = 6, the networkx clique cross-check and the end-to-end sanity runs are in
reproduce_window_general.py.
"""

from fractions import Fraction as Fr
from pathlib import Path
import contextlib
import hashlib
import io
import json
import math
import os
import subprocess
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cr_check  # noqa: E402
import gd_coins  # noqa: E402
import gd_table  # noqa: E402

TABLES = {
    "lp_3_cube_rho_256.json": "c43d5de37ae33e1f71156290939dd2555e62e1704ad1a81bcc35cc4053f8d149",
    "lp_2_box_rho.json": "8e23d624fe8fe11f152eda6d590d185da0988ae3374c06578a32f467f23e93fa",
}

# Octahedral table as displayed: d: (p_-, p_+, q, N_Q, log2 K_d, Delta_d >=, g_d >=, beta_*(d) >=,
# g_d >= with p_- = 1/(2d), p_+ = 7/20). p_+ is exact; p_- is rounded to 6 decimals (exact values below).
TABLE_R8 = {
    3: ("0.211027", "0.347298", 312, 312, "23.723", "2.3067e-10", "4.6135e-10", "10.317", "2.3957e-10"),
    4: ("0.146971", "0.278746", 1536, 1024, "27.385", "3.7128e-12", "7.4257e-12", "12.200", "2.0797e-12"),
    5: ("0.1", "0.228327", 5200, 2600, "30.206", "1.5528e-13", "3.1056e-13", "13.595", "1.0464e-14"),
    6: ("0.083333", "0.192101", 13920, 5568, "32.334", "1.3267e-14", "2.6535e-14", "14.734", "6.0968e-17"),
    7: ("0.071429", "0.168116", 31752, 10584, "34.184", "1.6140e-15", "3.2281e-15", "15.710", "3.9579e-19"),
    8: ("0.0625", "0.148514", 64512, 18432, "35.734", "2.7123e-16", "5.4247e-16", "16.535", "2.7947e-21"),
    10: ("0.05", "0.120502", 208800, 46400, "38.278", "1.4374e-17", "2.8748e-17", "17.892", "1.6890e-25"),
    20: ("0.025", "0.062039", 7539200, 793600, "45.872", "2.0603e-21", "4.1207e-21", "21.971", "1.0330e-45"),
    100: ("0.005", "0.012719", 26009280000, 525440000, "62.616", "5.4386e-30", "1.0877e-29", "31.042",
          "6.9823e-199"),
}
P_MINUS_EXACT = {3: Fr(10000, 47387) - Fr(1, 10 ** 6), 4: Fr(10000, 68040) - Fr(1, 10 ** 6)}
# Axis-star coins: d: (q, q_c, N_Q, N_c, g_d >= (plaquette coins), g_d >= (axis-star coins), gain)
TABLE_R8_COINS = {
    3: (312, 226, 312, 366, "4.6135e-10", "6.3614e-10", None),
    4: (1536, 914, 1024, 1392, "7.4257e-12", "1.2473e-11", "1.680"),
    5: (5200, 2726, 2600, 4060, "3.1056e-13", "5.9232e-13", "1.907"),
    6: (13920, 6688, 5568, 9876, "2.6535e-14", "5.5224e-14", "2.081"),
    7: (31752, 14320, 10584, 21084, "3.2281e-15", "7.1576e-15", "2.217"),
    8: (64512, 27734, 18432, 40816, "5.4247e-16", "1.2618e-15", "2.326"),
    10: (208800, 83904, 46400, 123720, "2.8748e-17", "7.1540e-17", "2.489"),
    12: (540672, 207658, 98304, 307104, "2.7337e-18", "7.1178e-18", "2.604"),
}
# Finer admissible bracket p_+ = hi_60 + 10^-6 (displayed, not used elsewhere), and beta_*(5) with p_c >= 1/9.
FINE_G = {4: "7.4258e-12", 5: "3.1057e-13", 20: "4.1209e-21", 100: "1.0878e-29"}


def atanh_bounds(z, terms=60):
    lo = sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(terms))
    return lo, lo + z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))


def ln_bounds(x):
    """Enclosure of ln x for rational x > 0; x = 2^e r, r in [1, 2) cut to 80 bits in both directions."""
    x = Fr(x)
    e = x.numerator.bit_length() - x.denominator.bit_length()
    if Fr(2) ** e > x:
        e -= 1
    r = x / Fr(2) ** e
    r_lo, r_hi = Fr(math.floor(r * 2 ** 80), 2 ** 80), Fr(math.ceil(r * 2 ** 80), 2 ** 80)
    l2_lo, l2_hi = (2 * t for t in atanh_bounds(Fr(1, 3)))
    a_lo = 2 * atanh_bounds((r_lo - 1) / (r_lo + 1))[0]
    a_hi = 2 * atanh_bounds((r_hi - 1) / (r_hi + 1))[1]
    return (e * (l2_lo if e >= 0 else l2_hi) + a_lo, e * (l2_hi if e >= 0 else l2_lo) + a_hi)


def log2_bounds(x):
    lo, hi = ln_bounds(x)
    l2_lo, l2_hi = (2 * t for t in atanh_bounds(Fr(1, 3)))
    return lo / l2_hi, hi / l2_lo


def unit(shown):
    """One unit of the last displayed digit of a decimal string such as '2.3067e-10' or '10.317'."""
    mant, _, exp = shown.partition("e")
    places = len(mant.split(".")[1]) if "." in mant else 0
    return Fr(10) ** ((int(exp) if exp else 0) - places)


def is_floor(x, shown):
    """shown <= x < shown + one unit of its last digit."""
    return Fr(shown) <= x < Fr(shown) + unit(shown)


def rounds_to(lo, hi, shown):
    """An enclosure [lo, hi] lies within half a unit of the last digit of shown (rounding to nearest)."""
    half = unit(shown) / 2
    return Fr(shown) - half <= lo and hi <= Fr(shown) + half


def log2_in(x, shown):
    return rounds_to(*log2_bounds(x), shown)


def beta_floor(pc_lo, delta, shown):
    """beta_* = (1/2) ln(pc_lo / delta): shown <= beta_* < shown + one unit (decided with an enclosure)."""
    lo, hi = ln_bounds(pc_lo / delta)
    return Fr(shown) <= lo / 2 and hi / 2 < Fr(shown) + unit(shown)


def capture(fn):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        val = fn()
    return buf.getvalue(), val


def run_script(name, args):
    """Run a checker in a child process of this interpreter (keeps the peak memory of each run separate)."""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, str(HERE / name), *map(str, args)], cwd=HERE, env=env, check=True,
                          capture_output=True, text=True).stdout


def expect(text, name):
    want = (HERE / "expected_outputs" / name).read_text()
    assert text == want, f"output differs from expected_outputs/{name}"


def run_cr_check(table, name):
    saved = sys.argv
    sys.argv = ["cr_check.py", str(HERE / "params" / table), "5", "6", "7"]
    try:
        text, res = capture(cr_check.main)
    finally:
        sys.argv = saved
    expect(text, name)
    return res


def dec_floor(x, places):
    q = math.floor(x * 10 ** places)
    return f"{q // 10 ** places}.{q % 10 ** places:0{places}d}"


def check_cube(report):
    t = time.monotonic()
    r = run_cr_check("lp_3_cube_rho_256.json", "cr_check_d3.out")
    assert (r["d"], r["p_lo"], r["p_hi"]) == (3, Fr(211, 1000), Fr(3473, 10000))
    assert (r["q"], r["NQ"], r["coins"][:2]) == (180, 180, (78, 234))
    assert r["K"] == r["Ks"]["B"] == r["ep_vals"]["B"][0] >= max(r["Ks"]["N"], r["Ks"]["BN"], r["Ks"]["S"])
    assert log2_in(r["K"], "20.704819") and log2_in(r["ep_vals"]["B"][1], "20.704816")
    assert log2_in(r["Ks"]["N"], "15.2143") and abs(r["K"] - Fr("1709115.617")) < Fr(1, 1000)
    # Type S: the realized keys of an output contribute at most 6 x (worst key) (six edges per plaquette in
    # the cube regions). 6 x worst = 2^14.64083 > 2^14.6408, so the safe display is K_S <= 2^14.6409.
    assert log2_in(r["worstS"], "12.0559")
    lo, hi = log2_bounds(6 * r["worstS"])
    assert Fr("14.6408") < lo and hi < Fr("14.6409")
    assert log2_bounds(r["Ks"]["S"])[1] < Fr("14.6408")
    assert r["Delta"] == cr_check.window(180, 180, r["K"], r["p_hi"])
    assert r["Delta_coins"] == cr_check.window(78, 234, r["K"], r["p_hi"])
    assert is_floor(r["Delta"], "3.232324e-9") and is_floor(2 * r["Delta"], "6.464648e-9")
    assert is_floor(r["Delta_coins"], "7.405526e-9") and is_floor(2 * r["Delta_coins"], "1.481105e-8")
    pc3 = Fr(10000, 47387)              # p_c(Z^3) >= 1/mu_3 >= 1/4.7387
    assert beta_floor(pc3, r["Delta_coins"], "8.5826") and beta_floor(pc3, r["Delta"], "8.9971")
    assert 2 / ln_bounds(pc3 / r["Delta_coins"])[0] < Fr("0.1166")
    print(f"PASS d=3 cube certificate (plaquette and vertex coins): K = K_B = 2^20.704819, Delta_3 >= 3.232324e-9, "
          f"vertex coins >= 7.405526e-9, K_S <= 2^14.6409 ({time.monotonic() - t:.1f}s)", flush=True)
    report["cube_d3"] = {
        "p_bracket": ["211/1000", "3473/10000"], "q": 180, "N_Q": 180, "q_c": 78, "N_c": 234,
        "log2_K_B": "20.704819", "log2_K_N": "15.2143", "log2_K_S_bound": "14.6409",
        "Delta_plaquette_coins_floor": "3.232324e-9", "g_plaquette_coins_floor": "6.464648e-9",
        "Delta_vertex_coins_floor": "7.405526e-9", "g_vertex_coins_floor": "1.481105e-8",
        "beta_star_vertex_coins_floor": "8.5826", "beta_star_plaquette_coins_floor": "8.9971",
    }


def check_box(report):
    t = time.monotonic()
    r = run_cr_check("lp_2_box_rho.json", "cr_check_d2.out")
    assert (r["d"], r["p_lo"], r["p_hi"]) == (2, Fr(499999, 10 ** 6), Fr(500001, 10 ** 6))
    assert (r["q"], r["NQ"], r["coins"][:2]) == (32, 64, (20, 76))
    assert r["K"] == r["Ks"]["B"] == r["ep_vals"]["B"][1] >= max(r["Ks"]["N"], r["Ks"]["BN"], r["Ks"]["S"])
    assert log2_in(r["K"], "19.344193") and log2_in(r["ep_vals"]["B"][0], "19.344182")
    assert log2_in(r["Ks"]["N"], "12.8191") and abs(r["K"] - Fr("665552.744")) < Fr(1, 1000)
    # Type S: at most 8 x (worst key) = 2^12.8074 <= 2^12.81; computed sum over all keys 2^11.7748.
    assert log2_in(r["worstS"], "9.8074") and log2_bounds(8 * r["worstS"])[1] < Fr("12.81")
    assert log2_in(r["Ks"]["S"], "11.7748")
    assert r["Delta"] == cr_check.window(32, 64, r["K"], r["p_hi"])
    assert r["Delta_coins"] == cr_check.window(20, 76, r["K"], r["p_hi"])
    assert is_floor(r["Delta"], "4.552624e-8") and is_floor(2 * r["Delta"], "9.105248e-8")
    assert is_floor(r["Delta_coins"], "7.153991e-8") and is_floor(2 * r["Delta_coins"], "1.430798e-7")
    assert beta_floor(Fr(1, 2), r["Delta_coins"], "7.8799")        # p_c(Z^2) = 1/2
    assert 2 / ln_bounds(Fr(1, 2) / r["Delta_coins"])[0] < Fr("0.1270")
    print(f"PASS d=2 box certificate (plaquette and domino coins): K = K_B = 2^19.344193 at p_hi, Delta_2 >= 4.552624e-8, "
          f"domino coins >= 7.153991e-8 ({time.monotonic() - t:.1f}s)", flush=True)
    report["box_d2"] = {
        "p_bracket": ["499999/1000000", "500001/1000000"], "q": 32, "N_Q": 64, "q_c": 20, "N_c": 76,
        "log2_K_B": "19.344193", "log2_K_N": "12.8191", "log2_K_S_bound": "12.81",
        "Delta_plaquette_coins_floor": "4.552624e-8", "g_plaquette_coins_floor": "9.105248e-8",
        "Delta_domino_coins_floor": "7.153991e-8", "g_domino_coins_floor": "1.430798e-7",
        "beta_star_domino_coins_floor": "7.8799",
    }


def check_octahedral(report):
    t = time.monotonic()
    for args, name in (((3, 5, 1, 1, "old"), "oct_check_d3_L5_old.out"), ((3, 5, 1, 1, "R8"), "oct_check_d3_L5_R8.out")):
        expect(run_script("oct_check.py", args), name)
    print(f"PASS octahedral rule, d=3, L=5 (rule validity, costs, key multiplicity; brackets old and R8) ({time.monotonic() - t:.1f}s)",
          flush=True)
    report["octahedral_rule_checked"] = {"d": 3, "L_S": [5], "brackets": ["old", "R8"]}


def check_closed_form(report):
    t = time.monotonic()
    expect(run_script("gd_table.py", []), "gd_table.out")
    expect(run_script("gd_coins.py", [12]), "gd_coins.out")
    rows = {}
    for d, (pm, pp, q, nq, log2k, delta, g, beta, g_univ) in TABLE_R8.items():
        q0, nq0, p_lo, p_hi, K, D = gd_table.window(d)
        assert (q0, nq0) == (q, nq) == gd_table.q_NQ(d)
        assert p_hi == Fr(pp) and p_lo == P_MINUS_EXACT.get(d, Fr(1, 2 * d))
        assert abs(p_lo - Fr(pm)) <= unit("0.000001") / 2
        assert rounds_to(*log2_bounds(K), log2k)
        assert is_floor(D, delta) and is_floor(2 * D, g)
        assert beta_floor(p_lo, D, beta)
        assert is_floor(2 * gd_table.window(d, universal=True)[5], g_univ)
        if d in FINE_G:
            assert is_floor(2 * gd_table.window(d, fine=True)[5], FINE_G[d])
        rows[str(d)] = {"p_minus": str(p_lo), "p_plus": str(p_hi), "q": q, "N_Q": nq, "Delta_floor": delta,
                        "g_floor": g, "beta_star_floor": beta, "g_floor_universal_bracket": g_univ}
    assert beta_floor(Fr(1, 9), gd_table.window(5)[5], "13.648")
    # observation about the tabulated range (all rows of gd_table.py): 1.0e-5 <= d^12 g_d <= 2.5e-4
    assert all(Fr("1.0e-5") <= d ** 12 * 2 * gd_table.window(d)[5] <= Fr("2.5e-4")
               for d in list(range(3, 13)) + [15, 20, 30, 50, 100])
    report["octahedral_table"] = rows
    coin_rows = {}
    for d, (q, qc, nq, nc, g, gc, gain) in TABLE_R8_COINS.items():
        q1, qc1, nc1 = gd_coins.counts(d)
        q0, nq0, p_lo, p_hi, K, D = gd_table.window(d)
        assert (q1, qc1, nq0, nc1) == (q, qc, nq, nc) and q0 == q
        Dc = cr_check.window(qc, nc, K, p_hi)
        assert is_floor(2 * D, g) and is_floor(2 * Dc, gc)
        if gain:
            assert rounds_to(Dc / D, Dc / D, gain)
        coin_rows[str(d)] = {"q_c": qc, "N_c": nc, "g_floor_axis_star": gc}
    report["axis_star_coins"] = coin_rows
    print(f"PASS octahedral closed form: covering counts, the octahedral table and the axis-star-coin table, "
          f"every displayed entry ({time.monotonic() - t:.1f}s)", flush=True)


def main():
    start = time.monotonic()
    report = {"rule_table_sha256": dict(TABLES)}
    for table, digest in TABLES.items():
        assert hashlib.sha256((HERE / "params" / table).read_bytes()).hexdigest() == digest
    print("PASS rule tables T_3 and T_2: SHA-256 as recorded", flush=True)
    check_cube(report)
    check_box(report)
    check_octahedral(report)
    check_closed_form(report)
    print(f"PASS window_general light check ({time.monotonic() - start:.1f}s)")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

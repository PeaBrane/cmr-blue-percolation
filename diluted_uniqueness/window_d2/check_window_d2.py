"""Default light check of the d = 2 exact local-ratio window g_2 = 2 g* (standard library only).

Decides, with exact integers and rationals:

 1. Constants. The frozen inputs params/cert.json satisfy (C0); g* is the displayed rational;
    g_2 = 2 g*; the displayed enclosures of g*, g_2 and beta_*(2) = (1/2) ln(1/(2 g*)), the gain over
    g_2 = 2^-988/1130, and the parameter range [1/2 - eta, 1/2 + eta] of both paths.
 2. Geometry. The spacing-3 pattern (stars disjoint, 4 P-edges and 8 stubs per cell), the weighted
    coverage identity of the two Path-1 families on all 1250 edges of [-12,12]^2, and the engine spec
    files: params/specs/*.spec are regenerated and compared byte for byte.
 3. Special boxes near S_L: the S_L-centred intersection box (25 data, a superset of the realizable
    ones; max N_Gc/N_T = 3/4) and the box census for L = 9, 12, 15, 30.
 4. The committed record certificates certificates/*.canon.xz. Each file's SHA-256 equals the
    canonical digest recorded with the certificate, the records are strictly sorted, and the families
    are complete: all 218,790 side, 92,378 origin and 12,012 intersection boundary conditions, and
    for the cell family a complete system of D4-orbit representatives (26,584 partitions, 496,272
    records, sum over orbits of |orbit| C(k,2) = 3,879,876).
 5. (C1-side), (C1-o), (C1-int) on every Path-1 record and (C2) at u = 0 and u = 1/2 on every
    reduced cell record, exactly. With the D4 invariance of the counts (proved by hand), (C2) then
    holds for all 3,879,876 cell boundary conditions. Also the displayed maxima, the 20 binding cell
    boundary conditions, the binding representative, and the 69,292 cell boundary conditions with a
    pivotal grid edge and D = 0.
 6. Record correctness, partly. An independent Python implementation of the record definition
    recomputes the whole intersection family (its canonical digest must equal the committed one) and
    every record of a deterministic subset of the side and origin partitions (the partitions of the
    tightest records and every 32nd partition).

Not decided here: that the side and origin records outside that subset, and the cell records, are the
counts that their definition prescribes. The full reproduction (reproduce_window_d2.py: the C engine on
all families, including all 208,012 cell partitions, plus sampled literal brute force) decides that.
With --with-cc this check also compiles engine/jbrute.c and brute-forces the binding reduced cell records
and a few fixed-seed random ones (about 2 s each), deciding those cell records independently.
Runs in about 40 s with CPython 3.13 and needs 150-200 MB.
usage: python check_window_d2.py [--with-cc] [--brute-random N]
"""

from fractions import Fraction as F
from pathlib import Path
import contextlib
import argparse
import io
import json
import os
import random
import subprocess
import sys
import tempfile
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geom  # noqa: E402
import lrm_common as L  # noqa: E402
import sl_box  # noqa: E402

G_STAR = F(1059411733898907047504026797, 90604240507898863756433027512000)
A_P = F(1044434688759, 52158068414761700)
DIGESTS = {
    "side": "575c9bbdf87bc8282c1d36d0d39983df15defd829ba63a2f885fe011d58ea2aa",
    "sideo": "42f054a28ea9c56adf68561e3acbb67f70617e2ca2a84f8ca3080725e4335a4b",
    "int": "946bc6e254d49eefa36b5b11a4fdb6b1db91b285163648a752d5f3b8b5f76a2f",
    "cell_reps": "d51fc9f291e210791aac390acedba4fc6e7260ce895fbee7fa177253e8242ab7",
}
RECORDS = {"side": 218790, "sideo": 92378, "int": 12012, "cell_reps": 496272}
SUBSET_STRIDE = 32


def rounds_to(x, shown):
    """x lies within half a unit of the last displayed digit of `shown` (a decimal string)."""
    mant = shown.split("e")[0].lstrip("-")
    exp = int(shown.split("e")[1]) if "e" in shown else 0
    places = len(mant.split(".")[1]) if "." in mant else 0
    return abs(x - F(shown)) <= F(1, 2) * F(10) ** (exp - places)


def step(msg, start):
    print(f"PASS {msg} ({time.monotonic() - start:.1f}s)", flush=True)


def check_constants(cert, report):
    t = time.monotonic()
    L.check_c0(cert)
    g, aP, aT, eta = cert["g"], cert["aP"], cert["aT"], cert["eta"]
    assert g == G_STAR and aP == A_P
    g2 = 2 * g
    assert F("1.1692738970716e-5") <= g <= F("1.1692738970717e-5")
    assert F("2.3385477941433e-5") <= g2 <= F("2.3385477941434e-5")
    lo, hi = L.log_bounds(1 / g2)
    beta_lo, beta_hi = lo / 2, hi / 2
    assert F("5.3316976645232") <= beta_lo <= beta_hi <= F("5.3316976645233")
    assert F("0.187555") < 1 / beta_hi <= 1 / beta_lo < F("0.187565")
    gA = F(1, 2 ** 988 * 1130)
    l_lo, l_hi = L.log_bounds(g2 / gA)
    ten_lo, ten_hi = L.log_bounds(10)
    assert F("295.835") <= l_lo / ten_hi <= l_hi / ten_lo < F("295.845")   # displayed as 10^295.84
    assert -eta <= min(-g, -aT) and g <= eta
    for x, shown in ((aP, "2.0024412723e-5"), (aT, "3.6043942901e-5"),
                     (cert["c_side"][0], "-9.241729e-7"), (cert["c_side"][2], "-8.331674e-6"),
                     (cert["c_side"][3], "-7.994479e-6"), (cert["c_int"][0], "1.261691e-5"),
                     (cert["c_int"][3], "-1.635672e-5")):
        assert rounds_to(x, shown), shown
    report["constants"] = {
        "g_star": f"{g.numerator}/{g.denominator}",
        "g_star_bracket": [L.sci(g, 14, False), L.sci(g, 14, True)],
        "g_2": "2 g_star",
        "g_2_bracket": [L.sci(g2, 14, False), L.sci(g2, 14, True)],
        "beta_star_2": [L.dec(beta_lo, 13, False), L.dec(beta_hi, 13, True)],
        "log10_gain_over_2^-988/1130": [L.dec(l_lo / ten_hi, 4, False), L.dec(l_hi / ten_lo, 4, True)],
    }
    step("constants: (C0), g* exact, g_2 = 2g*, beta_*(2) in [5.3316976645232, 5.3316976645233]", t)


def check_geometry(report):
    t = time.monotonic()
    with contextlib.redirect_stdout(io.StringIO()) as log, tempfile.TemporaryDirectory() as tmp:
        saved, sys.argv = sys.argv, ["geom.py", tmp]
        try:
            geom.main()
        finally:
            sys.argv = saved
        for name in ("cell", "side", "sideo", "int"):
            fresh = Path(tmp, f"{name}.spec").read_bytes()
            assert fresh == (HERE / "params" / "specs" / f"{name}.spec").read_bytes(), name
    text = log.getvalue()
    assert "coverage side family (Gc,Gm,P,T): 1250 edges tested in [-12,12]^2, failures 0" in text
    assert "coverage intersection family (Gc,T): 1250 edges tested in [-12,12]^2, failures 0" in text
    report["geometry"] = {"coverage_edges_tested": 1250, "spec_files_regenerated": 4}
    step("geometry: pattern, weighted coverage identity (1250 edges per family), 4 spec files regenerated", t)


def check_special(cert_json, cert, report):
    t = time.monotonic()
    with contextlib.redirect_stdout(io.StringIO()) as log:
        assert sl_box.sl_centred(cert_json)
        for Lb in (9, 12, 15, 30):
            assert sl_box.census(Lb)
    text = log.getvalue()
    assert ("(4) census L=9: boxes {('int', 'S_L-centred'): 20, ('int', 'corner'): 4, ('int', 'generic'): 25, "
            "('int', 'irrelevant'): 72, ('side', 'generic'): 59, ('side', 'irrelevant'): 158, "
            "('side', 'origin'): 1, ('side', 'straddling'): 24}") in text
    assert "25 boundary conditions" in text
    eta, ci = cert["eta"], cert["c_int"]
    assert -ci[3] * (1 - 2 * eta) ** 4 > F("1.295") * ci[0] * (1 + 2 * eta) ** 4
    report["special_boxes"] = {"S_L_centred_data": 25, "max_NGc_over_NT": "3/4", "census_L": [9, 12, 15, 30]}
    step("special boxes: S_L-centred box (25 data, max N_Gc/N_T = 3/4), census L = 9, 12, 15, 30", t)


def load(family):
    npos = L.FAMILIES[family][1]
    canon = L.Canon.from_xz(HERE / "certificates" / f"{family}.canon.xz", npos, RECORDS[family])
    assert canon.sha256 == DIGESTS[family], family
    assert canon.strictly_sorted(), family
    return canon


def check_path1(cert, report):
    t = time.monotonic()
    shown = {"side": ("0.651213", 218790, 16796), "sideo": ("0.978415", 92378, 16796),
             "int": ("0.993495", 12012, 1430)}
    canons = {}
    for fam, (ratio, nrec, npart) in shown.items():
        nb, _, origin, _ = L.FAMILIES[fam]
        canon = load(fam)
        ok, nparts, npairs = L.completeness(canon, nb, origin)
        assert ok and nparts == npart and npairs == nrec == canon.n, fam
        k_sum = sum(L.narayana(nb, k) * (k if origin else k * (k - 1) // 2) for k in range(1, nb + 1))
        assert k_sum == nrec
        res = L.check_c1(canon.rows(), cert, fam)
        assert res["violations"] == 0 and res["nonzero_cell_slots"] == 0, fam
        assert rounds_to(res["max_ratio"], ratio), (fam, float(res["max_ratio"]))
        report[f"C1_{fam}"] = {"records": canon.n, "canonical_sha256": canon.sha256,
                               "max_ratio": L.dec(res["max_ratio"], 9, True)}
        canons[fam] = (canon, res)
    step("(C1-side), (C1-o), (C1-int) exactly on all 218,790 + 92,378 + 12,012 records; digests and "
         "completeness", t)
    return canons


def orbit_table():
    """D4 orbits of the 208,012 non-crossing partitions of the 12 cell-boundary positions."""
    assert L.validate_generator(9)
    ncs = L.gen_noncrossing(12)
    assert len(ncs) == len(set(ncs)) == 208012
    all_bcs = sum(L.comb(max(p) + 1, 2) for p in ncs)
    perms = L.d4_permutations()
    orbit_of, sizes = {}, []
    for p in ncs:
        key = bytes(p)
        if key in orbit_of:
            continue
        imgs = {bytes(L.image(p, g)[0]) for g in perms}
        for q in imgs:
            orbit_of[q] = len(sizes)
        sizes.append(len(imgs))
    del ncs
    assert len(orbit_of) == 208012
    return all_bcs, perms, orbit_of, sizes


def check_cells(cert, report):
    t = time.monotonic()
    canon = load("cell_reps")
    all_bcs, perms, orbit_of, sizes = orbit_table()
    npairs = {}
    for lab, pairs in canon.groups():
        tup = tuple(lab)
        assert tup == L.rgs(tup) and L.noncrossing(tup)
        k = max(tup) + 1
        assert pairs == [(x, y) for x in range(k) for y in range(x + 1, k)]
        npairs[lab] = len(pairs)
    reps = list(npairs)
    hit = [orbit_of[lab] for lab in reps]
    # the one-block partition is an orbit of its own and has no block pair, hence no record
    trivial = orbit_of[bytes(12)]
    assert sizes[trivial] == 1 and bytes(12) not in npairs
    assert len(hit) == len(set(hit)) == len(sizes) - 1 == 26583 and trivial not in hit
    total = sum(sizes[orbit_of[lab]] * npairs[lab] for lab in reps)
    assert total == all_bcs == 3879876
    assert canon.n == 496272
    step("reduced cell list: one representative of each of the 26,584 D4 orbits (26,583 with a block pair), "
         "496,272 records, covering 3,879,876 boundary conditions", t)

    t = time.monotonic()
    res = L.check_c2(canon.rows(), cert)
    assert res["violations"] == 0 and res["D0_with_star_pivot"] == 0
    assert rounds_to(res["max_ratio_u0"], "0.993503") and rounds_to(res["max_ratio_u_half"], "0.993729")
    assert 1 - res["max_ratio_u_half"] >= F("6.2e-3")
    mult = {lab: sizes[orbit_of[lab]] for lab in reps}
    del orbit_of
    binding = {}
    for arg in ("argmax_u0", "argmax_u_half"):
        bc = sum(mult[lab] for lab, _, _, _ in res[arg])
        by_d = {}
        for lab, _, _, c in res[arg]:
            by_d[c[0]] = by_d.get(c[0], 0) + mult[lab]
        binding[arg] = {"reduced_records": len(res[arg]), "boundary_conditions": bc,
                        "by_D": {str(k): v for k, v in sorted(by_d.items())}}
        assert len(res[arg]) == 4 and bc == 20 and by_d == {276: 8, 552: 8, 828: 4}, binding[arg]
    grid_d0 = sum(mult[lab] for lab, x, y, c in canon.rows()
                  if c[0] == 0 and c[1] + c[2] + c[5] + c[6] > 0)
    assert grid_d0 == 69292
    # binding representative: pi = 0,1,2,...,9,2,0; X = {(0,0),(0,1)}, Y = {(1,0)}
    pi = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 0)
    found = []
    for g in perms:
        img, bmap = L.image(pi, g)
        if bytes(img) in npairs:
            x, y = sorted((bmap[0], bmap[1]))
            found.append(canon.counts(canon.find(bytes(img) + bytes((x, y)))))
    assert found and all(c == found[0] for c in found)
    c = found[0]
    assert (c[0], c[1], c[3], c[4], c[7], c[8]) == (552, 4272736, 498192, 3512736, 498192, 3514464)
    report["C2_cells"] = {
        "reduced_records": canon.n, "boundary_conditions_covered": total, "canonical_sha256": canon.sha256,
        "max_ratio_u0": L.dec(res["max_ratio_u0"], 9, True),
        "max_ratio_u_half": L.dec(res["max_ratio_u_half"], 9, True),
        "binding": binding, "grid_pivot_with_D0": grid_d0,
    }
    step("(C2) exactly on all 496,272 reduced cell records (u = 0, 1/2); maxima 0.993503, 0.993729; "
         "20 binding boundary conditions", t)
    return canon, res


def check_brute(canon, res, cc, nrand, report):
    """Optional (--with-cc): the literal brute force engine/jbrute.c, which enumerates all 2^24 inside
    configurations of the cell for one boundary condition, on the binding reduced cell records and on nrand
    fixed-seed random reduced records."""
    t = time.monotonic()
    rng = random.Random(20260927)
    chosen = {(bytes(lab), x, y) for lab, x, y, _ in res["argmax_u0"] + res["argmax_u_half"]}
    chosen |= {(canon.label(i), canon.xs[i], canon.ys[i]) for i in rng.sample(range(canon.n), nrand)}
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp, "jbrute")
        subprocess.run([cc, "-O2", "-o", str(exe), str(HERE / "engine" / "jbrute.c")], check=True)
        for lab, x, y in sorted(chosen):
            c = canon.counts(canon.find(lab + bytes((x, y))))
            want = [c[0], 2 * c[3], 2 * c[4], 2 * (c[1] + c[2]), 2 * c[7], 2 * c[8], 2 * (c[5] + c[6])]
            out = subprocess.run([str(exe), "cell", "".join("0123456789abcdef"[v] for v in lab), str(x), str(y)], check=True,
                                 capture_output=True, text=True).stdout.split()
            assert [int(v) for v in out] == want, (lab, x, y)
    report["brute_force_cells"] = {"records": len(chosen), "binding": len(chosen) - nrand, "random": nrand}
    step(f"literal brute force (engine/jbrute.c) equals {len(chosen)} reduced cell records "
         f"({len(chosen) - nrand} binding, {nrand} random)", t)


def check_engine(canons, report):
    t = time.monotonic()
    specs = HERE / "params" / "specs"
    eng = L.PathOneEngine(L.parse_spec(specs / "int.spec"))
    items = []
    for lab in L.gen_noncrossing(8):
        for (x, y), c in eng.records(lab).items():
            items.append(bytes(lab) + bytes((x, y)) + L.COUNTS.pack(*c))
    data = L.canonical_stream(items, 8)
    assert L.Canon(data, 8).sha256 == DIGESTS["int"]
    step(f"Python engine: intersection family recomputed, canonical digest equal ({len(items)} records)", t)
    summary = {"int": {"partitions": 1430, "records": len(items)}}
    for fam in ("side", "sideo"):
        t = time.monotonic()
        canon, res = canons[fam]
        nb, npos = L.FAMILIES[fam][:2]
        eng = L.PathOneEngine(L.parse_spec(specs / f"{fam}.spec"))
        labels = sorted({canon.label(i) for i in range(canon.n)})
        chosen = set(labels[::SUBSET_STRIDE]) | {lab for lab, _, _, _ in res["argmax"]}
        want = {}
        for i in range(canon.n):
            if canon.label(i) in chosen:
                want[(canon.label(i), canon.xs[i], canon.ys[i])] = canon.counts(i)
        got = {}
        for lab in chosen:
            for (x, y), c in eng.records(tuple(lab[:nb])).items():
                got[(lab, x, y)] = c
        assert got == want, fam
        summary[fam] = {"partitions": len(chosen), "records": len(got)}
        step(f"Python engine: {fam} records of {len(chosen)} partitions recomputed, all {len(got)} equal", t)
    report["python_engine"] = summary


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--with-cc", action="store_true",
                    help="also compile engine/jbrute.c (C compiler CC, default cc) and brute-force the binding and "
                         "some random reduced cell records (about 2 s each)")
    ap.add_argument("--brute-random", type=int, default=8, help="random reduced cell records for --with-cc")
    args = ap.parse_args()
    start = time.monotonic()
    cert_path = HERE / "params" / "cert.json"
    with open(cert_path) as fh:
        cert_json = json.load(fh)
    cert = L.load_cert(cert_path)
    report = {}
    check_constants(cert, report)
    check_geometry(report)
    check_special(cert_json, cert, report)
    canons = check_path1(cert, report)
    check_engine(canons, report)
    del canons
    canon, res = check_cells(cert, report)
    if args.with_cc:
        check_brute(canon, res, os.environ.get("CC", "cc"), args.brute_random, report)
    print(f"PASS window_d2 light check ({time.monotonic() - start:.1f}s)")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

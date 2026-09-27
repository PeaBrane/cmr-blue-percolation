"""Full reproduction of the d = 2 exact local-ratio certificate from scratch (C engine + standard library).

Steps:
 1. regenerate the engine spec files (geom.py) and compare them with params/specs;
 2. compile engine/lrm.c and recompute every record family: the Path-1 side, origin and intersection
    boxes (all non-crossing boundary partitions), the reduced cell list (the D4-orbit representatives
    listed in certificates/cell_reps.canon.xz) and the cell family on ALL 208,012 non-crossing
    partitions (3,879,876 records);
 3. put each dump in canonical form and compare its SHA-256 with the canonical digest recorded with the
    certificate (for the four committed families this is the SHA-256 of the committed file);
 4. check (C1) on every regenerated Path-1 record and (C2) on all 3,879,876 cell records exactly,
    completeness of the all-partition family, and that every one of its records equals the record of
    its D4 image in the reduced list (Lemma 6.15 of the diluted-model manuscript, checked record by record);
 5. compile the literal brute-force programs engine/jbrute.c and engine/jbrute1.c, which enumerate
    all inside configurations for one boundary condition at a time and share no code with the engine,
    and compare them with sampled records: the binding and tightest records of every inequality and
    fixed-seed random ones.

Needs a C compiler (CC, default cc) with pthreads and Python 3.10+ (standard library only). The cell
family needs about 1.3 GB of memory in the engine and about 1 GB in this script. Reference run: Linux,
gcc 13.3.0, CPython 3.13.5, 6 threads: see README.md for the measured times.
usage: python reproduce_window_d2.py [--threads N] [--workdir DIR] [--brute-samples N] [--skip-brute]
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
import time

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_window_d2 as light  # noqa: E402
import lrm_common as L  # noqa: E402

ALL_CELLS_DIGEST = "7a69009a9acc7e31ff09a42f2e663e0b7f8bf72e2fc13037f5e3b660dfbdcb32"
HEX = "0123456789abcdef"


def run(cmd, **kw):
    print("$", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def canon_of_dump(path, npos):
    data = L.canonical_stream(L.raw_items(str(path), npos), npos)
    return L.Canon(data, npos)


def brute_force(work, cc, canons, threads, nsamples, cert):
    """Literal per-boundary-condition brute force on sampled records."""
    run([cc, "-O2", "-o", work / "jbrute", HERE / "engine" / "jbrute.c"])
    run([cc, "-O2", "-o", work / "jbrute1", HERE / "engine" / "jbrute1.c"])
    rng = random.Random(20260927)
    jobs = []
    # cells: the binding records (u = 0 and u = 1/2) and random reduced records
    res = L.check_c2(canons["cell_reps"].rows(), cert)
    cell = {(bytes(lab), x, y): c for lab, x, y, c in res["argmax_u0"] + res["argmax_u_half"]}
    reps = canons["cell_reps"]
    for i in rng.sample(range(reps.n), nsamples):
        cell[(reps.label(i), reps.xs[i], reps.ys[i])] = reps.counts(i)
    for (lab, x, y), c in sorted(cell.items()):
        expect = [c[0], 2 * c[3], 2 * c[4], 2 * (c[1] + c[2]), 2 * c[7], 2 * c[8], 2 * (c[5] + c[6])]
        jobs.append(("cell", [work / "jbrute", "cell", "".join(HEX[t] for t in lab), x, y], expect))
    for fam in ("side", "sideo", "int"):
        canon = canons[fam]
        res = L.check_c1(canon.rows(), cert, fam)
        idx = {canon.find(bytes(lab) + bytes((x, y))) for lab, x, y, _ in res["argmax"]}
        idx |= set(rng.sample(range(canon.n), max(1, nsamples // 2)))
        for i in sorted(idx):
            c = canon.counts(i)
            jobs.append((fam, [work / "jbrute1", fam, "".join(HEX[t] for t in canon.label(i)), canon.xs[i],
                               canon.ys[i]], list(c[1:5])))

    def one(job):
        kind, cmd, expect = job
        out = subprocess.run([str(c) for c in cmd], check=True, capture_output=True, text=True).stdout.split()
        return kind, [int(t) for t in out] == expect

    with ThreadPoolExecutor(threads) as ex:
        results = list(ex.map(one, jobs))
    summary = {}
    for kind, ok in results:
        s = summary.setdefault(kind, [0, 0])
        s[0] += 1
        s[1] += ok
    for kind, (n, ok) in summary.items():
        assert n == ok, (kind, n, ok)
    return {k: v[0] for k, v in summary.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument("--workdir", default=str(HERE / "build"))
    ap.add_argument("--brute-samples", type=int, default=48)
    ap.add_argument("--skip-brute", action="store_true")
    args = ap.parse_args()
    cc = os.environ.get("CC", "cc")
    work = Path(args.workdir)
    (work / "specs").mkdir(parents=True, exist_ok=True)
    cert = L.load_cert(HERE / "params" / "cert.json")
    L.check_c0(cert)
    t0 = time.monotonic()
    report = {"compiler": subprocess.run([cc, "--version"], capture_output=True, text=True).stdout.splitlines()[0],
              "python": sys.version.split()[0], "threads": args.threads}

    saved, sys.argv = sys.argv, ["geom.py", str(work / "specs")]
    try:
        light.geom.main()
    finally:
        sys.argv = saved
    for name in ("cell", "side", "sideo", "int"):
        assert (work / "specs" / f"{name}.spec").read_bytes() == (HERE / "params" / "specs" / f"{name}.spec").read_bytes()
    run([cc, "-O2", "-Wall", "-Wextra", "-o", work / "lrm", HERE / "engine" / "lrm.c", "-lpthread"])
    lrm = work / "lrm"
    times = {}
    for fam in ("int", "side", "sideo"):
        t = time.monotonic()
        run([lrm, work / "specs" / f"{fam}.spec", "all", work / f"{fam}.dump", args.threads, 18])
        times[fam] = round(time.monotonic() - t, 1)
    committed = {fam: L.Canon.from_xz(HERE / "certificates" / f"{fam}.canon.xz", L.FAMILIES[fam][1])
                 for fam in ("int", "side", "sideo", "cell_reps")}
    reps = sorted({committed["cell_reps"].label(i) for i in range(committed["cell_reps"].n)} | {bytes(12)})
    assert len(reps) == 26584
    (work / "cell_reps.list").write_bytes(b"".join(reps))
    t = time.monotonic()
    run([lrm, work / "specs" / "cell.spec", f"list:{work / 'cell_reps.list'}", work / "cell_reps.dump",
         args.threads, 23])
    times["cell_reps"] = round(time.monotonic() - t, 1)
    t = time.monotonic()
    run([lrm, work / "specs" / "cell.spec", "all", work / "cell_all.dump", args.threads, 23])
    times["cell_all"] = round(time.monotonic() - t, 1)
    report["engine_seconds"] = times

    canons = {}
    digests = {}
    for fam in ("int", "side", "sideo", "cell_reps"):
        canon = canon_of_dump(work / f"{fam}.dump", L.FAMILIES[fam][1])
        assert canon.sha256 == light.DIGESTS[fam] == committed[fam].sha256, fam
        digests[fam] = canon.sha256
        canons[fam] = canon
        if fam != "cell_reps":
            nb, _, origin, _ = L.FAMILIES[fam]
            assert L.completeness(canon, nb, origin)[0]
            res = L.check_c1(canon.rows(), cert, fam)
            assert res["violations"] == 0 and res["nonzero_cell_slots"] == 0
    print("PASS Path-1 families and reduced cell list regenerated; canonical digests equal the recorded ones",
          flush=True)

    # the all-partition cell family
    t = time.monotonic()
    items = sorted(L.raw_items(str(work / "cell_all.dump"), 12))
    h = hashlib.sha256()
    h.update(b"".join(r[:12] for r in items))
    h.update(bytes(r[12] for r in items))
    h.update(bytes(r[13] for r in items))
    for k in range(0, len(items), 65536):
        h.update(b"".join(r[14:] for r in items[k:k + 65536]))
    digests["cell_all"] = h.hexdigest()
    assert digests["cell_all"] == ALL_CELLS_DIGEST
    assert len(items) == 3879876 and all(items[k][:14] < items[k + 1][:14] for k in range(len(items) - 1))
    rows = ((r[:12], r[12], r[13], L.COUNTS.unpack_from(r, 14)) for r in items)
    res = L.check_c2(rows, cert)
    assert res["records"] == 3879876 and res["violations"] == 0 and res["D0_with_star_pivot"] == 0
    for arg in ("argmax_u0", "argmax_u_half"):
        by_d = {}
        for _, _, _, c in res[arg]:
            by_d[c[0]] = by_d.get(c[0], 0) + 1
        assert by_d == {276: 8, 552: 8, 828: 4}, by_d
    ncs = L.gen_noncrossing(12)
    expected = {bytes(p): L.comb(max(p) + 1, 2) for p in ncs if max(p) > 0}
    seen = {}
    grid_d0 = 0
    perms = L.d4_permutations()
    reduced = canons["cell_reps"]
    rep_counts = {reduced.key(i): c for i, c in enumerate(L.COUNTS.iter_unpack(reduced.cs))}
    rep_set = {reduced.label(i) for i in range(reduced.n)}
    to_rep = {}
    mismatches = 0
    for r in items:
        lab, x, y = r[:12], r[12], r[13]
        seen[lab] = seen.get(lab, 0) + 1
        c = L.COUNTS.unpack_from(r, 14)
        if c[0] == 0 and c[1] + c[2] + c[5] + c[6] > 0:
            grid_d0 += 1
        m = to_rep.get(lab)
        if m is None:
            for g in perms:
                img, bmap = L.image(tuple(lab), g)
                if bytes(img) in rep_set:
                    m = to_rep[lab] = (bytes(img), bmap)
                    break
        img, bmap = m
        a, b = sorted((bmap[x], bmap[y]))
        if rep_counts[img + bytes((a, b))] != c:
            mismatches += 1
    assert seen == expected and len(seen) == 208011
    assert grid_d0 == 69292 and mismatches == 0
    report["cell_all"] = {"records": len(items), "canonical_sha256": digests["cell_all"],
                          "max_ratio_u0": L.dec(res["max_ratio_u0"], 9, True),
                          "max_ratio_u_half": L.dec(res["max_ratio_u_half"], 9, True),
                          "binding_boundary_conditions": len(res["argmax_u0"]),
                          "grid_pivot_with_D0": grid_d0, "D4_covariance_mismatches": mismatches,
                          "python_seconds": round(time.monotonic() - t, 1)}
    del items
    print("PASS all 3,879,876 cell records: canonical digest, (C2) exactly, completeness, D4 covariance with the "
          "reduced list", flush=True)
    report["canonical_sha256"] = digests
    if not args.skip_brute:
        t = time.monotonic()
        report["brute_force_samples"] = brute_force(work, cc, canons, args.threads, args.brute_samples, cert)
        report["brute_force_seconds"] = round(time.monotonic() - t, 1)
        print("PASS literal brute force equals every sampled record", flush=True)
    report["total_seconds"] = round(time.monotonic() - t0, 1)
    print(f"PASS window_d2 full reproduction ({report['total_seconds']}s)")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

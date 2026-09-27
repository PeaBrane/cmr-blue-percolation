"""FULL reproduction of the macrostep engine certificates from scratch (needs NumPy and Numba; hours of CPU at most).

For every selected run of params.json this program
  1. runs the engine, the a priori error program of Section 5.5 of the first manuscript (runs with code_base "R0":
     engine/ind_certify.py; "R0t": engine/ind_certify_t.py, whose forward-tail caches engine/tail_tight.py builds
     first if they are missing) from the frozen exact inputs;
  2. converts the output to the portable format (engine/export_certificate.py);
  3. compares the SHA-256 digest of the recomputed near data (etabar, Mbar, the near Green rows and eta_F) with the
     digest stored in certificates/<label>.json;
  4. re-decides the recomputed certificate exactly twice: checks/recheck_certs.py on the engine pickle (an
     independent exact rechecker) and certify_macrostep.py --cert-dir on the exported certificate (standard library).
With --l2 it also reruns the directed-rounding program (engine_l2/) at the instances of the selected runs of
engine/ind_certify.py and compares it entrywise with the fresh output of the engine.

Digest agreement means that the stored matrices are recomputed bit for bit. The test vectors w and the scalars
lambda, H come from a floating-point search (LAPACK solves and a bisection); they are recorded but not required to
agree, because any vector that passes the exact checks gives a valid certificate.

    python reproduce_engine.py [--runs theorems|all|LABEL,...] [--out DIR] [--nproc N] [--tails DIR] [--l2]

Default: the seven near sets of Table 5 of the first manuscript and the two dimension-7 runs with tightened forward
tails (THEOREM_RUNS below), with at most 6 worker processes.
Measured on a 24-core x86-64 Linux machine with 6 worker processes: see README.md.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

HERE = Path(__file__).resolve().parent
ENGINE = HERE / "engine"
ENGINE_L2 = HERE / "engine_l2"
PARAMS = json.loads((HERE / "params.json").read_text())
THEOREM_RUNS = ["R0-d9-246", "R0-d8-246", "R0-d8-82", "R0-d8-46", "R0-d8s-46", "R0-d7a-804", "R0-d7b-804",
                "R0t-d7a-246", "R0t-d7b-246"]
# point keys of the directed-rounding program (engine_l2/common.py POINTS) of the instances
L2_POINT = {"inst1-d9-t7/50": "7/50", "inst1-d8-t3/20": "3/20", "inst2-d8-t3/20": "nc-3/20",
            "inst2-d8-t31/200": "nc-31/200", "inst3-d7-t3/20": "nc-3/20", "inst3-d7-t31/200": "nc-31/200"}
L2_ZMAX = {"1/10": 22, "3/25": 28}


def env(nproc):
    e = dict(os.environ)
    e.update(NPROC=str(nproc), NUMBA_NUM_THREADS="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
             MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1", RS=str(PARAMS["engine"]["R_S"]),
             R3=str(PARAMS["engine"]["R3"]), NF=str(PARAMS["engine"]["N_F"]))
    return e


def log(msg, fh):
    print(msg, flush=True)
    fh.write(msg + "\n")
    fh.flush()


def start_tails(fs, tails, fh):
    """Start building the forward-tail caches tt_f{f}.pkl of engine/ind_certify_t.py (one background process per f)."""
    nex = PARAMS["engine"]["forward_tails"]["R0t"]["NEX"]
    procs = []
    for f in sorted(fs):
        out = tails / f"tt_f{f}.pkl"
        if out.exists():
            log(f"tail cache {out} exists; reused", fh)
            continue
        cmd = [sys.executable, str(ENGINE / "tail_tight.py"), str(f), str(nex[str(f)]), str(out)]
        procs.append((f, subprocess.Popen(cmd, stdout=open(tails / f"tt_f{f}.log", "w"), stderr=subprocess.STDOUT,
                                          env=env(1))))
    return procs


def wait_tails(procs, tails, fh):
    for f, pr in procs:
        assert pr.wait() == 0, f"tail_tight.py f={f} failed"
        log((tails / f"tt_f{f}.log").read_text().strip().splitlines()[-1], fh)


def run_engine(run, out, tails, nproc, fh):
    I = PARAMS["instances"][run["instance"]]
    fam = PARAMS["engine"]["family"][str(I["d"])]
    script = "ind_certify.py" if run["code_base"] == "R0" else "ind_certify_t.py"
    pkl = out / "engine" / f"{run['label']}.pkl"
    args = [str(I["d"]), I["p"], I["g_K"], I["g_h"], I["xbar"], str(fam["f"]), str(PARAMS["engine"]["c"]),
            fam["y"], str(run["R_A"]), str(run["R_D"]), str(PARAMS["engine"]["N0"]), str(pkl)]
    e = env(nproc)
    if run["code_base"] == "R0t":
        e["TT_DIR"] = str(tails)
    t0 = time.time()
    with open(out / "engine" / f"{run['label']}.log", "w") as lf:
        rc = subprocess.run([sys.executable, str(ENGINE / script), *args], stdout=lf, stderr=subprocess.STDOUT, env=e)
    assert rc.returncode == 0, f"engine failed for {run['label']}"
    return pkl, time.time() - t0


def export(run, pkl, out):
    stem = out / "export" / run["label"]
    subprocess.run([sys.executable, str(ENGINE / "export_certificate.py"), str(pkl), str(stem), "--label", run["label"],
                    "--code-base", run["code_base"], "--instance", run["instance"]], check=True,
                   stdout=subprocess.DEVNULL, env=env(1))
    return json.loads(Path(str(stem) + ".json").read_text())


def run_l2(run, out, nproc, fh):
    """Directed-rounding program at the run's instance and near set; returns the path of its output pickle."""
    I = PARAMS["instances"][run["instance"]]
    fam = PARAMS["engine"]["family"][str(I["d"])]
    d2 = out / "l2"
    cache = d2 / f"lat_c3_y{fam['y'].replace('/', '_')}_N{PARAMS['engine']['N0']}.pkl"
    e = env(nproc)
    if not cache.exists():
        subprocess.run([sys.executable, str(ENGINE_L2 / "lateral_cache.py"), "3", fam["y"], str(PARAMS["engine"]["N0"]),
                        str(L2_ZMAX[fam["y"]]), str(cache)], check=True, cwd=ENGINE_L2, env=e,
                       stdout=subprocess.DEVNULL)
    e["LATCACHE"] = str(cache)
    pkl = d2 / f"{run['label']}.L2.pkl"
    with open(d2 / f"{run['label']}.L2.log", "w") as lf:
        subprocess.run([sys.executable, str(ENGINE_L2 / "certify.py"), str(I["d"]), L2_POINT[run["instance"]],
                        str(fam["f"]), "3", fam["y"], str(run["R_A"]), str(PARAMS["engine"]["N0"]), str(pkl),
                        str(run["R_D"])], check=True, cwd=ENGINE_L2, env=e, stdout=lf, stderr=subprocess.STDOUT)
        subprocess.run([sys.executable, str(ENGINE_L2 / "certify2.py"), str(pkl)], check=True, cwd=ENGINE_L2, env=e,
                       stdout=lf, stderr=subprocess.STDOUT)
    return pkl


def compare_l2(l2pkl, r0pkl):
    """Entrywise comparison of the directed-rounding program (keys "L2") and the engine (keys "R0") on the same
    near set (both upper bounds of the same quantities)."""
    import pickle
    import numpy as np
    L = pickle.load(open(l2pkl, "rb"))
    Lc = pickle.load(open(str(l2pkl) + ".cert.pkl", "rb"))
    R = pickle.load(open(r0pkl, "rb"))
    key = lambda s: (tuple(int(v) for v in s[0]), tuple(int(v) for v in s[1]))
    Ls, Rs = [key(s) for s in L["states"]], [key(s) for s in R["states"]]
    assert set(Ls) == set(Rs), "different near sets"
    perm = [Ls.index(s) for s in Rs]
    ML = np.array(L["M"])[np.ix_(perm, perm)]
    mask = ML > 1e-14
    GL = np.array([L["Gn"][s] for s in L["states"]])[np.ix_(perm, perm)]
    rel = lambda a, b: [float((a / b - 1).min()), float((a / b - 1).max())]
    return {"eta_R0_over_L2_minus_1": rel(np.array(R["eta"]), np.array(L["eta_bar"])[perm]),
            "M_R0_over_L2_minus_1": rel(np.array(R["M"])[mask], ML[mask]),
            "Gnear_R0_over_L2_minus_1": rel(np.array(R["Gnear"]), GL),
            "etaF_L2": float(Lc["etaF"]), "etaF_R0": float(R["far"]["etaF"]),
            "lambda_min_L2": float(Lc["lam"]), "lambda_best_L2": float(Lc["lam_best"]),
            "theta_best_L2": float(Lc["theta_best"]), "L2_exact_recheck": list(Lc["exact_ok"])}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", default="theorems", help="'theorems' (default), 'all', or comma-separated labels")
    ap.add_argument("--out", default=None, help="output directory (default: a new temporary directory)")
    ap.add_argument("--nproc", type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument("--tails", default=None, help="directory of the tail caches of engine/ind_certify_t.py (default: OUT/tails)")
    ap.add_argument("--l2", action="store_true", help="also rerun the directed-rounding program (engine_l2/) for the selected runs of engine/ind_certify.py")
    ap.add_argument("--no-recheck", action="store_true", help="skip the final exact re-check")
    a = ap.parse_args()
    runs = {r["label"]: r for r in PARAMS["runs"]}
    labels = (THEOREM_RUNS if a.runs == "theorems" else list(runs) if a.runs == "all" else a.runs.split(","))
    out = Path(a.out or tempfile.mkdtemp(prefix="macrostep-repro-")).resolve()
    for sub in ("engine", "export", "l2"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    tails = Path(a.tails).resolve() if a.tails else out / "tails"
    tails.mkdir(parents=True, exist_ok=True)
    fh = open(out / "reproduction.log", "a")
    log(f"reproduce_engine: {len(labels)} runs, nproc={a.nproc}, out={out}, python {sys.version.split()[0]}", fh)
    T0 = time.time()
    fs = {PARAMS["engine"]["family"][str(PARAMS["instances"][runs[l]["instance"]]["d"])]["f"]
          for l in labels if runs[l]["code_base"] == "R0t"}
    tail_procs = start_tails(fs, tails, fh)            # built in the background while the other runs proceed
    nproc_r0 = max(1, a.nproc - len(tail_procs))
    summary = []
    for lab in labels:
        run = runs[lab]
        if run["code_base"] == "R0t" and tail_procs:
            t0 = time.time()
            wait_tails(tail_procs, tails, fh)
            log(f"tail caches ready (waited {time.time() - t0:.0f}s)", fh)
            tail_procs, nproc_r0 = [], a.nproc
        pkl, dt = run_engine(run, out, tails, nproc_r0, fh)
        h = export(run, pkl, out)
        rc = subprocess.run([sys.executable, str(HERE / "checks" / "recheck_certs.py"), str(pkl)],
                            stdout=open(out / "engine" / f"{lab}.recheck.log", "w"), stderr=subprocess.STDOUT,
                            env=env(1)).returncode
        stored = json.loads((HERE / "certificates" / f"{lab}.json").read_text())
        same_near = h["near_data_sha256"] == stored["near_data_sha256"]
        same_all = h["payload"]["sha256_uncompressed"] == stored["payload"]["sha256_uncompressed"]
        entry = {"label": lab, "engine_seconds": round(dt, 1), "near_data_sha256": h["near_data_sha256"],
                 "near_data_equal_to_stored": same_near, "payload_equal_to_stored": same_all,
                 "payload_sha256": h["payload"]["sha256_uncompressed"], "recheck_certs_pass": rc == 0,
                 "lambda_min_fresh": float.fromhex(h["certificates"]["min_lambda"]["lambda"]),
                 "lambda_min_stored": float.fromhex(stored["certificates"]["min_lambda"]["lambda"])}
        if a.l2 and run["code_base"] == "R0":
            t0 = time.time()
            l2 = run_l2(run, out, nproc_r0, fh)
            entry["L2"] = compare_l2(l2, pkl)
            entry["L2"]["seconds"] = round(time.time() - t0, 1)
        summary.append(entry)
        log(f"{lab}: engine {dt:.0f}s; near data {'IDENTICAL to' if same_near else 'DIFFERENT from'} the stored "
            f"certificate ({h['near_data_sha256'][:16]}); whole payload {'identical' if same_all else 'differs'}; "
            f"recheck_certs {'PASS' if rc == 0 else 'FAIL'}"
            + (f"; directed-rounding lambda {entry['L2']['lambda_min_L2']:.6f}" if 'L2' in entry else ""), fh)
    ok = all(e["near_data_equal_to_stored"] and e["recheck_certs_pass"] for e in summary)
    if not a.no_recheck:
        rc = subprocess.run([sys.executable, str(HERE / "certify_macrostep.py"), "--cert-dir", str(out / "export"),
                             "--runs", ",".join(labels), "--certificates-only"],
                            stdout=open(out / "recheck.log", "w"), stderr=subprocess.STDOUT)
        log(f"exact re-check of the recomputed certificates: {'PASS' if rc.returncode == 0 else 'FAIL'} "
            f"(recheck.log)", fh)
        ok = ok and rc.returncode == 0
    (out / "reproduction.json").write_text(json.dumps({"runs": summary, "seconds": round(time.time() - T0, 1)},
                                                      indent=1) + "\n")
    log(f"{'PASS' if ok else 'FAIL'}: {sum(e['near_data_equal_to_stored'] for e in summary)}/{len(summary)} runs "
        f"reproduce the stored near data bit for bit [{time.time() - T0:.0f}s]", fh)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

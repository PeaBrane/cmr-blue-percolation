"""Convert an engine output pickle (ind_certify.py or ind_certify_t.py) into the portable certificate format.

    python export_certificate.py RUN.pkl OUT_STEM --label LABEL --code-base R0|R0t --instance KEY

writes OUT_STEM.json (header) and OUT_STEM.bin.xz (payload). The payload is the concatenation of IEEE-754
binary64 little-endian arrays, compressed with LZMA (xz container):

    eta[n], w_min_lambda[n], w_best_theta[n], M[n*n] (row-major), Gnear[n*n] (row-major)

Every binary64 number is an exact dyadic rational, so the payload is an exact record of the stored certificate.
Scalars (etaF, lambda, H) are stored in the header as float.hex() strings, which are exact as well. The header
also records SHA-256 digests of each array and of the "near data" (eta, M, Gnear, etaF), which the full
reproduction compares; the test vectors w, lambda and H come from a floating-point search and may legitimately
differ between platforms, while any vector that passes the exact checks gives a valid certificate.

Needs NumPy (to read the pickle); the portable files are read with the standard library alone.
"""
import argparse
import hashlib
import json
import lzma
import pickle
import sys

import numpy as np


def le_bytes(a):
    return np.ascontiguousarray(np.asarray(a, dtype="<f8")).tobytes()


def hexf(x):
    return float(x).hex()


def near_digest(eta, M, Gnear, etaF):
    h = hashlib.sha256()
    for part in (le_bytes(eta), le_bytes(M), le_bytes(Gnear), hexf(etaF).encode()):
        h.update(part)
    return h.hexdigest()


def export(pkl, stem, label, code_base, instance):
    o = pickle.load(open(pkl, "rb"))
    cert = o["cert"]
    if not cert.get("ok"):
        raise SystemExit(f"{pkl}: no certificate stored")
    n = len(o["eta"])
    arrays = [("eta", o["eta"], n), ("w_min_lambda", cert["w"], n), ("w_best_theta", cert["w_th"], n),
              ("M", o["M"], n * n), ("Gnear", o["Gnear"], n * n)]
    blobs = []
    layout = []
    for name, arr, size in arrays:
        b = le_bytes(np.asarray(arr).reshape(-1))
        assert len(b) == 8 * size, name
        blobs.append(b)
        layout.append({"name": name, "count": size, "sha256": hashlib.sha256(b).hexdigest()})
    payload = b"".join(blobs)
    comp = lzma.compress(payload, preset=9 | lzma.PRESET_EXTREME)
    with open(stem + ".bin.xz", "wb") as fh:
        fh.write(comp)
    far = o["far"]
    header = {
        "format": "macrostep-certificate/1",
        "label": label,
        "code_base": code_base,
        "instance": instance,
        "inputs": {"d": o["d"], "p": str(o["p"]), "g_K": str(o["gK"]), "g_h": str(o["gh"]), "xbar": str(o["xbar"]),
                   "f": o["f"], "c": o["c"], "y": str(o["y"]), "R_A": o["RA"], "R_D": o["RD"], "N0": o["N0"]},
        "engine_local_constants": {"rho_minus": str(o["rho"]), "rho_c": str(o["rho_c"]), "kappa": str(o["kappa"]),
                                   "t": str(o["t"])},
        "n_states": n,
        "states": [[list(A), list(D)] for (A, D) in o["states"]],
        "etaF": hexf(far["etaF"]),
        "certificates": {
            "min_lambda": {"lambda": hexf(cert["lam"]), "H": hexf(cert["H"]), "w": "w_min_lambda"},
            "best_theta": {"lambda": hexf(cert["lam_th"]), "H": hexf(cert["H_th"]), "w": "w_best_theta"},
        },
        "payload": {"file": stem.rsplit("/", 1)[-1] + ".bin.xz", "encoding": "IEEE-754 binary64 little-endian, LZMA/xz",
                    "bytes_uncompressed": len(payload), "bytes_compressed": len(comp),
                    "sha256_uncompressed": hashlib.sha256(payload).hexdigest(), "layout": layout},
        "near_data_sha256": near_digest(o["eta"], o["M"], o["Gnear"], far["etaF"]),
        "engine_report": {  # floating summaries printed by the engine; informational, never used by the checker
            "rho_Mbar_float": float(cert["rhoM"]), "kappa_maxB": float(o["xmax"]), "h_in": float(far["h_in"]),
            "eta_out": float(far["eta_out"]), "Gamma_w_min_lambda": float(cert["Gam"]), "C_fin": float(cert["Cfin"]),
            "theta_min_lambda": float(cert["theta"]), "theta_best": float(cert["theta_best"]),
            "bound_best": float(cert["bound_best"])},
    }
    with open(stem + ".json", "w") as fh:
        json.dump(header, fh, indent=1)
        fh.write("\n")
    return header


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pkl")
    ap.add_argument("stem")
    ap.add_argument("--label", required=True)
    ap.add_argument("--code-base", required=True, choices=["R0", "R0t"])
    ap.add_argument("--instance", required=True)
    a = ap.parse_args()
    h = export(a.pkl, a.stem, a.label, a.code_base, a.instance)
    print(f"{a.label}: n={h['n_states']} payload {h['payload']['bytes_uncompressed']} -> "
          f"{h['payload']['bytes_compressed']} bytes; near data sha256 {h['near_data_sha256']}")


if __name__ == "__main__":
    sys.exit(main())

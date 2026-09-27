# Macrostep interaction-matrix engine (Theorems Q9, Q8, Q8♯, Q7)

This directory holds the certificates for the macrostep route in dimensions 9, 8 and 7. The route replaces the
global step of the overlap-revealed route (uniform oriented paths) with macrostep paths that make monotone lateral
excursions. Theorem M reduces the torus bound (★) `P^aux_L(E_n) >= theta_*` to an exact finite criterion
(Theorem G). The criterion has two conditions on a near set `C'` of orbit types `z = (A, D)`:

```text
(i)   (Mbar w)(z) + H eta_F etabar(z) <= lambda w(z)   for every z in C',
(ii)  max_z (Gnear w)(z) + H eta_F      <= lambda H,
      w > 0,  H > 0,  0 < eta_F < lambda < 1,
theta_* = (99/100) / (C_fin (1 + Lambda lambda H / (1 - lambda))),   Lambda = max(max_z etabar/w, 1/H).
```

`Mbar`, `etabar`, `Gnear` and `eta_F` are rigorous upper bounds computed by the engine. The vector `w` and the
scalars `lambda`, `H` come from a floating-point search. Each binary64 number is read as the exact dyadic rational
it represents.

Section references (§) below are to the research write-up of the macrostep route, which the manuscript
supersedes: the statements are in §2, the computations in §8, the tables in §9 and the scripts in §13.

## Theorems and certificates

| Theorem | (d, tanh β) | local inputs | certificate(s) used | λ ≤ (minimal) | θ\* ≥ (θ\*-maximising) |
|---|---|---|---|---|---|
| Q9 | (9, 7/50) | instance 1: tangent bond floor and Holley line of the overlap-revealed route | `R0-d9-246` | 0.779209 | 0.001908 |
| Q8 | (8, 3/20) | instance 1 | `R0-d8-246` (also `R0-d8-82`, `R0-d8-46`) | 0.892491 | 0.0006603 |
| Q8♯ | (8, 3/20) | instance 2: Lemmas N, P, S (noisy cavity) | `R0-d8s-46` | 0.858128 | 0.001937 |
| Q7 | (7, 31/200) | instance 3: Lemmas N, P, S | `R0-d7b-804` | 0.924670 | 0.0003267 |
| Q7 | (7, 3/20) | instance 3 | `R0-d7a-804` | 0.935275 | 0.0002274 |
| Q7 addendum | (7, 3/20), (7, 31/200) | instance 3 and Lemma R.1 (tightened forward tails) | `R0t-d7a-246`, `R0t-d7b-246` | 0.966950, 0.964284 | 0.0002132, 0.0002888 |

The values in the last two columns are printed by `certify_macrostep.py`. They are exact values rounded in the safe
direction, with λ rounded up and θ\* rounded down. The theorems quote some of them less sharply; for example Q9
quotes `lambda <= 0.779215`, which is the value of code base L2. The other 16 stored certificates are the remaining rows of
§9.2 and Table 9.5b (R0 and R0t at 46, 82, 246 and 804 types). No theorem uses them, but every one is re-decided.

## Files

| Path | Content |
|---|---|
| `certify_macrostep.py` | **Default checker** (standard library). |
| `params.json` | Frozen exact inputs: the six instances `(d, t, p, g_K, g_h, xbar)`, the instance-1 tangent points and Holley weights, the engine constants (`c = 3`, `N0 = 24`, `N_F = 60`, `R_S = 4`, families), the 25 runs, the theorem map, and 276 displayed numbers with their sources in the research write-up. |
| `local_sharpened.json` | Instance 2 and 3 local certificate data: `wc'`, `psi_U(x)`, Holley tangent points and `Lambda'_k`. It is identical, field by field, to the research pickles `paircert_*.pkl` and `holley_*.pkl` of code base L1. |
| `certificates/<label>.json`, `<label>.bin.xz` | The 25 stored certificates in portable form (format below). 23 MB in total; the largest file is 4.7 MB. |
| `engine/` | Code base R0 (`ind_tables.py`, `ind_family.py`, `ind_kernel.py`, `ind_green.py`, `ind_certify.py`), code base R0t (`ind_certify_t.py` and the forward-tail cache builder `tail_tight.py`), code base B (`indep_global.py`: local constants, standard library), and `export_certificate.py` (pickle to portable format). |
| `engine_l2/` | The first implementation (code base L2), copied unchanged from the research code apart from the path of code base B. `check_boost.py` and `check_green.py` need SciPy. |
| `reproduce_engine.py` | **Full reproduction** (NumPy and Numba). `reproduction_record.json` is the summary it wrote in the recorded run below. |
| `checks/recheck_certs.py` | The research exact rechecker. It re-decides both certificates of an engine pickle in exact integer arithmetic and recomputes `C_fin` and θ\*. It shares no code with the producers. `reproduce_engine.py` runs it on every fresh pickle. |
| `checks/check_configs_exact.py` | Optional (NumPy, Numba). Evaluates the defining formula of the per-macrostep weight in exact arithmetic on random configurations and checks that both kernels (R0 and L2) are upper bounds (§8.3 item 4). With the three argument lists in its docstring it gives the ratios in the research logs, for example R0/exact in `[1.000000000001, 1.000000732859]` over 8640 configurations at d = 9. The runs took 80 s, 57 s and 36 s (d = 9, 8, 7) on the Linux machine and used about 210 MB each. |
| `checks/cw_bounds.py` | Optional (NumPy). Exact Collatz–Wielandt enclosure of `rho(Mbar)` of a computed near matrix. It is used only for the robustness statements of §9.4 and Table 9.5c. |
| `search/freeze_instance1.py` | Optional (SciPy). Repeats the floating-point choice of the instance-1 tangent points and weights, with the same Nelder–Mead procedure as the research code base B. It is not part of any verification. The frozen values are reproduced exactly on Apple Silicon (SciPy 1.16.0, NumPy 2.2.6). On x86-64 Linux (SciPy 1.18.1, NumPy 2.5.3) the search ends at a slightly different point, and the tangent points differ from about the eighth significant digit. Any positive tangent points are admissible, and the frozen ones are the ones certified. |

Instances 2 and 3 are noisy-cavity instances. Their local certificates use the deciders in
`../overlap_revealed/noisy_cavity/local_certificates.py`. Instance 1 uses the functions of
`../overlap_revealed/certify.py`.

## Default verification

```sh
python macrostep/certify_macrostep.py
```

This needs only the standard library and makes every decision in exact rational or integer arithmetic.

1. **Local inputs.** Instance 1: the bond floor `p <= p_B` (all vertex sums `Delta_j <= 0`) and the Holley line
   in all `2d+1` environments. Instances 2 and 3: Lemma P (point and chord conditions, class sums `Psi(j) <= 0`)
   and Lemma S (every two-block class `Gamma >= 0`), together with the noisy-class constraint on `wc'`.
2. **Hypotheses of Theorem M.** `g_h < 1 <= g_K`, `g_K g_h <= 1`, `g_K^(2d) < 2718/1000`, the mean-field root
   test at the frozen `xbar`, `alpha < 1`, `q_* < 1`, `rho_-^2 >= c_cov/4` and `rho_c p < rho_-`. It also
   recomputes `rho_-`, `rho_c`, `kappa`, `t`, `b(1)`, `psi_3(0,0)` and `C_fin`, and the torus sizes `L_0(n)` at
   `n = 1, 10, 100` from rational logarithm enclosures.
3. **Stored certificates.** It checks the payload digests and that the certificate inputs and the engine's local
   constants equal those recomputed in step 2. It checks that the state list is exactly the orbit-type set
   `{|D|_1 <= 2R_D, |A|_1 <= R_A}` with `R_A >= 2c`. It re-decides (i) and (ii) for both certificates of each run
   (minimal λ and θ\*-maximising), and computes Λ and θ\* exactly.
4. **Displays.** It asserts every one of the 276 numbers in `params.json["displays"]`.

The JSON printed at the end equals `../expected/macrostep.json`. The elapsed times appear only on the `PASS`
lines. The JSON records two digests of the exact values: `local_vector_sha256` and `certificate_vector_sha256`.
With CPython 3.13.5 on an Apple Silicon laptop the run takes about 11 s and uses about 105 MB of memory. The
804-type payloads dominate both.

Options: `--runs LABEL,...` checks a subset; `--cert-dir DIR --certificates-only` re-decides freshly exported
certificates, skipping the local inputs and the displays.

**What the default run does not check.** It takes `Mbar`, `etabar`, `Gnear` and `eta_F` as stored. The
certificate is valid only if these are upper bounds of the quantities in Theorem G. That rests on the engine code
and its a-priori floating-point error analysis (§8.2 of the research write-up). The full reproduction below recomputes them, and
`engine_l2/` gives an independent second computation.

## Certificate format

`<label>.json` is the header. It contains `format = "macrostep-certificate/1"`, the run label, code base and
instance; the exact `inputs` (`d, p, g_K, g_h, xbar, f, c, y, R_A, R_D, N0`); the engine's local constants
(`rho_minus`, `rho_c`, `kappa`, `t`) as exact rationals; the ordered state list; `etaF`; and for each of
`min_lambda` and `best_theta` the scalars `lambda` and `H` as `float.hex` strings together with the name of its
`w` array. It also holds the payload description and SHA-256 digests. `engine_report` contains floating-point
diagnostics for information only.

`<label>.bin.xz` is an LZMA/xz stream of IEEE-754 binary64 little-endian arrays, concatenated in this order:
`eta[n]`, `w_min_lambda[n]`, `w_best_theta[n]`, `M[n*n]` and `Gnear[n*n]` (row-major). `near_data_sha256` hashes
`eta`, `M`, `Gnear` and `hex(etaF)`, which are the engine outputs. The test vectors `w` and the scalars `lambda`
and `H` come from a float search. They may differ between platforms, and any values that pass (i) and (ii) give a
valid certificate.

## Full reproduction

```sh
uv pip install --python .venv/bin/python numpy numba      # the reproduction used numpy 2.5.3, numba 0.67.0
.venv/bin/python macrostep/reproduce_engine.py --nproc 6 --l2 --out repro-out
```

For each selected run the script:

- runs the engine from `params.json` (R0: `engine/ind_certify.py`; R0t: `engine/ind_certify_t.py`, after
  building the exact forward-tail caches with `engine/tail_tight.py`);
- exports the output to the portable format and compares the digests with `certificates/`;
- re-decides the fresh certificate with `checks/recheck_certs.py` (from the pickle) and with
  `certify_macrostep.py --cert-dir` (from the export).

With `--l2` it also reruns the L2 pipeline at the same instance and near set and compares the two implementations
entry by entry. By default it runs the nine theorem runs above; `--runs all` runs all 25.

**Recorded reproduction** (27 September 2026, shared 24-core x86-64 Linux machine, CPython 3.13.5, numpy
2.5.3, numba 0.67.0, `--nproc 6 --l2`). The run started from a clean copy of this directory placed over the
repository at commit `97b57ed`. All nine theorem runs reproduced the stored near data **bit for bit**. The whole
payload, test vectors included, was identical too. Both exact rechecks (`checks/recheck_certs.py` and
`certify_macrostep.py --cert-dir`) passed on every fresh certificate. The script ended with
`PASS: 9/9 runs reproduce the stored near data bit for bit [1262s]`. The digests agree with `certificates/*.json`:

| run | engine time | `near_data_sha256` (prefix) | payload `sha256_uncompressed` (prefix) | L2 rerun |
|---|---|---|---|---|
| R0-d9-246 | 38 s | `160006a5add4baa4` | `e1cbaeb8e6fec1c9` | 156 s |
| R0-d8-246 | 33 s | `950eee026fc1fa49` | `760b580de3cabf03` | 92 s |
| R0-d8-82 | 24 s | `0b8ada9342c6eae4` | `12049b7a501c6860` | 33 s |
| R0-d8-46 | 18 s | `bbc102843ce532f4` | `d85f8e5bf65c5167` | 22 s |
| R0-d8s-46 | 18 s | `e8249c52f76d1039` | `6f31c47440cee1ee` | 22 s |
| R0-d7a-804 | 88 s | `05209e95fe5b2efb` | `35a76e01def3a64b` | 290 s |
| R0-d7b-804 | 81 s | `876d293165a332a7` | `503b9d011d24d353` | 250 s |
| R0t-d7a-246 | 38 s | `2695d35195949c7b` | `82ec3b074c4128d8` | – |
| R0t-d7b-246 | 41 s | `7b0672f87279a1c9` | `464e6d6a34b18856` | – |

The whole reproduction took 21 min of wall time (1262 s). Most of it went to the L2 reruns and to the exact
`f = 4` tail cache (`V_4(n)` for `n <= 6000`, 498 s on one core, built in parallel with the R0 runs). An earlier
identical run took 1199 s. Without `--l2` it should take about 10 min; that estimate comes from these timings
and was not measured separately. At every instance L2 agrees entrywise with R0, with R0 slightly larger as its
error margins require:

- `etabar`: `R0/L2 - 1` lies in `[2.0e-8, 3.9e-8]`;
- `Mbar`: `[2.2e-8, 5.5e-8]` on the entries above `1e-14`;
- near Green rows: `[1.09e-8, 1.10e-8]`.

The L2 exact recheck passed at every instance. Its minimal λ were 0.779214, 0.892511, 0.934026, 0.962153,
0.858693, 0.935275 and 0.924671.

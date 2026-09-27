# Exact computational certificates for CMR blue percolation

Companion code for Yan Ru Pei, *Towards the spin-glass transition in finite
dimensions via blue percolation* (2026), and for the companion manuscript on
the diluted ±1 Edwards–Anderson model, *Exact frustration cancellation and
uniqueness at all positive temperatures for diluted ±J spin glasses slightly
above the percolation threshold* (2026).

The revised manuscript (v2) proves CMR blue percolation of both overlap signs,
for iid fair unit couplings, in these cases:

- **dimension 10**, at `tanh(beta) = 27/200` and `13/100`;
- **dimension 11**, at `13/100` and `3/25`;
- **dimension 12**, at `3/25`, `11/100`, `23/200`, `1/8` and `13/100`.

The proof is the overlap-revealed route of Section 4. Its certificate is in
[overlap_revealed/](overlap_revealed/). Appendix A of the revised manuscript
adds a single-floor oriented comparison: it reproves dimension 22, gives
dimensions 16–20 and gives a closed-form result in dimension 25. Its
certificate is in [single_floor/](single_floor/).

Two later arguments extend the overlap-revealed route to lower dimensions,
with the same conclusion (both overlap signs percolate):

- the **macrostep route** replaces the uniform oriented paths of the global
  step by macrostep paths. Its certificates cover **dimension 9** at
  `tanh(beta) = 7/50`, **dimension 8** at `3/20` (with the local inputs of the
  overlap-revealed route, and again with the sharpened noisy-cavity inputs)
  and **dimension 7** at `31/200` and `3/20` (with the noisy-cavity inputs).
  They are in [macrostep/](macrostep/);
- the **noisy-cavity local inputs** sharpen the two local inputs of the
  overlap-revealed route (the pair certificate for the bond floor and the
  symmetrised Holley line) and keep its oriented global step. With them the
  route reaches **dimension 9** at `tanh(beta) = 3/20`, `29/200`, `31/200`,
  `4/25` and `7/50`. Their certificate is in
  [overlap_revealed/noisy_cavity/](overlap_revealed/noisy_cavity/). The same
  local deciders supply the dimension-8 and dimension-7 inputs of the
  macrostep route.

The directory READMEs use the theorem labels of the research write-ups. Both
routes call one of their results "Q9" (the macrostep result at `7/50` and the
noisy-cavity result at `3/20`); the manuscripts number them differently.

The first version's certificates remain at the top level. The dimension-22
fresh-star certificate at `tanh(beta) = 11/125` belongs to the first
version's Section 4 and Appendix A; the revised Appendix A describes that
argument but no longer uses it. The dimension-180 example of simultaneous
blue percolation and finite spin-glass susceptibility belongs to Appendix C
in both versions. The window constants and finite checks of the diluted-model
manuscript are in [diluted_uniqueness/](diluted_uniqueness/). Its explicit
windows are certified in
[diluted_uniqueness/window_d2/](diluted_uniqueness/window_d2/) (the
two-dimensional exact local-ratio window `g_2 = 2g* ≈ 2.3385e-5`) and
[diluted_uniqueness/window_general/](diluted_uniqueness/window_general/) (the
sparse-insulation windows in `d = 2, 3` and the octahedral closed form for
`d >= 3`).

This repository is private during manuscript preparation. The intended public
location is <https://github.com/PeaBrane/cmr-blue-percolation>, to be made public
with the paper. Each manuscript cites a specific commit so that its
computational inputs remain identifiable after later development.

## Reproduce

The default and supplementary programs need Python 3.10 or later and **only
the standard library**. They were validated with CPython 3.13.5. They need no
Peapods installation, downloaded input data, numerical libraries or network
access.

```sh
uv venv --python 3.13
.venv/bin/python run_all.py
```

To add the independent finite-model checks and crosschecks, which also use
only the standard library:

```sh
.venv/bin/python run_all.py --supplementary
```

To add the optional checks that need extra packages (mpmath, SymPy, NumPy,
networkx), install them and add `--optional-deps`:

```sh
uv pip install --python .venv/bin/python mpmath sympy numpy networkx
.venv/bin/python run_all.py --supplementary --optional-deps
```

A check whose package cannot be imported is skipped with a `SKIPPED` line.
None of these checks is an input to a proof.

The runner uses its own interpreter for every child process. All entrypoints
reject optimized Python execution, because assertions implement the decisive
comparisons. Successful runs end with one of:

- `PASS: 10 verification scripts completed.` (default run);
- `PASS: 18 verification scripts completed.` (with `--supplementary`);
- up to ten more scripts with `--optional-deps` (`PASS: 28 verification
  scripts completed.` when every optional package is installed).

Typical times with CPython 3.13.5 on an Apple Silicon laptop: about 95 s for
the default run, of which the four programs of the later arguments take about
75 s (`window_d2/check_window_d2.py` about 40 s, `certify_noisy_cavity.py`
14 s, `window_general/check_window_general.py` 12 s and
`certify_macrostep.py` 11 s); about 100 s more for `--supplementary`; and
about 5.5 min more for `--optional-deps`, mostly `crosscheck_global.py`
(about 75 s), `pairing_exact.py` (about 105 s) and the two brute-force
checks of the noisy-cavity certificates (about 40 s and 60 s). No default or
supplementary program needs more than about 210 MB of memory; the largest are
`window_d2/check_window_d2.py` (190–210 MB on macOS, about 150 MB on Linux)
and `macrostep/certify_macrostep.py` (about 105 MB). The optional
`pairing_exact.py` needs about 210 MB.

### Full reproductions

Three programs recompute stored certificate data from scratch. They are not
part of `run_all.py`, because they need a C compiler or compiled numerical
packages, several cores and minutes to tens of minutes:

```sh
uv pip install --python .venv/bin/python numpy numba networkx
.venv/bin/python macrostep/reproduce_engine.py --nproc 6 --l2 --out repro-macrostep
.venv/bin/python diluted_uniqueness/window_d2/reproduce_window_d2.py --threads 6 --workdir repro-window-d2
.venv/bin/python diluted_uniqueness/window_general/reproduce_window_general.py --jobs 6 --e2e
```

- `macrostep/reproduce_engine.py` (NumPy and Numba; recorded with NumPy 2.5.3
  and Numba 0.67.0) reruns the macrostep engine for the nine certificates
  used by the theorems. It checks that the recomputed engine outputs equal the
  stored ones bit for bit and re-decides the fresh certificates exactly. With
  `--l2` it also reruns the independent first implementation and compares the
  two entry by entry. It takes about 21 min of wall time with 6 worker
  processes on a 24-core x86-64 Linux machine; see
  [macrostep/README.md](macrostep/README.md).
- `diluted_uniqueness/window_d2/reproduce_window_d2.py` needs a C compiler
  with POSIX threads. It regenerates every record of the `d = 2` certificate
  with the C engine, including all 3,879,876 cell boundary conditions, and
  checks the canonical digests (about 10 min of wall time and 49 CPU minutes
  with 6 threads, 1.3 GB of memory, about 360 MB of engine dumps in the work
  directory).
- `diluted_uniqueness/window_general/reproduce_window_general.py` (networkx
  for one cross-check) reruns the heavier octahedral checks in `d = 3, 4, 5`,
  the clique cross-check and, with `--e2e`, the recorded end-to-end sanity
  runs, and compares every output byte for byte with
  `window_general/expected_outputs/` (about 18 min of wall time with `--e2e`,
  6.5 min without; up to 3.8 GB per `d = 4` process).

## What is checked

Default run:

| Program | Mathematical scope |
|---|---|
| [overlap_revealed/certify.py](overlap_revealed/certify.py) | Explicit dimensions 10, 11 and 12, nine rows. For each row it decides the conditions (C1)–(C6) of the revised Section 4: the bond floor `p <= p_B` (vertex sums), the Holley line in all `2d+1` environments, the hypotheses on the Ising comparison field, the mean-field test for the plus-density `rho_-`, the walk-sum ratio `lambda_* < 1`, and the second-moment criterion (`Score < 1`, or `Score_c < 1` at `d = 10`) with its constant θ\*. It asserts every entry of the table of certified constants, for example `Score <= 0.907068` and `theta_* >= 0.036292` at `d = 12`, `t = 3/25` and `Score_c <= 0.994576` at `d = 10`, `t = 27/200`, and the Section 4 text displays, among them the Holley slack, the Green-function bounds and the certified margins. It also checks that the parameter `kappa_2d` of the centered susceptibility bound exceeds 1 at every row, as the manuscript states. |
| [single_floor/verify_single_floor.py](single_floor/verify_single_floor.py) | Revised Appendix A: `p_0(2d, beta) = p_* > F_d` at `d = 16`–`20` and `22`, through the closed form at `k = 1` and the block reduction; the signed refinements at `d = 22` (floor `0.0586`) and `d = 16` (floor `0.0689`); the closed-form floor at `d = 25` and `26`; the appendix tables and displays. |
| [verify_signed_star.py](verify_signed_star.py) | First version's complete d22 fresh-star certificate: five signed-residual cases with 299 reduced cavity orientations, 205 other overlap/failure cases, five failure envelopes, all 63 proper triangular direction subsets, positive applicable increments, and the strict triangular critical inequality. |
| [verify_signed_star_independent.py](verify_signed_star_independent.py) | Separately implemented direct rational enumeration of the five signed d22 cases. Its accumulation uses absolute cavity fields. It does not check the other 205 cases or the global comparison. |
| [verify_physical_separation.py](verify_physical_separation.py) | d180 at `tanh(beta)=1/56`: the MNS site parameter exceeds `23/2000` and the cited site threshold; `kappa < 5/8` gives `chi_SG < 8/3`. |
| [diluted_uniqueness/gd_constants.py](diluted_uniqueness/gd_constants.py) | Window constants of the diluted model: `K_F`, `K_Q`, `N_1`; `g_2 = 2^-988/1130` in `(3.38, 3.39)·10^-301`; `g_2 > 10^-347` from `p_c >= 1/3` alone; `g_3 > 10^-7973` and `g_4 > 10^-152290`; and the classical bounds `beta_*(2)`, `beta_*(3)`. |
| [diluted_uniqueness/window_d2/check_window_d2.py](diluted_uniqueness/window_d2/check_window_d2.py) | Diluted model, `d = 2`: `g* = 1059411733898907047504026797/90604240507898863756433027512000`, `g_2 = 2g*` in `[2.3385477941433e-5, 2.3385477941434e-5]`, `beta_*(2)` in `[5.3316976645232, 5.3316976645233]`. It decides (C0), the geometry and coverage identity, the boxes near `S_L`, (C1-side), (C1-o) and (C1-int) on every committed Path-1 record, and (C2) at `u = 0, 1/2` on every committed reduced cell record; through the `D4` invariance these cover all 3,879,876 cell boundary conditions. It recomputes the intersection records and a subset of the side and origin records with an independent Python engine. That the committed cell records are the counts their definition prescribes is decided by the full reproduction, not by this program. |
| [diluted_uniqueness/window_general/check_window_general.py](diluted_uniqueness/window_general/check_window_general.py) | Diluted model, sparse insulation: the rule tables `T_3` and `T_2` (design hypotheses, exact costs, cliques, `K`), `Delta_3 >= 3.232324e-9` (`7.405526e-9` with vertex coins), `Delta_2 >= 4.552624e-8` (`7.153991e-8` with domino coins) and the `beta_*` thresholds; the octahedral rule in `d = 3` at `L = 5`; every displayed entry of the octahedral closed-form table (`d = 3..8, 10, 20, 100`) and of the axis-star coin table (`d = 3..8, 10, 12`), each as an exact floor. |
| [macrostep/certify_macrostep.py](macrostep/certify_macrostep.py) | Macrostep route in dimensions 9, 8 and 7: six instances and 25 stored certificates, nine of them used by the theorems. For each instance it decides the local inputs (tangent bond floor and Holley line at `d = 9, 8`; the noisy-cavity Lemmas P and S at `d = 8, 7`) and the hypotheses of the macrostep theorem, and recomputes `rho_-`, `rho_c`, `kappa`, `C_fin` and the torus sizes `L_0(n)`. For each certificate it checks the payload digests and the state list, and decides both conditions of the interaction-matrix criterion (near matrix, and near Green function against the far-region bound) in exact dyadic arithmetic, with the exact θ\*; for example `lambda <= 0.779209`, `theta_* >= 0.001908` at `d = 9`, `t = 7/50` and `lambda <= 0.924670`, `theta_* >= 0.0003267` at `d = 7`, `t = 31/200`. It asserts 276 displayed numbers. |
| [overlap_revealed/noisy_cavity/certify_noisy_cavity.py](overlap_revealed/noisy_cavity/certify_noisy_cavity.py) | Noisy-cavity proof of dimension 9 at five temperatures: the noisy-class constraint of Lemma N, the pair certificate of Lemma P (162 point conditions, 1377 chord polynomials, 18 class sums), the Holley line of Lemma S (1330 two-block classes in 19 environments), the strict comparison of the plain and noisy witnesses, and the oriented second-moment criterion with Score_c (at `t = 7/50`, Score′_c from a second implementation of the covariance lemma); for example `Score_c <= 0.9906982`, `theta_* >= 0.003173` at `t = 3/20`. It asserts 155 displayed numbers. |

Added by `--supplementary`:

| Program | Mathematical scope |
|---|---|
| [checks/smallstar_audit.py](checks/smallstar_audit.py) | Central-spin and activation sums, the nonnegative residual coefficient and orientation reduction on degrees 1–7 at three rational temperatures. |
| [checks/numerator_audit.py](checks/numerator_audit.py) | Independent tilted-numerator, pointwise maximum and weighted-Jensen checks, including rational cavity mixtures, on degrees 1–7. |
| [checks/check_susceptibility.py](checks/check_susceptibility.py) | 2,916 one-edge identities/inequalities and complete disorder/spin sums on four small graphs; 3,136 odd-set-pair bounds. |
| [overlap_revealed/crosscheck/crosscheck_local.py](overlap_revealed/crosscheck/crosscheck_local.py) | The first research code base, an independent implementation of the overlap-revealed local inputs. It reruns its own `(12, 3/25)` certificate, and reproduces all 639 local rationals of `certify.py` exactly (checked by digest). |
| [single_floor/check_lemmas.py](single_floor/check_lemmas.py) | The balanced-tilt lemma at 256 cases with the complete `h = 0` rows, the single-ratio identity, the monotone-pieces lemma, the signed enumeration against brute force, and exact agreement of the dimension-22 signed residuals with `verify_signed_star.py`. |
| [single_floor/check_side_results.py](single_floor/check_side_results.py) | Side results: exact witnesses that `k = 1` does not minimise the floor uniformly in `beta`, the same-cavity counterexample, and `d = 15` for overlap counts `k >= 2` (remark "Where the route stops"). |
| [diluted_uniqueness/check_typeS_rule.py](diluted_uniqueness/check_typeS_rule.py) | Exhaustive local check of the type-S surgery rule: all configurations in `d = 2` for `L = 10, 11, 13` (20,820 of them) and in `d = 3` for `L = 10` (1,780,110 of them). |
| [overlap_revealed/noisy_cavity/crosscheck/nesting_exact.py](overlap_revealed/noisy_cavity/crosscheck/nesting_exact.py) | Nesting of the noisy cavity classes, exactly in 108 small instances, and strictness of the plain/noisy witness comparison at the five dimension-9 temperatures. |

Added by `--optional-deps`, crosschecks and finite sanity checks rather than
proof inputs:

| Program | Package | Scope |
|---|---|---|
| [overlap_revealed/crosscheck/crosscheck_global.py](overlap_revealed/crosscheck/crosscheck_global.py) | mpmath | The first code base's global criterion in 200-bit interval arithmetic for all nine rows. It meets every displayed `rho_-`, `eta'`, Score, Score_c and θ\*. |
| [diluted_uniqueness/optional/](diluted_uniqueness/optional/) | mpmath, SymPy, NumPy, networkx | Interval recomputation of the window constants, symbolic decimation identities and exact Gibbs-reduction enumeration on small boxes, each with the counts displayed in the diluted-model manuscript, and the surgery map on (+)-pivotal edges of seeded `d = 2` samples; `surgery_counts.py --slow` adds `d = 3` and `--full`, which takes hours, reproduces the manuscript's totals. |
| [overlap_revealed/noisy_cavity/crosscheck/crosscheck_global_A.py](overlap_revealed/noisy_cavity/crosscheck/crosscheck_global_A.py) | mpmath | The first code base at the five noisy-cavity rows, including the covariance lemma for any `2dK' < 1`. It meets every displayed number that is the worse of the two code bases. |
| [overlap_revealed/noisy_cavity/crosscheck/](overlap_revealed/noisy_cavity/crosscheck/): `verify_pair_bruteforce.py`, `verify_holley_bruteforce.py`, `lemmaN_exact.py`, `pairing_exact.py` | NumPy | Independent implementations of the pair certificate and the Holley line with brute-force enumeration; exact enumeration checks of the noisy-cavity reduction and of the pairing identity and inequality chain of the pair certificate. |

[overlap_revealed/search/](overlap_revealed/search/) holds the floating-point
searches that chose the frozen inputs of `overlap_revealed/params.json`. They
need NumPy and SciPy and are not part of any verification run.
[macrostep/search/](macrostep/search/) repeats the floating-point choice of the
macrostep instance-1 tangent points (SciPy). That search is platform
dependent: the values it finds may differ from the frozen ones in later
digits, and any admissible values give a valid certificate.
[macrostep/checks/](macrostep/checks/) holds optional exact rechecks of engine
output (NumPy, Numba); they are not part of `run_all.py`.

## Arithmetic and expected results

Every decisive comparison is an assertion on exact integers or
`fractions.Fraction` values. Decimal values in the output are summaries,
rounded in the safe direction by the newer programs. No floating-point
comparison proves an assertion.

- The d22 signed cases use exact kernels. Other d22 cases use outward
  rational enclosures on a `2^-80` grid, with signs determining the direction
  of rounding.
- The overlap-revealed certificate uses rational enclosures of `exp`,
  `atanh`, square roots and `pi`, and rounds upper bounds outward on a
  `10^-50` grid. Its floating-point searches only chose frozen rational
  inputs; any admissible input gives a valid certificate.
- In the single-floor certificate every quantity is rational at rational
  `tanh(beta)`.
- The diluted constants are integers plus rational series enclosures of
  logarithms.
- The macrostep certificates read every stored binary64 number as the exact
  dyadic rational it represents and decide the criterion on integers. The
  near matrix, near Green rows and far-region bound are engine outputs,
  computed in rounded floating point with an a-priori error analysis. The
  default run takes them as stored. The full reproduction recomputes them bit
  for bit, and an independent implementation (`macrostep/engine_l2/`) agrees
  with them entrywise to about `5·10^-8` relative, on the upper side as its
  error margins require.
- `window_d2` commits its records as canonical xz streams. The SHA-256 of
  each decompressed stream is an order-independent digest of the record
  family, and the default check decides every inequality on them. The
  all-partition cell family (3,879,876 records) is regenerated by the full
  reproduction; its canonical digest is
  `7a69009a9acc7e31ff09a42f2e663e0b7f8bf72e2fc13037f5e3b660dfbdcb32`.
  `window_general/check_window_general.py` also compares the printed output
  of its component checkers with `window_general/expected_outputs/` as a
  regression check; its decisions are assertions on recomputed exact values.

The expected d22 score is approximately **1.000816029521908**, strictly larger
than one by exact rational comparison. [expected/dimension22.json](expected/dimension22.json)
records the exact score, selected parameters and output-vector digest:

```text
d7f1367cac657b8821c60093f30aac8932e2941175cb4dceac57b0450da49f7d
```

This digest hashes the ordered rational bound vector, not the program source.
In the same way:

- [expected/overlap_revealed.json](expected/overlap_revealed.json) records the
  directed bounds of every row and two digests: the local rationals, and the
  local rationals together with `rho_-`, `eta'`, Score, Score_c and θ\*;
- [expected/single_floor.json](expected/single_floor.json) records the
  single-floor bounds and a digest of the exact values;
- [expected/diluted_gd_constants.json](expected/diluted_gd_constants.json)
  records the window constants;
- [expected/diluted_window_d2.json](expected/diluted_window_d2.json) and
  [expected/diluted_window_general.json](expected/diluted_window_general.json)
  record the window constants, the canonical record digests and the table
  floors;
- [expected/macrostep.json](expected/macrostep.json) records the directed
  bounds of every instance and certificate and two digests: the local
  rationals, and the certificate values (`eta_F`, λ, `H`, `Gamma_w` and θ\*
  of both certificates of each run);
- [expected/noisy_cavity.json](expected/noisy_cavity.json) records the
  directed bounds of the five rows and a digest of the exact values.

Each file equals the JSON printed at the end of the corresponding program.
The Git commit identifies the source. The verifiers recompute every inequality
from raw parameters and do not read the expected files as oracles. The one
exception is `crosscheck_local.py`: it compares its own digest with the
recorded one to confirm that the two code bases agree. Elapsed times in the
outputs are not reproducibility targets.

## Relation to the proofs

- **Revised manuscript, Section 4 (overlap-revealed route).** The manuscript
  proves the gauge representation and the overlap-revealed star law, the
  extremal partition and the bond floor after failures, overlap odds and
  one-sided Holley domination, the mean-field density lemma, the Ising pair
  ratio, the weighted second moment with its renewal bound, uniformity in the
  torus size, transfer to selected limits and almost-sure coexistence.
  `overlap_revealed/certify.py` certifies the finite inequalities (C1)–(C6)
  that those statements reduce to at each certified row. Its frozen inputs
  are in `overlap_revealed/params.json`.
- **Revised Appendix A (fresh-star route and single-floor comparison).** The
  appendix proves the star estimates, the single-floor theorem, the collision
  Green function bounds, the single-ratio identity for the row `h = 0`, the
  block reduction to `k = 1` and the closed-form star floor.
  `single_floor/verify_single_floor.py` certifies the resulting rational
  inequalities.
- **First version, Section 4 and Appendix A.** The cumulative comparison,
  with the fresh-star estimates (exact normalized masses, signed and
  nonnegative residual bounds, and finite sums), gives the dimension-22
  theorem certified by `verify_signed_star.py`. The revised manuscript
  reproves that dimension by the single-floor comparison.
- **Appendix C (both versions).** The d180 integer inequalities of the
  susceptibility example are checked by `verify_physical_separation.py`.
- **Macrostep route (dimensions 9, 8, 7).** The manuscript proves the
  geometry of macrostep paths, the weighted second moment for macrostep
  paths, the interaction-matrix criterion with a far region, and the transfer
  to the torus bound of Section 4; the sharpened local inputs at `d = 8` and
  the inputs at `d = 7` come from the noisy-cavity lemmas.
  `macrostep/certify_macrostep.py` certifies the finite inequalities: local
  inputs, hypotheses, and conditions (i)–(ii) of the criterion on the stored
  engine outputs.
- **Noisy-cavity proof of dimension 9.** The manuscript proves the
  noisy-cavity reduction, the pair certificate and the symmetrised Holley
  certificate. `overlap_revealed/noisy_cavity/certify_noisy_cavity.py`
  certifies their finite inequalities and the unchanged oriented criterion at
  the five temperatures.
- **Diluted-model manuscript.** The manuscript proves the decimation lemma,
  the Gibbs-state reduction and the strict inequality. `gd_constants.py`
  computes the displayed window constants; the other programs directly under
  `diluted_uniqueness/` are finite sanity checks. The two-dimensional window
  theorem reduces to the finite exact inequalities (C0)–(C2) on explicit
  boxes, which `window_d2/` certifies. The sparse-insulation theorems reduce
  to finite statements about the explicit rule tables and to exact
  evaluations of the closed form, which `window_general/` certifies. The hand
  proofs (planarity reduction, curve, Farkas split, `D4` invariance, design
  lemma, counting lemmas, octahedral rule) are in the manuscript.

The programs certify finite computations only. The analytic conditioning,
probability and infinite-volume arguments are supplied by the manuscripts. The
supplementary and optional small-graph and small-box checks are independent
finite evidence for identities used there.

The repository contains the computational dependencies of the explicit
certificates in both manuscripts. The parameter searches in
`overlap_revealed/search/` and `macrostep/search/` are included because they
chose frozen inputs. Historical dimension searches, exploratory geometry,
private worklogs and reference PDFs are maintained separately.

## Cite

Use [CITATION.cff](CITATION.cff) and include the commit hash recorded by the
manuscript. This is a computational companion to the manuscripts, not a formal
proof-assistant development.

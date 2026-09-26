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

The first version's certificates remain at the top level. The dimension-22
fresh-star certificate at `tanh(beta) = 11/125` belongs to the first
version's Section 4 and Appendix A; the revised Appendix A describes that
argument but no longer uses it. The dimension-180 example of simultaneous
blue percolation and finite spin-glass susceptibility belongs to Appendix C
in both versions. The window constants and finite checks of the diluted-model
manuscript are in [diluted_uniqueness/](diluted_uniqueness/).

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

- `PASS: 6 verification scripts completed.` (default run);
- `PASS: 13 verification scripts completed.` (with `--supplementary`);
- up to five more scripts with `--optional-deps`.

Typical times with CPython 3.13.5 on an idle Apple Silicon laptop: about
16 s for the default run, about 100 s more for `--supplementary`, and about
95 s more for `--optional-deps` (about 75 s of it `crosscheck_global.py`).
On a loaded machine an independent reproduction measured about 40 s, 150 s
and 200 s respectively, with no single program above about 60 s. No program
needs more than 80 MB of memory.

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

Added by `--optional-deps`, a crosscheck and finite sanity checks rather than
proof inputs:

| Program | Package | Scope |
|---|---|---|
| [overlap_revealed/crosscheck/crosscheck_global.py](overlap_revealed/crosscheck/crosscheck_global.py) | mpmath | The first code base's global criterion in 200-bit interval arithmetic for all nine rows. It meets every displayed `rho_-`, `eta'`, Score, Score_c and θ\*. |
| [diluted_uniqueness/optional/](diluted_uniqueness/optional/) | mpmath, SymPy, NumPy, networkx | Interval recomputation of the window constants, symbolic decimation identities and exact Gibbs-reduction enumeration on small boxes, each with the counts displayed in the diluted-model manuscript, and the surgery map on (+)-pivotal edges of seeded `d = 2` samples; `surgery_counts.py --slow` adds `d = 3` and `--full`, which takes hours, reproduces the manuscript's totals. |

[overlap_revealed/search/](overlap_revealed/search/) holds the floating-point
searches that chose the frozen inputs of `overlap_revealed/params.json`. They
need NumPy and SciPy and are not part of any verification run.

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
  records the window constants.

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
- **Diluted-model manuscript.** The manuscript proves the decimation lemma,
  the Gibbs-state reduction and the strict inequality. Only the displayed
  window constants are computed here; the other checks are finite sanity
  checks.

The programs certify finite computations only. The analytic conditioning,
probability and infinite-volume arguments are supplied by the manuscripts. The
supplementary and optional small-graph and small-box checks are independent
finite evidence for identities used there.

The repository contains the computational dependencies of the explicit
certificates in both manuscripts. The parameter searches in
`overlap_revealed/search/` are included because they chose frozen inputs.
Historical dimension searches, exploratory geometry, private worklogs and
reference PDFs are maintained separately.

## Cite

Use [CITATION.cff](CITATION.cff) and include the commit hash recorded by the
manuscript. This is a computational companion to the manuscripts, not a formal
proof-assistant development.

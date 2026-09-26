# Overlap-revealed route: dimensions 10, 11 and 12

These programs certify the finite inequalities behind the explicit-dimension
theorem of the revised manuscript (v2, Section 4, "An explicit dimension
through overlap-revealed exploration"). The theorem covers iid fair ±1
couplings and every selected periodic joint limit. Both overlap signs have
root blue-percolation probability at least θ\* > 0 and infinite blue
components almost surely, in these cases:

- `d = 12` at `t = tanh β = 3/25`, `11/100`, `23/200`, `1/8` and `13/100`;
- `d = 11` at `t = 13/100` and `3/25`;
- `d = 10` at `t = 27/200` and `13/100`.

For each row of the manuscript's table of certified constants, the inputs are
exact rationals: a bond floor `p`, a Holley line `g_K = e^{2K'}`,
`g_h = e^{2h}`, and a mean-field point `mbar`.

## Certificate

[certify.py](certify.py) uses only the standard library. It reads the frozen
inputs in [params.json](params.json) and decides, for every row, the
conditions (C1)–(C6) listed in Section 4.6 of the manuscript ("The certified
constants"):

| Condition | Manuscript result | Exact check | Size per row |
|---|---|---|---|
| (C1) `p <= p_B` | Proposition "A finite certificate for the floor", after the theorem on the extremal partition | vertex sums `Delta_j <= 0`, `j = 0..2d-1`, with frozen tangent points `c_m`; also `p` lies below the aligned frozen-cavity witness | `2d` double sums of at most `d(d+1)` terms |
| (C2) Holley line in all `2d+1` environments | Lemmas "Numerator" and "Denominator", Proposition "Certified Holley line" | `C^S Num_k / Den_k >= g_h g_K^S` for `S = 2k-2d`, `k = 0..2d`; `Num_k` with frozen tangent weights `w_k`; `Den_k` as the maximum over every frozen class, asserted to sit at the class `floor((2d-k)/2)` named in the lemma | `sum_k (k+1)` numerator double sums (325 at `d = 12`) and `sum_k (2d-k+1)` denominator convolutions |
| (C3) hypotheses | hypotheses of the second-moment subsection | `g_h < 1 <= g_K`, `g_K g_h <= 1`, `g_K^{2d} < 2718/1000 < e` (so `4dK' < 1`); also `g_K^d < 2718/1000` and `d >= 6` | rational powers |
| (C4) Ising plus-density | Lemma "Mean-field lower bound on the plus-density" | `tanh(abs(h) + 2dK' mbar) < mbar`, giving `rho_- = (1 - mbar)/2` | three `atanh` enclosures |
| (C5) walk sums | Lemma "Walk-sum bounds" | `lambda_* = d tbar/(1 - 2d tbar) < 1` with `tbar = tanh K' = (g_K-1)/(g_K+1)` exactly, and `B(m+1) <= lambda_* B(m)` for `m < 86` | closed forms |
| (C6) second moment | Theorem "Second-moment criterion", Lemmas "Collision Green function", "Cycle weights" and "Correlation gain along shared edges" | `eta' < 1`; `Score < 1` for `d = 11, 12`; `rho_-^2 >= c_cov/4` and `Score_c < 1` for `d = 10` | `u_k` for `k < 4d`, `B(m)` for `m <= 86`, `psibar(r)` and `abs(A_r)` for `r <= 40`, closed-form tails |
| θ\* | root-probability constant, Corollary "Uniformity in the torus size" | `theta_* = e^{-1/100} rho_- (1-eta')(1-Score)/(1-1/d)`, with `e^{-1/100} >= 99/100`, at least the displayed value | one product |

Then `certify.py` asserts every number that the manuscript displays for these
rows. Each row's `published` block in `params.json` lists them, keyed
`<quantity>_<direction>`, where the direction is `upper`, `lower`, `exact`,
or `loglower` (a lower bound on a logarithm, checked as `value >= e^shown`
with an upper bound on `e^shown`):

- the table: the inputs `p`, `g_K`, `g_h` as its terminating decimals, and
  `rho_-`, `eta'`, Score, Score_c and θ\* for every row;
- the certified margins: `1 - Score` is about 9% at `(12, 3/25)` and ranges
  from about 7.3% to 9.5% in dimension 12, it is about 4.8% and 3.9% in
  dimension 11, and `1 - Score_c` is 0.54% and 0.27% in dimension 10, where
  the shared-edge gain `rho_c - rho_-` is about 0.006;
- at `(12, 3/25)`: `p_A = 75/196`, `tanh 2β = 75/317`, the closed floor
  `p_A/2 = 75/392` and the Score of about 1.02 that it would give; the largest
  vertex sum `Delta_20 ≈ -1.95e-4`; the witness `p_B <= 0.2210000` and the
  frozen-success values `0.2186899`, `0.2163712`, `0.2140597`, `0.2117739`
  (remark "Observed successes are not harmless"); `K' ≈ 0.013000` and
  `h ≈ -0.090000`; the Holley slack, at least `0.001827` at `S = 24`,
  `0.00192` at `S = -24` and `0.0057` in the other 23 environments; the worst
  conditional plus probability of about 0.31; `G_12` in
  `[1.1011537095, 1.1011552715]`, `F_12 - 1/12 <= 0.0085295258` and
  `G^(2)_12 <= 1.2272959943`; the three Score terms `<= 0.886231`,
  `0.019603`, `0.001235`; `kappa <= 1.29818` and `alpha_D <= 0.311972`;
- `F_11 - 1/11 <= 0.0104053880` and `F_10 - 1/10 <= 0.0130075189`;
- the boost term of Score, at most `0.0019` in dimension 12 and `0.0034` in
  dimension 10;
- in dimension 10, basic Score bounds of at least `1.0077` and `1.0098`: the
  basic criterion certifies neither row;
- at every row, the bond floor is "about `0.9 tanh 2β`": `p / tanh 2β >= 0.85`
  and the witness satisfies `witness / tanh 2β <= 0.95`;
- the parameter `kappa_2d = 2d B A^(2d-1)` of the manuscript's centered
  susceptibility bound, with `A, B = (cosh 4β ± cosh 2β)/2`: it is about 5.5
  at `(12, 3/25)` (remark "Scope and consequences"), and `certify.py` asserts
  `kappa_2d > 1` at every row, as the proof of the dimension-180 separation
  corollary states, so that bound does not apply at these points.

Finally, it evaluates the global criterion at the research row `d = 10`,
`t = 3/25`, which the manuscript reports as not certified, and asserts that
the Score_c bound lies in `[1.0192, 1.0193]`, above one.

The floating-point observations of the manuscript's remark "Dimension nine;
numerical observations" come from scans in the research notes
(`assembly/explore_t.py`, `revision/d9_headroom.py`); they are not
certificates and are not reproduced here.

The run ends by printing a JSON summary equal to
[../expected/overlap_revealed.json](../expected/overlap_revealed.json). It
records the directed decimal bounds of every row and two SHA-256 digests: of
the 639 local rationals (`Delta_j`, witness, `Num_k`, `Den_k`), and of those
together with `rho_-`, `eta'`, Score, Score_c and θ\*.

| d | t | p | g_K | g_h | rho_- ≥ | eta' ≤ | Score ≤ | Score_c ≤ | theta\* ≥ |
|---|---|---|---|---|---|---|---|---|---|
| **12** | **3/25** | 2161/10000 | 102634/100000 | 83527/100000 | 0.4351284 | 0.168984 | **0.907068** | 0.897470 | **0.036292** |
| 12 | 11/100 | 2023/10000 | 255577/250000 | 873143/10^6 | 0.4540507 | 0.127091 | 0.926888 | 0.918156 | 0.031296 |
| 12 | 23/200 | 2094/10000 | 1024283/10^6 | 855173/10^6 | 0.4453800 | 0.146071 | 0.913716 | 0.904575 | 0.035441 |
| 12 | 1/8 | 2224/10000 | 1028469/10^6 | 407843/500000 | 0.4240888 | 0.196297 | 0.905149 | 0.895074 | 0.034915 |
| 12 | 13/100 | 2283/10000 | 515339/500000 | 794291/10^6 | 0.4111645 | 0.230516 | 0.910357 | 0.899716 | 0.030630 |
| **11** | **13/100** | 2310/10000 | 6443/6250 | 25553/31250 | 0.4252374 | 0.215864 | **0.952104** | 0.940625 | **0.017392** |
| 11 | 3/25 | 2180/10000 | 513251/500000 | 427717/500000 | 0.4454995 | 0.162436 | 0.960915 | 0.950461 | 0.015882 |
| **10** | **27/200** | 2399/10000 | 516713/500000 | 102809/125000 | 0.4279573 | 0.236224 | 1.007739 | **0.994576** | **0.001950** |
| 10 | 13/100 | 2336/10000 | 1031101/10^6 | 840843/10^6 | 0.4379919 | 0.205637 | 1.009832 | 0.997261 | 0.001048 |

These are the manuscript's displays. Each is the worse of two independent
code bases, rounded in the safe direction, so the exact bounds of
`certify.py` (code base B) are at least as good; they sometimes differ in the
last digit, for example `rho_- >= 0.4540508` and `theta_* >= 0.036293`.

## Arithmetic

All algebraic quantities are exact `fractions.Fraction` values. The
transcendental quantities have rational enclosures:

- `exp`: a 31-term Taylor sum plus the remainder `3 x^31/31!`, for `0 <= x <= 1`;
- `atanh`: 26 series terms plus a geometric tail. This gives
  `K' = atanh((g_K-1)/(g_K+1))` and `abs(h) = atanh((1-g_h)/(1+g_h))`;
- square roots: integer square roots with 30 decimal digits;
- `pi >= 3.14159265358979`, checked with Machin's formula; `2718/1000 < e`
  and `e^{-1/100} >= 99/100` are checked from Taylor sums.

The mean-field test is decided in the logarithmic form
`abs(h) + 2dK' mbar < atanh(mbar)` with these enclosures. The manuscript's
equivalent power form `g_h^{-N_2} g_K^{N_1} <= ((1+mbar)/(1-mbar))^{N_2}`,
where `2d mbar = N_1/N_2`, has exponents in the millions at the frozen `mbar`
(`N_1 = 3892293`, `N_2 = 1250000` at `(12, 3/25)`); `crosscheck_local.py`
decides the power form at a four-decimal point `x >= mbar`.

The global criterion uses `tanh K' = (g_K-1)/(g_K+1)` exactly. `c_cov` uses
the identities `e^{2M} = g_K^{2d-1}/g_h` and `2 sinh 2K' = g_K - 1/g_K`, and
only `e^{K'} = sqrt(g_K)` is enclosed; a lower bound for `c_cov` may replace
it in the shared-edge lemma. Upper bounds are rounded outward to a `10^-50`
grid. Floating point never decides a comparison.

## Frozen parameters

`params.json` stores exact decimal or fraction strings:

- `p`, `g_K`, `g_h` and `mbar`;
- the bond-floor tangent points, either as `c_m` (`"direct"`) or as
  multipliers `rho_m` with `c_m = rho_m * omega(m)` (`"rho_times_omega"`,
  row `(12, 3/25)` only; these are the multipliers displayed in the
  manuscript), where `omega(m) = w^((1+m)/2) + w^(-(1+m)/2)` and
  `w = (1+t)/(1-t)`;
- the Holley tangent weights `w_0, ..., w_2d`;
- the criterion used, `basic` or `refined_c`;
- the displayed numbers (`published`);
- the uncertified row `(10, 3/25)` (`uncertified_rows`).

Floating-point searches chose the tangent points, the weights and `mbar`.
Any positive tangent point and any weight in `[0, 1]` are admissible, and any
`mbar` that passes the exact root test is valid, so the search cannot affect
soundness. At `(12, 3/25)` the tangent points and weights are those of the
first research code base, for which the manuscript displays the floor
multipliers and the Holley slack; the other rows use the searches of code
base B. `provenance` records how each value was found and is not read by the
verifier.

## Independent crosscheck (`crosscheck/`)

This directory holds code base A of the research notes: the local checkers
(`certify_local.py`, `certify_pB.py`, `certify_odds.py`) and the global
criterion (`criterion.py`, `green.py`, `covariance.py`, `evaluate_part1.py`).
They were written independently of `certify.py`, which follows code base B.
The files are copies with two changes: each rejects `python -O`, and
`certify_pB.certify` accepts explicit tangent points.

- [crosscheck_local.py](crosscheck/crosscheck_local.py) uses only the
  standard library and runs with `run_all.py --supplementary`. It reruns the
  research-log certificate for `(12, 3/25)`
  ([certify_local.py](crosscheck/certify_local.py), whose output matches the
  research log `results_headline.txt`) and confirms that its tangent points
  and weights are the frozen ones. It then checks every row with the frozen
  inputs, including the root test in the exact power form
  `g_h^{-N_2} g_K^{N_1} < ((1+x)/(1-x))^{N_2}` of the manuscript, at
  `x = ceil(mbar * 10^4)/10^4`. Finally it confirms that code base A's 639
  local rationals hash to `local_vector_sha256` in the expected JSON, which is
  exact agreement with code base B.
- [crosscheck_global.py](crosscheck/crosscheck_global.py) needs mpmath and
  runs with `run_all.py --optional-deps`. It encloses `K'` and `abs(h)` in
  200-bit outward-rounded intervals, finds its own `mbar`, takes `tbar = K'`
  and sums to radius 60. It checks that its bounds meet every displayed
  `rho_-`, `eta'`, Score, Score_c and θ\*. It also evaluates the sharper
  Score'_c of the research notes (global.md Lemma 3.6, which needs `L >= 40`
  and `2dK' <= 0.35`); no theorem of the manuscript uses it. Row labels as
  arguments restrict it to those rows, for example
  `python crosscheck/crosscheck_global.py d12-t3/25`.

## Parameter search (`search/`, optional)

These programs need NumPy and SciPy and are not needed for verification.

- [freeze_params.py](search/freeze_params.py) re-derives the searched
  values of `params.json` (tangent points, weights and `mbar`) with code base
  B's searches, rebuilds the rest from its tables, including the `published`
  blocks, and compares the result with the frozen file; `--write PATH` writes
  the fresh file. It takes about 35 CPU seconds and 75 MB of memory.
- [choose_rows.py](search/choose_rows.py) holds the searches that chose
  `(p, g_K, g_h)`: the first code base's procedure (research notes
  `scripts/certify_all.py`, whose output is the research log
  `results_nearby.txt`) at its nine research points, eight certified rows and
  the uncertified `(10, 3/25)`, and the code base B procedure of
  `assembly/certify_new_point.py` for `(10, 27/200)`. It ends by comparing
  its choices with `params.json`. At `(12, 3/25)` only `p` comes from this
  search; the line there is the exploration-phase line `K' = 0.013`,
  `h = -0.09` with `e^{2K'}` and `e^{2h}` rounded down to five decimals. Its
  floating line search tries 30001 values of `K'` per row, so a full run takes
  tens of CPU minutes.
- `odds_float.py`, `pB_family.py` and `explore_pB_opt2.py` are the first code
  base's floating helpers, without their exploration drivers.

## Relation to the proof

The programs certify the finite inequalities only. Section 4 of the
manuscript supplies the analytic steps:

- the gauge representation and the overlap-revealed star law;
- the extremal partition and the bond floor after failures;
- overlap odds, the single-site bound and one-sided Holley domination;
- the mean-field density lemma (GHS, GKS) and the Ising pair ratio (FKG,
  Dobrushin comparison);
- the weighted Paley–Zygmund bound, the renewal argument with the Schur test,
  and uniformity in the torus size;
- the overlap-revealed exploration, transfer to selected limits and
  almost-sure coexistence.

Runtimes with CPython 3.13.5 on an idle Apple Silicon laptop, in CPU seconds
(roughly doubled on a loaded machine): `certify.py`
about 8 s, `crosscheck_local.py` about 5 s, `crosscheck_global.py` about
75 s (6–10 s per row), each under 40 MB of memory.

## Research-note sources

The research notes behind the manuscript use other file names. In this
repository:

- `certify.py` replaces code base B, `assembly/indep_local.py`,
  `assembly/indep_global.py` and `assembly/certify_new_point.py`, with the
  searches moved to `search/`;
- `crosscheck/` holds code base A: `scripts/certify_local.py`,
  `scripts/certify_pB.py`, `scripts/certify_odds.py`, `criterion.py`,
  `green.py`, `covariance.py` and `evaluate_part1.py`. `crosscheck_local.py`
  and `crosscheck_global.py` run it at every row, including
  `(10, 27/200)`, which `assembly/crosscheck_new_point.py` checked. They also
  take the place of `assembly/directed_table.py` and
  `assembly/theta_table.py`: the table is the worse of the two code bases, so
  each code base must meet every displayed bound, which `certify.py` and
  `crosscheck_global.py` assert;
- `search/` holds the searches of `scripts/certify_all.py`,
  `scripts/explore_pB_opt2.py`, `scripts/pB_family.py`,
  `scripts/odds_float.py` and `assembly/certify_new_point.py`.

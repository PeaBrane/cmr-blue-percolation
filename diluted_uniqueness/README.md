# Diluted ±1 EA: uniqueness at every temperature slightly above p_c

These checks accompany the companion manuscript on the diluted ±1
Edwards–Anderson model, Yan Ru Pei, *Exact frustration cancellation and
uniqueness at all positive temperatures for diluted ±J spin glasses slightly
above the percolation threshold* (2026). Its couplings are iid, equal to `0`
with probability `1-p` and to `±1` with probability `p/2` each.

**Main theorem.** For `d >= 2` and `p < p_c + g_d/2`:

- the model has a unique Gibbs state at every inverse temperature;
- `E<sigma_0 sigma_x>^2` decays exponentially, uniformly in `beta`;
- the thermodynamic Edwards–Anderson response vanishes.

For `p > p_c` the `|J|`-graph percolates. The mechanism is the exact
cancellation of frustrated isolated diamonds.

## Decisive constants

[gd_constants.py](gd_constants.py) uses only the standard library and runs in
the default `run_all.py` (under 1 s). With surgery radius 4 the window uses

    K_F = 12 d 11^(d-1),  K_Q = C(d,2) 14^2 13^(d-2),  N_1 = d (4 * 11^d + 9^d),
    g_d = (p_c/2)^K_F / (N_1 2^(K_F + K_Q)).

The program checks, with exact integers and rationals, the values in the
manuscript's appendix on explicit constants:

- for `d = 2, 3, 4`, `K_F` equals the number of edges meeting the box `B_5`
  and `K_Q` the number of plaquettes with a corner in `B_6`, both counted by
  enumeration as in the manuscript's lemma on the sizes of the surgery sets
  (`N_1` is a union bound there); and
  `(K_F, K_Q, N_1, 3K_F + K_Q) = (264, 196, 1130, 988)` for `d = 2` and
  `(4356, 7644, 18159, 20712)` for `d = 3`;
- `p_c >= 1/(2d-1)` gives `g_2 > 10^-347` and `g_3 > 10^-7973`, each at the
  stated power of ten (for `d = 3` the integer `10^K_F 2^(K_F+K_Q) N_1` has
  7973 digits); the program also checks `g_4 > 10^-152290`, a value of the
  research notes that the manuscript does not display;
- with `p_c(Z^2) = 1/2`, `g_2 = 4^-264/(1130 * 2^460) = 2^-988/1130`, the
  window width is `g_2/2 = 2^-989/1130`, `log10(1/g_2) = 300.4707...`, and
  `3.38e-301 < g_2 < 3.39e-301`;
- `beta_*(d) = (1/2)[(3K_F + K_Q) ln 2 + ln N_1]` lies in
  `[345.929693652, 345.929693653]` for `d = 2` and in
  `[7183.13566267, 7183.13566268]` for `d = 3`, so
  `1/beta_*(2) < 2.8908e-3` and `1/beta_*(3) < 1.3922e-4`. Throughout the
  window, `beta_*(d)` is a lower bound for the classical uniqueness threshold
  `beta_l(p) = (1/2) ln(p/(p - p_c))`.

The program encloses logarithms with the series `ln x = 2 atanh((x-1)/(x+1))`
and a geometric tail bound. The percolation inputs `1/(2d-1) <= p_c <= 1/2`
and `p_c(Z^2) = 1/2` come from the literature cited in the manuscript. The run
ends by printing a JSON summary equal to
[../expected/diluted_gd_constants.json](../expected/diluted_gd_constants.json).

## Supplementary sanity check

[check_typeS_rule.py](check_typeS_rule.py) uses only the standard library and
runs with `run_all.py --supplementary` (about 40 s). It checks the local
facts (S1)–(S6) of the type-S surgery rule, in original coordinates, for every
configuration `(z, x', x)`, and asserts the counts displayed in the
manuscript:

| d | L | configurations |
|---|---|---|
| 2 | 10 | 5,820 |
| 2 | 11 | 6,660 |
| 2 | 13 | 8,340 |
| 3 | 10 | 1,780,110 |

It is a sanity check, not a proof input.

## Optional finite sanity checks

These scripts need extra packages and run with `run_all.py --optional-deps`
when the package can be imported. They are finite checks of statements that
the manuscript proves by hand, and they are not proof inputs. Each asserts the
counts that the manuscript's appendix displays; for the surgery counts this
needs `surgery_counts.py --full` (see below). Together they take about 20 CPU
seconds and under 80 MB of memory.

| Script | Package | What it checks |
|---|---|---|
| [optional/check_classical_window.py](optional/check_classical_window.py) | mpmath | `g_2`, `g_3`, `beta_*(2)` and `beta_*(3)` in 50-digit outward interval arithmetic, inside the displayed intervals: an independent recomputation of `gd_constants.py` |
| [optional/check_decimation_exact.py](optional/check_decimation_exact.py) | SymPy | the isolated-diamond decimation identities as Laurent-polynomial identities in `e^beta`: 16 sign patterns (the 8 frustrated ones cancel), 36 single-path cases, 256 mixed-magnitude patterns (48 cancel) |
| [optional/check_gibbs_reduction.py](optional/check_gibbs_reduction.py) | NumPy | exact integer enumeration on free-boundary boxes of `Z^2` with at most 16 spins and `e^beta` in `{2, 3}` (160 box–temperature pairs, fixed seed): all 2,794 predicted zero correlations vanish, and 812 of them are nonzero-coupling connected, so they vanish only through cancellation |
| [optional/surgery_counts.py](optional/surgery_counts.py) | networkx | reruns [optional/check_surgery.py](optional/check_surgery.py), a literal implementation of the surgery map, which enumerates every (+)-pivotal edge of each sampled state and checks that the output plaquette is s-pivotal, together with every intermediate property |

The manuscript reports 6,546 (+)-pivotal instances in `d = 2` (739 resolved in
Step 1; 2,285, 1,892 and 1,630 of types N, B and S) and 2,197 in `d = 3`
(442, 793, 384 and 578). They are the sums over five runs of
`check_surgery.py` (`python optional/check_surgery.py d samples seed
[random|path] [L]`). The research logs did not record the seeds. Each seed
was recovered by rerunning the program and matching the first line of its
log, the cumulative counts after 10 samples; `surgery_counts.py` lists the
five runs with their seeds (13, 12, 11, 15 and 14).

- By default it runs the first 10 samples of the three `d = 2` runs, 10 to
  20 CPU seconds, and compares the counts with the first lines of their logs.
- `--slow` adds the first 10 samples of the two `d = 3` runs. This takes
  minutes: the first 2 samples at `L = 10` alone take about 15 CPU seconds.
- `--full` completes all five runs, compares their final counts with the last
  lines of the logs and checks the manuscript's totals; this takes hours.

## Relation to the proof

The manuscript proves the decimation lemma, the Gibbs-state reduction and the
strict inequality. The strict inequality follows the Aizenman–Grimmett and
Grimmett–Stacey template, and the surgery is new. Of the programs here, only
`gd_constants.py` supplies displayed constants; the rest are sanity checks.
Runtimes were measured with CPython 3.13.5 on Apple Silicon.

The programs keep the research-note file names. `gd_constants.py` now also
recounts `K_F` and `K_Q` and checks, in rational arithmetic, the displays that
the research notes computed with mpmath in `check_classical_window.py`;
`surgery_counts.py` is new.

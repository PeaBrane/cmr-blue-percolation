# Sparse-insulation windows: d = 2 box, d = 3 cube, octahedral rule for d >= 3

These programs certify the explicit windows of the sparse pendant-arc
surgery in the diluted-model manuscript (Section 5; Theorem 1.2(ii) and (iii),
Theorems 5.26, 5.31 and 5.32, Table 1 and Remark 5.27). With `Delta_d = g_d/2` the width of
the window (`p < p_c + Delta_d`), the certified values are:

| d | construction | `Delta_d >=` | `g_d = 2 Delta_d >=` | `beta_*(d) >=` |
|---|---|---|---|---|
| 3 | unit-cube region, randomized rule `T_3`, plaquette coins (Theorem 5.31) | `3.232324e-9` | `6.464648e-9` | `8.9971` |
| 3 | same, vertex coins | `7.405526e-9` | `1.481105e-8` | `8.5826` (`T < 0.1166`) |
| 2 | `3 x 3` box, randomized rule `T_2`, plaquette coins (Theorem 5.32) | `4.552624e-8` | `9.105248e-8` | |
| 2 | same, domino coins | `7.153991e-8` | `1.430798e-7` | `7.8799` (`T < 0.1270`) |
| >= 3 | octahedral hand rule, closed form (Theorem 5.26) | Table 1 | e.g. `g_3 >= 4.6135e-10`, `g_4 >= 7.4257e-12`, `g_5 >= 3.1056e-13`, `g_10 >= 2.8748e-17` | Table 1 |
| 4..12 | same, axis-star coins (Remark 5.27) | coin table | e.g. `g_4 >= 1.2473e-11`, `g_5 >= 5.9232e-13`, `g_10 >= 7.1540e-17` | |

Here `beta_*(d) = (1/2) ln(2 p_c^- / g_d)` with a rigorous lower bound
`p_c^- <= p_c`: `1/2` for `d = 2`, `10000/47387` (`1/mu_3`, `mu_3 <= 4.7387`)
for `d = 3`, and the bracket end `p_-` of Table 1 for the octahedral rule. Every
displayed window is a floor of the exact rational value.

The proofs of all lemmas are by hand (manuscript). The computer-assisted
parts are finite checks of explicit rule tables with exact rational
arithmetic: for each covering configuration, designs (sets of insulated
vertices and arcs) satisfying the hypotheses (H1)–(H6) of Definition 5.5
(the design lemma is Lemma 5.6),
their exact key costs `c_kappa(p) = rho_hi^n p^-|O| (1-p)^-|W*\O|`, the
compatibility graph of keys and its maximal cliques, which give the
multiplicity constant `K = max(K_B, K_N, K_BN, K_S)`, the counts `q`, `N_Q` (or
`q_c`, `N_c` for coins), and

    Delta = (1 - 2^-(q+1)) / ((q+1) (K + N_Q / (2 (1 - p_hi)))).

For general `d` the octahedral hand rule gives the closed form (5.4)

    K_d = 2d [ max( c(7,10d-12) / (2(d-2)), c(6,8d-9) ) + c(6,8d-9) ],
    c(a,b) = max over p in {p_-, p_+} of p^-a (1-p)^-b,
    q = 4 d^2 (d-1) (2d^2-3d+4) / 3,   N_Q = 8 d^2 (2d^2-3d+4) / 3,

with `p_+` from the Gomes–Pereira–Sanchis bound (exact bisection plus
`10^-6`) and `p_- = 1/mu_d - 10^-6` (`d = 3, 4`) or `1/(2d)` (`d >= 5`).

## Files

| path | role |
|---|---|
| [check_window_general.py](check_window_general.py) | **default light check** (standard library) |
| [reproduce_window_general.py](reproduce_window_general.py) | **full reproduction** (standard library + networkx) |
| [cr_check.py](cr_check.py) | independent checker of a rule table: design hypotheses, exact costs, weights, compatibility with all anchor variants, own Bron–Kerbosch cliques, `K_B`, `K_N`, `K_BN`, the type-S rule and its bound, `q`, `N_Q`, coin counts, `Delta` |
| [oct_check.py](oct_check.py) | independent check of the octahedral hand rule in a given `d`: every configuration (types B, N, and S at a given `L`) gets designs with the costs of its case; same-edge key pairs; realizable-key sum `<= K_d` per labelled plaquette; exact bulk clique value |
| [gd_table.py](gd_table.py) | the closed form: counts by formula and by enumeration (`d = 2..7`), exact brackets, the octahedral table |
| [gd_coins.py](gd_coins.py) | axis-star coin counts `q_c`, `N_c` by exact enumeration and the coin table (`d <= 12`) |
| [audit_extra.py](audit_extra.py) | sanity: own cliques versus `networkx.find_cliques`, the anchor-variant diagnosis, type-N designs contain `0`, mutation tests |
| [e2e_check.py](e2e_check.py) | sanity: literal end-to-end test of the surgery map on random states |
| [params/lp_3_cube_rho_256.json](params/lp_3_cube_rho_256.json) | frozen rule table `T_3` (SHA-256 `c43d5de37ae33e1f71156290939dd2555e62e1704ad1a81bcc35cc4053f8d149`; 168 B and 168 N configurations, 813 options) |
| [params/lp_2_box_rho.json](params/lp_2_box_rho.json) | frozen rule table `T_2` (SHA-256 `8e23d624fe8fe11f152eda6d590d185da0988ae3374c06578a32f467f23e93fa`; 112 B and 128 N configurations, 339 options) |
| [expected_outputs/](expected_outputs/) | the outputs of the reference runs, compared line by line (`e2e/`: the recorded end-to-end logs), and the JSON summary of the reference full reproduction |

`cr_check.py`, `oct_check.py`, `gd_table.py`, `gd_coins.py`, `audit_extra.py`
and `e2e_check.py` were written from the definitions independently of the
program that produced the rule tables (only the JSON tables are read).
`cr_check.main()` also returns its exact values (and tracks the worst type-S
key over all `L`) for `check_window_general.py`, `gd_table.py` encloses the
logarithm in `beta_*` with a rational atanh series (the printed floors are the
same as with mpmath intervals), and `gd_coins.py` stores plaquettes as compact
byte strings. The bracket names of `oct_check.py` are `table` (the exact
bracket of Table 1), `coarse` (a coarser admissible bracket) and `univ`
(`[1/(2d), 7/20]`). In the output of `audit_extra.py`, "one-anchor
bookkeeping" is the bookkeeping of the program that produced the tables.

## Default light check

```sh
python diluted_uniqueness/window_general/check_window_general.py
```

Standard library only; about 12–14 s with CPython 3.13.5, at most about
150 MB per process. It decides, in exact rational arithmetic:

1. both rule tables have the recorded SHA-256;
2. `cr_check.py` on `T_3` and `T_2` (type S at `L = 5, 6, 7`): every option
   satisfies the design hypotheses, the configuration lists are complete,
   the weights are positive and sum to 1, and the exact constants:
   - cube: `q = N_Q = 180`, `q_c = 78`, `N_c = 234`,
     `K = K_B = 2^20.704819` (attained at `p_lo = 211/1000`),
     `K_N = 2^15.2143`, type-S realized-key bound `6 x 2^12.0559 <= 2^14.6409`,
     and the two cube windows and thresholds above;
   - box (`p` in `[1/2 - 10^-6, 1/2 + 10^-6]`): `q = 32`, `N_Q = 64`,
     `q_c = 20`, `N_c = 76`, `K = K_B = 2^19.344193` (attained at `p_hi`),
     `K_N = 2^12.8191`, type-S bound `8 x 2^9.8074 <= 2^12.81`, and the two box
     windows and threshold above;
3. `oct_check.py` for `d = 3` at `L = 5` on both brackets (`coarse` and `table`):
   the rule designs (B 918, N 972, S 14,409 configurations) satisfy the
   hypotheses with the costs of their cases; all 35,778 same-edge key pairs
   are incompatible except N-int pairs at different terminals; the
   realizable-key sum per labelled plaquette is at most `K_d` (ratio exactly 1
   at the worst); the exact bulk clique value is at most `K_d`;
4. the closed form: `q(d)`, `N_Q(d)` by formula and by enumeration for
   `d = 2..7` (Lemma 5.25); every displayed entry of Table 1
   (`d = 3..8, 10, 20, 100`: exact brackets, `q`, `N_Q`, `log2 K_d`, the
   floors of `Delta_d`, `g_d`, `beta_*(d)` and of `g_d` with the universal
   bracket `[1/(2d), 7/20]`, the finer-bracket floors, and
   `1.0e-5 <= d^12 g_d <= 2.5e-4` on the tabulated range); every entry of the
   axis-star coin table (Remark 5.27) (`d = 3..8, 10, 12`: `q_c`, `N_c`, floors, gains).
   A displayed floor `x` of an exact value `v` is checked as
   `x <= v < x + one unit of the last digit`.

The printed outputs of the four programs are also compared with
`expected_outputs/`; that comparison is a regression check against the
reference runs, the decisions are the assertions on recomputed values. The run
ends with a JSON summary equal to
[../../expected/diluted_window_general.json](../../expected/diluted_window_general.json).

## Full reproduction

```sh
python diluted_uniqueness/window_general/reproduce_window_general.py --jobs 6 --e2e
```

Needs networkx for `audit_extra.py`. It runs the default check and then, in
parallel child processes, compares every output line with
`expected_outputs/`:

- `oct_check.py` for `d = 3` at `L = 6` (bracket `table`; S 22,977), `d = 4` at
  `L = 5` (brackets `table` and `coarse`; S 383,076; about 45 s and 3.8 GB each),
  `d = 5` at `L = 5` (bracket `table`, without the bulk clique; S 7,811,935;
  about 6.5 min and 0.4 GB): validity of the rule with the costs of each
  case, 880,808 (`d = 4`) and 1,848,890 (`d = 5`) same-edge key pairs, the
  realizable-key bound, and on the stated brackets the bulk clique values
  `2^27.100618 <= K_4 = 2^27.384688` (`table`), `2^29.003687 <= 2^29.220844`
  (`coarse`), `2^29.862414 <= K_5 = 2^30.205576`;
- `audit_extra.py` on both tables (sanity, not a proof input);
- with `--e2e`, the end-to-end sanity runs whose parameters were recorded
  (not a proof input): `d = 2` at `L = 5, 6, 7` (600 states each; seeds 13,
  11, 12; `p = 0.35, 0.40, 0.45`; `s = 0.3, 0.5, 1.0`) and `d = 3` at `L = 5`
  (seed 21, `p = 0.20`, `s = 0.5`, first 180 states). Other `d = 3` sanity
  runs, whose `p` and `s` were not recorded, are not repeated.

Reference run (clean copy of the repository, x86-64 Linux, CPython 3.13.5,
6 parallel jobs, shared machine): 18 min 19 s wall and 26 CPU-minutes with
`--e2e`; the `d = 3` end-to-end run alone takes 18 min on one core, and
without `--e2e` the run takes about 6.5 min (the `d = 5` check). Peak memory
3.8 GB per `d = 4` process (two run in parallel with `--jobs >= 2`). Every
output equals the recorded one, including the four end-to-end logs. Its
summary is
[expected_outputs/reproduce_window_general.json](expected_outputs/reproduce_window_general.json)
(its labels follow the current bracket names; the recorded values and
times are those of the reference run).

## Relation to the proof

The manuscript proves the design lemma (Lemma 5.6), the counting lemmas
(Section 5.3: exact preimage sum, compatibility, arc reconstruction,
log-convexity on the bracket), the curve (Section 5.4), the coin partitions
(Section 5.5), the octahedral rule's validity, costs and key multiplicity, and
the covering counts (Section 5.6). The programs certify the finite
statements about the explicit tables `T_3`, `T_2` and the type-S rule, and
evaluate the constants exactly. The octahedral closed form is a hand proof;
`oct_check.py` confirms it exhaustively for `d <= 5`, and the table entries
are evaluated exactly.

Inputs from the literature (cited in the manuscript, not re-derived):
`mu_3 <= 4.7387`, `mu_4 <= 6.8040` (Pönitz–Tittmann), Hammersley's
`p_c >= 1/mu`, `p_c >= 1/(2d-1)`, the Gomes–Pereira–Sanchis upper bound,
`p_c(Z^3) <= 2 sin(pi/18)`, and Kesten's `p_c(Z^2) = 1/2`.

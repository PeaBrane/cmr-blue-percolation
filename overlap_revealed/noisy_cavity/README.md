# Noisy-cavity local inputs and the rows of Table 3 in dimension 9

This directory holds the sharpened local inputs of the overlap-revealed route (Section 4.4 of the first manuscript)
and the certificate of the five dimension-9 rows of Table 3 built on them. The proof keeps the oriented engine of
`../certify.py` (Theorem 4.33) unchanged and replaces two local inputs:

- **Lemma 4.15** (noisy-cavity reduction). For every cavity law, the cavity factor `2cosh(beta y)` may be replaced
  by `2cosh(beta_c' y)`, for any `beta_c'` in `[atanh(t b), beta]` with `b = tanh((2d-1) beta)`. The effect is
  a star-independent constant and a flip-invariant law. The certificate stores `wc' = e^{2 beta_c'}`.
- **Lemma 4.19** (pair certificate for the floor). Rationals `psi_U(x)` satisfying the chord conditions (31) and
  the class sums `Psi(j) <= 0` of (32) imply that every conditional bond probability is at least `p`.
- **Lemma 4.20** (symmetrized odds certificate) and **Corollary 4.21** (sharpened Holley line). If every two-block
  class has `Gamma(l+, l-) >= 0`, then the Holley line `Phi+/Phi- >= g_h g_K^S` holds in all `2d+1`
  environments.

The finite conditions are (N1) and (N2) of Section 4.7, together with (C3)–(C6). The same deciders also certify the
noisy-cavity rows of the macrostep route (Table 4, d = 8 and 7), in `../../macrostep/`.

## Rows

| t = tanh β | p (certified floor) | line (g_K, g_h) | score used | θ\* ≥ |
|---|---|---|---|---|
| 3/20 | 2689/10000 | (1038887/10^6, 831941/10^6) | Score_c ≤ 0.9906982 | 0.003173 |
| 29/200 | 1309/5000 | (2073/2000, 53049/62500) | Score_c ≤ 0.9945159 | 0.002004 |
| 31/200 | 2759/10000 | (1041383/10^6, 101831/125000) | Score_c ≤ 0.9910479 | 0.002804 |
| 4/25 | 707/2500 | (521969/500000, 99531/125000) | Score_c ≤ 0.9967978 | 0.000899 |
| 7/50 | 159/625 | (1034171/10^6, 864557/10^6) | Score′_c ≤ 0.9975456 (Lemma 4.28) | 0.000950 |

In `params.json` and in the printed JSON, the field `theorem` of a row names its row of Table 3 (`Table3-row1` to
`Table3-row5`, in the printed order of `t = 7/50, 29/200, 3/20, 31/200, 4/25`).

The global entries of Table 3 are the worse of two independent implementations, rounded in the safe direction: the
first implementation (`../crosscheck/`) and the second (`../certify.py`). For Score′_c the two implementations of
Lemma 4.28 are `../crosscheck/covariance.py` and `lemma36.py`.

## Files

| Path | Content |
|---|---|
| `certify_noisy_cavity.py` | **Default checker** (standard library). |
| `local_certificates.py` | Exact deciders for Lemmas 4.15, 4.19 and 4.20 (standard library). It decides the chord conditions by Bernstein subdivision and Sturm sequences, and computes the class sums and class counts from binomial formulas. |
| `lemma36.py` | Second implementation of Lemma 4.28 (`c'_cov`, `rho'_c`) and the least torus size `L_*` of its variant for any `2dK' < 1` (standard library). |
| `params.json` | The five rows `(t, p, g_K, g_h, xbar, score used)` and 155 asserted numbers, each with its locator in the manuscript or marked as a reference value. |
| `data/CERT_DATA_d9.json` | The certificate rationals: `wc'`, the 162 values `psi_U(x)`, the Holley tangent points and `Lambda'_k`. SHA-256 `6637d2a4594c94de42bb623e377451c8bd431f915e5e7a37f0655a7892e96edd`. |
| `crosscheck/` | Optional independent checks (below). |

## Default verification

```sh
python overlap_revealed/noisy_cavity/certify_noisy_cavity.py
```

For each of the five rows the checker decides, in exact rational arithmetic:

- **(N)** the noisy-class constraint `e^{2 atanh(t b)} <= wc' <= e^{2 beta}`;
- **(P)** the 162 point conditions, the 1377 chord polynomials of degree 5, and the 18 class sums `Psi(j) <= 0`;
- **(S)** that all 1330 two-block classes have `Gamma >= 0` in the 19 environments, with the stored `Lambda'_k`
  recomputed and the flip symmetry checked;
- **(W)** plain aligned-frozen witness `< p <=` noisy aligned-frozen witness, so the noisy class is strictly smaller
  than the plain class (Remark 4.17);
- **(H)** the hypotheses (C3) and the mean-field root test (C4) at the frozen `xbar`;
- **(G)** `rho_-`, `eta'`, Score and Score_c with `rho_-^2 >= c_cov/4`, using the functions of `../certify.py`;
- **(G′)** `c'_cov`, `rho'_c` and Score′_c from `lemma36.py`, and the condition `2d Kb <= 35/100` of Lemma 4.28
  where that lemma is used;
- **(T)** θ\* and θ′.

It then asserts the 155 listed numbers. The printed values are those of the second implementation and of
`lemma36.py`. They are slightly sharper than the worse-of-two entries of Table 3, for example Score_c ≤ 0.9906968
against the displayed 0.9906982.

One asserted value differs from the first implementation's. At `t = 3/20` the first implementation gives
`9 rho'_c p >= 1.0675481`; its own `mbar` gives a `rho_-` larger by about `4e-8`. The second implementation gives
`1.06754804...`, so `params.json` asserts the safe 7-digit value `1.0675480`. `crosscheck/crosscheck_global_A.py`
asserts the value `1.0675481` against the first implementation.

The JSON printed at the end equals `../../expected/noisy_cavity.json`, including `value_vector_sha256`, a digest
of the exact values. With CPython 3.13.5 on an Apple Silicon laptop the run takes about 14 s and uses about 50 MB.

## Crosschecks (optional)

Run from the repository root. Times are for the same laptop.

| Program | Needs | Scope | Time |
|---|---|---|---|
| `crosscheck/crosscheck_global_A.py` | mpmath | The first implementation (200-bit intervals, `tbar = K'`, its own `mbar`) at all five rows. It uses Lemma 4.28 or, at `t = 31/200` and `4/25`, its variant for any `2dK' < 1` (`crosscheck/covariance_36prime.py`; its `L_*` is 39 and 41, one more than the least values). It asserts every worse-of-two display, the first implementation's own reference values and its components at `t = 3/20`. | 31 s |
| `crosscheck/verify_pair_bruteforce.py` | NumPy | Lemma 4.19 by an implementation written from the statement (the second program of Section 4.7). It uses its own Sturm sequences and computes the class sums (32) by brute-force enumeration of all `2^17` sign vectors. | 38 s |
| `crosscheck/verify_holley_bruteforce.py` | NumPy | Lemma 4.20 with brute-force block enumeration and an integer convolution. | 57 s |
| `crosscheck/lemmaN_exact.py` | NumPy | Lemma 4.15 by exact rational enumeration on the torus `T_4^2` and on two synthetic `m = 6` stars. It also checks the geometry for `d <= 9`, `L = 4..7`, the `L = 3` counterexample and a negative control. | 2 s |
| `crosscheck/nesting_exact.py` | — | The nesting of the noisy classes (Remark 4.17) exactly in 108 instances; strictness at the five rows. | < 1 s |
| `crosscheck/pairing_exact.py` | NumPy | Pairing identity in the proof of Lemma 4.19 in 24 small instances, and the full chain of Lemma 4.19 at d = 9 for 27 laws per temperature (`t = 3/20` and `31/200`). | 105 s, about 210 MB |

The floating-point linear programs that chose the certificate data are not included; they play no role in the
verification.

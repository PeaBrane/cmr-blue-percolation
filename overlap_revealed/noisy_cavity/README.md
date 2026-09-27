# Noisy-cavity local inputs and the second d = 9 proof

This directory holds the sharpened local inputs of the overlap-revealed route and the second proof of dimension 9
built on them. The proof keeps the audited oriented engine of `../certify.py` (global.md Thm 4.5, code base B)
unchanged and replaces two local inputs:

- **Lemma N** (noisy-cavity reduction). For every cavity law, the cavity factor `2cosh(beta y)` may be replaced by
  `2cosh(beta_c' y)`, for any `beta_c'` in `[atanh(t b), beta]` with `b = tanh((2d-1) beta)`. The effect is
  a star-independent constant and a flip-invariant law. The certificate stores `wc' = e^{2 beta_c'}`.
- **Lemma P** (pair certificate for the bond floor). Rationals `psi_U(x)` satisfying the chord conditions (V) and
  the class sums `Psi(j) <= 0` (C) imply that every conditional bond probability is at least `p`.
- **Lemma S** (symmetrised coupled odds certificate). If every two-block class has `Gamma(l+, l-) >= 0`, then the
  Holley line `Phi+/Phi- >= g_h g_K^S` holds in all `2d+1` environments.

Section references (§) below are to the research write-up of the noisy-cavity route, which the manuscript
supersedes: the statements are in §2, the lemmas in §3 and the numbers in §7. The same deciders also certify instances 2 and 3
of the macrostep route (d = 8 and 7), in `../../macrostep/`.

## Theorems

The research notes call these Q9 and Q9′. The macrostep route has its own Q9, at `t = 7/50` with other local
inputs, so the manuscripts must keep the two apart.

| label | t = tanh β | p (certified floor) | line (g_K, g_h) | score used | θ\* ≥ |
|---|---|---|---|---|---|
| Q9 | 3/20 | 2689/10000 | (1038887/10^6, 831941/10^6) | Score_c ≤ 0.9906982 | 0.003173 |
| Q9′ | 29/200 | 1309/5000 | (2073/2000, 53049/62500) | Score_c ≤ 0.9945159 | 0.002004 |
| Q9′ | 31/200 | 2759/10000 | (1041383/10^6, 101831/125000) | Score_c ≤ 0.9910479 | 0.002804 |
| Q9′ | 4/25 | 707/2500 | (521969/500000, 99531/125000) | Score_c ≤ 0.9967978 | 0.000899 |
| Q9′ | 7/50 | 159/625 | (1034171/10^6, 864557/10^6) | Score′_c ≤ 0.9975456 (Lemma 3.6) | 0.000950 |

The displayed global numbers are the worse of two research code bases, rounded in the safe direction: code base A
(`../crosscheck/`) and code base B (`../certify.py`). For Score′_c the two implementations are code base A and the
second implementation `lemma36.py`.

## Files

| Path | Content |
|---|---|
| `certify_noisy_cavity.py` | **Default checker** (standard library). |
| `local_certificates.py` | Exact deciders for Lemmas N, P and S (standard library). It decides the chord conditions by Bernstein subdivision and Sturm sequences, and computes the class sums and class counts from binomial formulas. |
| `lemma36.py` | Second implementation of global.md Lemma 3.6 (`c'_cov`, `rho'_c`) and the least torus size `L_*` of Lemma 3.6′ (standard library). |
| `params.json` | The five rows `(t, p, g_K, g_h, xbar, score used)` and 155 displayed numbers with their sources in the research write-up. |
| `data/CERT_DATA_d9.json` | The certificate rationals: `wc'`, the 162 values `psi_U(x)`, the Holley tangent points and `Lambda'_k`. The file is identical to the research file, SHA-256 `6637d2a4594c94de42bb623e377451c8bd431f915e5e7a37f0655a7892e96edd`. |
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
  than the plain class;
- **(H)** the Step-0 hypotheses and the mean-field root test at the frozen `xbar`;
- **(G)** `rho_-`, `eta'`, Score and Score_c with `rho_-^2 >= c_cov/4`, using code base B;
- **(G′)** `c'_cov`, `rho'_c` and Score′_c from `lemma36.py`, and the Lemma 3.6 condition `2d Kb <= 35/100` where
  that lemma is used;
- **(T)** θ\* and θ′.

It then asserts the 155 displayed numbers. The printed values are code base B's and the second implementation's.
They are slightly sharper than the worse-of-two displays, for example Score_c ≤ 0.9906968 against the displayed
0.9906982.

One display differs from the research text. The notes print `9 rho'_c p >= 1.0675481` at `t = 3/20` for code base
A, whose own `mbar` gives a `rho_-` larger by about `4e-8`. The second implementation gives `1.06754804...`, so
`params.json` asserts the safe 7-digit value `1.0675480`. `crosscheck/crosscheck_global_A.py` asserts the
research value against code base A.

The JSON printed at the end equals `../../expected/noisy_cavity.json`, including `value_vector_sha256`, a digest
of the exact values. With CPython 3.13.5 on an Apple Silicon laptop the run takes about 14 s and uses about 50 MB.

## Crosschecks (optional)

Run from the repository root. Times are for the same laptop.

| Program | Needs | Scope | Time |
|---|---|---|---|
| `crosscheck/crosscheck_global_A.py` | mpmath | Code base A (200-bit intervals, `tbar = K'`, its own `mbar`) at all five rows. It uses Lemma 3.6 or, at `t = 31/200` and `4/25`, Lemma 3.6′ (`crosscheck/covariance_36prime.py`; its `L_*` is 39 and 41, one more than the least values). It asserts every worse-of-two display, the code-base-A columns of §7.2 and the code-base-A components at `t = 3/20`. | 31 s |
| `crosscheck/verify_pair_bruteforce.py` | NumPy | Lemma P by an implementation written from the statement. It uses its own Sturm sequences and computes the class sums (C) by brute-force enumeration of all `2^17` sign vectors (research script `r1_verify_pair.py`). | 38 s |
| `crosscheck/verify_holley_bruteforce.py` | NumPy | Lemma S with brute-force block enumeration and an integer convolution (research script `r1_verify_holley.py`). | 57 s |
| `crosscheck/lemmaN_exact.py` | NumPy | Lemma N by exact rational enumeration on the torus `T_4^2` and on two synthetic `m = 6` stars. It also checks the geometry for `d <= 9`, `L = 4..7`, the `L = 3` counterexample and a negative control. | 2 s |
| `crosscheck/nesting_exact.py` | — | Lemma N′ (nesting of the noisy classes) exactly in 108 instances; strictness at the five rows. | < 1 s |
| `crosscheck/pairing_exact.py` | NumPy | Pairing identity of Lemma P in 24 small instances, and the full chain of Lemma P at d = 9 for 27 laws per temperature (`t = 3/20` and `31/200`). | 105 s, about 210 MB |

The first implementation of the local certificates (the primary checker `exact_pair.py`/`exact_local.py` of code base L1, run in
a SciPy environment),
the second implementation of code base L1 (`verify_*_indep.py`), and the floating-point LPs that chose the data remain in the
research notes. They are not needed here.

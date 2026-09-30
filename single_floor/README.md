# Single-floor oriented comparison: dimensions 16–20, 22 and 25

These programs certify the finite inequalities of the single-floor oriented
comparison in Appendix C of the first manuscript ("The fresh-star route: star
estimates and a single-floor comparison"). The appendix reproves the
dimension-22 theorem, gives dimensions 16–20, and gives a closed-form result
in dimension 25.

**Theorem C.2 (single-floor oriented comparison).** Let `p` be a uniform
fresh-star floor. If `p > F_d := 1 - 1/G_d`, where `G_d` is the collision
Green function of two independent uniform oriented walks, then every selected
periodic joint limit satisfies

    nu(o <-> infinity blue, q_o = s) >= (1/2) (p - F_d) / (p (1 - F_d))   (s = +1, -1),

and both overlap signs have infinite blue components almost surely. The row
`h = 0` of the fresh-star floor, `p_0(2d, beta)`, is such a floor.

The results used are Theorem C.2, Corollary C.3 (dimensions 16 to 20 and
22), Lemma C.5 (collision Green function bounds), Lemma C.7 (the balanced tilt
maximizes), Lemma C.8 (the row `h = 0` as a single ratio), Lemma C.9
(monotone pieces), Proposition C.10 (reduction to `k = 1` by blocks),
Remark C.11 (signed refinements), Proposition C.13 (a closed-form star floor)
and Corollary C.14 (a closed-form certificate in dimensions 25 and 26).

## Main certificate

[verify_single_floor.py](verify_single_floor.py) uses only the standard
library and runs in the default `run_all.py` (about 7 CPU seconds). Its
checks:

| Manuscript result | Exact check |
|---|---|
| Lemma C.5 ("Collision Green function bounds") | `u_0, ..., u_3` in closed form; `u_n <= bar m_n` and `bar m_n` nonincreasing on the head; exact head sums over `n < 4d` plus the maximal-atom tail of part (e), which determine `F_d` to within `10^-8`; the eight-term closed form `hat G_d` of part (f); the bound `G_22 <= 1.0602931005 < 1.0602932` of part (g), whose tail term is below `8.8e-8` |
| Lemma C.8 ("The row h = 0 as a single ratio") | the closed form `p_* = sinh 2beta / Gamma_1`, with `E_1` in balanced-pair form, equals the nonnegative-residual value `p(1)`; `Gamma_1` lies in `[3.020, 3.093]` |
| Proposition C.10 ("Reduction to k = 1 by blocks") | the frozen block covers of `{2, ..., 2d}` satisfy `B(k_a, k_b) <= Gamma_1`, so `p(k) >= p(1)` for every overlap count and `p_0(2d, beta) = p_*`; the block slacks and sum counts of Table 7 |
| Corollary C.3 ("Dimensions 16 to 20 and 22") | `p_* > 1 - 1/hat G_d >= F_d` at the seven points below, with every entry of Table 6; `p(2)/p(1)` in `[1.0094, 1.0149]` |
| Remark C.11 ("Signed refinements"), `d = 22` | signed residuals `M_(0,1) <= -0.000374848573757` (44 orientations) and `M_(0,2) <= -0.019889903399663` (86 orientations) at `lambda = 1/0.0586`, `c = 119/100`; nonnegative residual `min_(3<=k<=44) p(k) >= 0.059010354650604`, while `p(2) = 0.0583248...` alone is below `0.0586`; `F_22 < 0.0568646 < 0.0586` and root constant at least `0.0157` |
| Remark C.11 ("Signed refinements"), `d = 16` | signed `M_(0,1) < -0.00136` (32 orientations) at `lambda = 10000/689`, `c = 6/5`, and `min_(k>=2) p(k) >= 0.0690100`; the floor `0.0689` exceeds `1 - 1/hat G_16` by a factor of at least `1.0214` |
| Proposition C.13 ("A closed-form star floor"), Corollary C.14 ("A closed-form certificate in dimensions 25 and 26") | `min(p_pair(1), p_pair(2d)) = p_pair(1)` exceeds `1 - 1/hat G_d` at `d = 25` and `26`, with no finite sums; also at `d = 28`, a reference value that the manuscript does not state |

| d | t | `p_0 = p_* >=` | `F_d <=` (`hat G_d`) | ratio (`hat G_d`) `>=` | θ\* `>=` | ratio (exact head) `>=` | smallest block slack `>=` |
|---|---|---|---|---|---|---|---|
| 16 | 13/125 | 0.0679971076 | 0.067455 | 1.0080 | 0.0042 | 1.0146 | 0.0267 |
| 16 | 51/500 | 0.0680421165 | 0.067455 | 1.0087 | 0.0046 | 1.0153 | 0.0065 |
| 17 | 1/10 | 0.0659393703 | 0.063115 | 1.0447 | 0.0228 | 1.0502 | 0.0073 |
| 18 | 12/125 | 0.0640504639 | 0.059313 | 1.0798 | 0.0393 | 1.0845 | 0.0052 |
| 19 | 7/75 | 0.0623013944 | 0.055951 | 1.1135 | 0.0539 | 1.1174 | 0.0064 |
| 20 | 1/11 | 0.0606881634 | 0.052955 | 1.1460 | 0.0672 | 1.1494 | 0.0044 |
| 22 | 11/125 | 0.0577808419 | 0.047843 | 1.2077 | 0.0903 | 1.2103 | 0.0099 |

Here θ\* `= (1/2)(p_* - F_d)/(p_* (1 - F_d))` with `F_d <= 1 - 1/hat G_d`,
as in Table 6. With the exact head bound on `F_d` it is at
least `0.0077`, `0.0080`, `0.0255`, `0.0414`, `0.0556`, `0.0686` and `0.0912`
in the same row order.

The certificate-free dimension is 25, at `t = 27/400`: `p_pair(1) >= 0.04196433`,
`p_pair(50) >= 0.05133661`, `F_25 <= 0.04180338`, ratio at least `1.0038` and
θ\* at least `0.0020`. At `d = 26`, `t = 33/500` the values are `0.04113940`,
`0.05037252` and `0.04011778`, with ratio at least `1.0254` and θ\* at least
`0.0129`; at `d = 28`, `t = 8/125` the ratio is at least `1.0673`.

Every displayed number above is asserted against the exact rationals. The run
ends by printing a JSON summary equal to
[../expected/single_floor.json](../expected/single_floor.json), including a
SHA-256 digest of the ordered exact values.

## Supplementary checks

Both programs use only the standard library and run with
`run_all.py --supplementary`.

- [check_lemmas.py](check_lemmas.py), about 30 CPU seconds:
  - Lemma C.7 ("The balanced tilt maximizes") at all 256 `(d, t, k)` triples,
    with the complete row `h = 0` evaluated over every cavity orientation
    index; each row is minimised at `k = 1`, independently of the block
    reduction (the manuscript's "independent check" of the full row);
  - the single-ratio identity of Lemma C.8 for every `k`, and Lemma C.9
    ("Monotone pieces") (i)–(v), at `(16, 13/125)` and `(17, 1/10)`;
  - the signed count enumeration against brute force of the definition in 18
    cases, and against a separate `k = 1` formula;
  - at `d = 22`, `M_(0,1)` and `M_(0,2)` equal, as exact rationals, the
    maxima of the dimension-22 certificate
    [../verify_signed_star.py](../verify_signed_star.py).
- [check_side_results.py](check_side_results.py), about 25 CPU seconds,
  covers statements that no corollary uses:
  - the minimum over `k` of the nonnegative-residual floor is not at `k = 1`
    uniformly in `beta`: exact witnesses at
    `(Delta, t) = (16, 1/2), (20, 3/10), (20, 1/3), (32, 1/5), (44, 3/20)`,
    and seven controls;
  - a same-cavity counterexample at `e^beta = 3` and its full-support
    perturbations;
  - Remark C.12 ("Where the route stops"): `F_15` in
    `[0.0718626355, 0.0718626462]`, and at `t = 27/250` the signed criterion
    at `k = 2` (`lambda = 1/0.0725`, `c = 6/5`, 58 orientations) gives
    `M_(0,2) < -0.0105`, while the nonnegative residual gives
    `min_(k>=3) p(k) >= 0.072807`. The remark's statement about `k = 1` rests
    on a floating-point scan and is not checked here.

## Relation to the proof

The appendix supplies the single-floor theorem, the collision Green function
bounds, the single-ratio identity, the two monotonicity lemmas, the block
reduction, the closed-form floor and the transfer to limits. It also proves
the fresh-star floor and its signed criterion, the star estimate of Appendix C
of the first version (arXiv:2609.22301v1). These programs certify only the finite rational
inequalities that those statements reduce to. Runtimes were measured with
CPython 3.13.5 on Apple Silicon; each program uses under 50 MB of memory.

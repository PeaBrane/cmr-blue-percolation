# Exact computational certificates for CMR blue percolation

Companion code for Yan Ru Pei, *Towards the spin-glass transition in finite
dimensions via blue percolation* (2026).

The main certificate verifies the finite inequalities used to prove blue
percolation in dimension **22**, at `tanh(beta) = 11/125`, for iid fair unit
couplings. A second certificate checks the paper's simultaneous percolation
and finite spin-glass susceptibility example in dimension **180**.

This repository is private during manuscript preparation. The intended public
location is <https://github.com/PeaBrane/cmr-blue-percolation>, to be made public
with the paper. The manuscript cites a specific commit so that its computational
inputs remain identifiable after subsequent development.

## Reproduce

The programs need Python 3.10 or later and **only the standard library**.
They were validated with CPython 3.13.5. No Peapods installation, downloaded
input data, numerical libraries or network access is needed to run them.

```sh
uv venv --python 3.13
.venv/bin/python run_all.py
```

To include all independent finite-model checks:

```sh
.venv/bin/python run_all.py --supplementary
```

An existing local virtual environment with a supported Python also works.
The runner uses that same interpreter for every child process. All entrypoints
reject optimized Python execution, because assertions implement the decisive
comparisons. Successful runs end with `PASS: 3 verification scripts completed.`
or `PASS: 6 verification scripts completed.` respectively.

## What is checked

| Program | Mathematical scope |
|---|---|
| [verify_signed_star.py](verify_signed_star.py) | Complete d22 parameter certificate: five signed-residual cases with 299 reduced cavity orientations, 205 other overlap/failure cases, five failure envelopes, all 63 proper triangular direction subsets, positive applicable increments, and the strict triangular critical inequality. |
| [verify_signed_star_independent.py](verify_signed_star_independent.py) | Separately implemented direct rational enumeration of the five signed cases. Its accumulation uses absolute cavity fields. It does not check the other 205 cases or global comparison. |
| [verify_physical_separation.py](verify_physical_separation.py) | d180 at `tanh(beta)=1/56`: the MNS site parameter exceeds `23/2000` and the cited site threshold; `kappa < 5/8` gives `chi_SG < 8/3`. |
| [checks/smallstar_audit.py](checks/smallstar_audit.py) | Central-spin and activation sums, the nonnegative residual coefficient and orientation reduction on degrees 1–7 at three rational temperatures. |
| [checks/numerator_audit.py](checks/numerator_audit.py) | Independent tilted-numerator, pointwise maximum and weighted-Jensen checks, including rational cavity mixtures, on degrees 1–7. |
| [checks/check_susceptibility.py](checks/check_susceptibility.py) | 2,916 one-edge identities/inequalities and complete disorder/spin sums on four small graphs; 3,136 odd-set-pair bounds. |

The primary certificate uses integers and `fractions.Fraction`. The signed
cases use exact kernels. Other d22 cases use outward rational enclosures on
a `2^-80` grid, with signs determining the direction of rounding. Decimal
values in output are summaries; no floating-point comparison proves an assertion.

The expected d22 score is approximately **1.000816029521908**, strictly larger
than one by exact rational comparison. [expected/dimension22.json](expected/dimension22.json)
records the exact score, selected parameters and output-vector digest:

```text
d7f1367cac657b8821c60093f30aac8932e2941175cb4dceac57b0450da49f7d
```

This digest hashes the ordered rational bound vector, not the program source.
The Git commit identifies the source. The verifier recomputes the inequalities
from raw parameters and does not read the expected-results file as an oracle.
Elapsed times in some independent-check outputs are not reproducibility targets.

## Relation to the proof

Section 4 of the manuscript supplies the physical exploration and cumulative
comparison. Appendix A defines the exact normalized masses, proves the signed
and nonnegative residual bounds, and reduces their maxima to the finite sums
checked here. Appendix C gives the d180 integer inequalities. The programs
certify those finite computations; the analytic conditioning, probability and
infinite-volume arguments are supplied by the manuscript. The supplementary
small-graph checks are independent finite evidence for identities used there.

The signed refinement is used only for `(h,k) = (0,1), (0,2), (1,1), (1,2),
(2,1)`, with tangent point `119/100` and `lambda_h = 1/p_h`. Every other
allowed overlap count uses the nonnegative residual. The full-star quenched
normalization retains the signs on failed incident edges in both calculations.

The repository contains the computational dependencies of the paper's two
explicit certificates. Historical dimension searches, exploratory geometry,
private worklogs and reference PDFs are maintained separately.

## Cite

Use [CITATION.cff](CITATION.cff) and include the commit hash recorded by the
manuscript. This is a computational companion to the paper, not a formal
proof-assistant development.

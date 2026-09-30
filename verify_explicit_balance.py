"""Exact finite-susceptibility check at the seventeen explicit pairs (Corollary 5.5, Table 8).

For unit couplings F(2 beta) = cosh(2 beta) = (1 + t^2) / (1 - t^2) with t = tanh(beta).
Proposition 5.1 needs F^(2d-1) < 3 and bounds chi_SG by
chi_* = 1 + (F^(2d) - 1) / (3 - F^(2d-1)). The program asserts F^(2d-1) < 12/5 and
chi_* < 7/2 at every pair of Table 1, chi_* < 7/4 at (7, 3/20), and every entry of
Table 8 (values rounded upward to four decimals).
"""

from fractions import Fraction as Fr
from math import ceil
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

# (d, t): (F^(2d-1), chi_*) as printed in Table 8
TABLE = {
    (12, "3/25"): ("1.9396", "1.9394"), (12, "11/100"): ("1.7448", "1.6274"),
    (12, "23/200"): ("1.8375", "1.7628"), (12, "1/8"): ("2.0520", "2.1784"),
    (12, "13/100"): ("2.1760", "2.5179"), (11, "13/100"): ("2.0338", "2.1422"),
    (11, "3/25"): ("1.8310", "1.7566"), (10, "27/200"): ("1.9990", "2.0721"),
    (10, "13/100"): ("1.9008", "1.8790"), (9, "7/50"): ("1.9474", "1.9740"),
    (9, "29/200"): ("2.0441", "2.1842"), (9, "3/20"): ("2.1493", "2.4673"),
    (9, "31/200"): ("2.2638", "2.8678"), (9, "4/25"): ("2.3884", "3.4749"),
    (8, "3/20"): ("1.9643", "2.0183"), (7, "31/200"): ("1.8679", "1.8478"),
    (7, "3/20"): ("1.7952", "1.7286"),
}


def up4(x: Fr) -> str:
    return f"{ceil(x * 10**4) / 10**4:.4f}"


def explicit_balance_certificate():
    chis = {}
    for (d, t_str), printed in TABLE.items():
        t = Fr(t_str)
        f = (1 + t * t) / (1 - t * t)
        y = f ** (2 * d - 1)
        chi = 1 + (f ** (2 * d) - 1) / (3 - y)
        assert y < Fr(12, 5) and chi < Fr(7, 2)
        assert (up4(y), up4(chi)) == printed, (d, t_str)
        chis[d, t] = chi
        print(f"d={d:2d}  t={t_str:7s}  F^(2d-1)<={up4(y)}  chi_*<={up4(chi)}")
    assert chis[7, Fr(3, 20)] < Fr(7, 4)
    print("PASS explicit balance at the 17 pairs of Table 1")
    print("  F^(2d-1) < 12/5 and chi_* < 7/2 everywhere; chi_* < 7/4 at (7, 3/20)")


if __name__ == "__main__":
    explicit_balance_certificate()

"""Exact d180 blue-percolation and physical-susceptibility certificate."""

from fractions import Fraction as F
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")


def physical_separation_certificate():
    t, degree = F(1, 56), 360
    L = 2 * t / (1 + t) ** 2
    R = (1 + 6 * t**2 + t**4) / (1 - t**4)
    assert L == F(112, 3249)
    assert R == F(9853313, 9834495)
    p = L / (1 + R**degree)
    assert 149273 * 9834495**360 > 74727 * 9853313**360
    assert p > F(23, 2000)
    assert 2 * 1977**60 < 2000**60
    A = (1 + 3 * t**2) / (1 - t**2) ** 2
    B = t**2 * (3 + t**2) / (1 - t**2) ** 2
    assert A == F(9843904, 9828225)
    assert B == F(9409, 9828225)
    kappa = degree * B * A ** (degree - 1)
    assert 8 * 360 * 9409 * 9843904**359 < 5 * 9828225**360
    assert kappa < F(5, 8)
    assert 1 / (1 - kappa) < F(8, 3)
    print("PASS physical separation d=180, tanh(beta)=1/56")
    print("  MNS site parameter > 23/2000 > rigorous site threshold")
    print("  kappa < 5/8; physical spin-glass susceptibility < 8/3")


if __name__ == "__main__":
    physical_separation_certificate()

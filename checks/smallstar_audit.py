"""Independent exhaustive rational Appendix A audit; stdlib only.

This file is derived from central Ising spins and the displayed formulas,
without reading the manuscript's existing verification program.
"""

from fractions import Fraction as F
from functools import cache
from itertools import product
from math import lcm
import json
import time
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")


START = time.monotonic()
COUNTS = dict(parameter_cases=0, star_cases=0, sign_weight_checks=0,
              central_spin_assignments=0, lambda_checks=0,
              residual_orientation_checks=0, residual_epsilon_terms=0,
              residual_maximum_checks=0)


def scaled_integers(values):
    denominator = lcm(*(q.denominator for q in values))
    return denominator, [q.numerator * (denominator // q.denominator)
                         for q in values]


for eta in (F(11, 13), F(3, 5), F(1, 3)):
    COUNTS['parameter_cases'] += 1
    alpha = 1 - eta**2
    big_c = (eta + 1 / eta) / 2
    v = (1 - eta**2) / (1 + eta**2)
    u = (1 + v) / 2

    @cache
    def psi(z):
        a = eta**abs(z)
        return 4 * a / (1 + a)**2

    @cache
    def tilted_expectation(n, ell, j, offset):
        # Enumerate each Bernoulli sign directly; no coefficient polynomial.
        ans = F(0)
        means = (v,) * j + (-v,) * (ell - j) + (F(0),) * (n - ell)
        for signs in product((-1, 1), repeat=n):
            probability = F(1)
            for sign, mean in zip(signs, means):
                probability *= (1 + sign * mean) / 2
            ans += probability * psi(offset + sum(signs))
        return ans

    for degree in range(1, 8):
        signs = tuple(product((-1, 1), repeat=degree))
        sums = tuple(range(-degree, degree + 1, 2))
        psi_den, psi_nums = scaled_integers([psi(z) for z in sums])
        kernel_nums = {z: a for z, a in zip(sums, psi_nums)}
        matrix = [[kernel_nums[sum(a * b for a, b in zip(eps, x))]
                   for eps in signs] for x in signs]

        for h in range(degree):
            for k in range(1, degree - h + 1):
                COUNTS['star_cases'] += 1
                s = degree - h - k
                failure_indices = range(k, k + h)
                overlaps = (1,) * (k + h) + (-1,) * s
                weights = []
                numerators = []
                for eps in signs:
                    # Four independent uniform cavity central replicas.
                    direct_w = F(0)
                    direct_a = F(0)
                    for center_sigma, center_tau in product((-1, 1), repeat=2):
                        exponent_twice = -sum(e * (center_sigma + center_tau * q)
                                              for e, q in zip(eps, overlaps))
                        assert exponent_twice % 2 == 0
                        mass = eta**(exponent_twice // 2) / 4
                        satisfies = [center_sigma * e == 1
                                     and center_tau * q * e == 1
                                     for e, q in zip(eps, overlaps)]
                        for i in failure_indices:
                            mass *= 1 - alpha * int(satisfies[i])
                        direct_w += mass
                        direct_a += mass * alpha * int(satisfies[0])
                        COUNTS['central_spin_assignments'] += 1
                    sum_u = sum(eps[:k])
                    sum_v = sum(eps[k + h:])
                    formula_w = (eta**h * (eta**sum_u + eta**(-sum_u))
                                 + eta**sum_v + eta**(-sum_v)) / 4
                    formula_a = eta**h * alpha * eta**(-eps[0] * sum_u) / 4
                    assert direct_w == formula_w, ('w', eta, degree, h, k, eps)
                    assert direct_a == formula_a, ('a', eta, degree, h, k, eps)
                    COUNTS['sign_weight_checks'] += 1
                    weights.append(direct_w)
                    numerators.append(direct_a)

                cs = F(1) if s % 2 == 0 else big_c
                lam = (1 + eta**(2 * k) + 2 * cs * eta**(k - h)) / alpha
                assert min(w / a for w, a in zip(weights, numerators)) == lam
                COUNTS['lambda_checks'] += 1
                residual = [w - lam * a for w, a in zip(weights, numerators)]
                assert min(residual) >= 0
                residual_den, residual_nums = scaled_integers(residual)
                common_den = len(signs) * psi_den * residual_den
                direct_values = []
                for x, kernel_row in zip(signs, matrix):
                    direct = F(sum(r * p for r, p in zip(residual_nums, kernel_row)),
                               common_den)
                    direct_values.append(direct)
                    y = tuple(a * x[0] for a in x)
                    j_u = sum(a == 1 for a in y[1:k])
                    j_v = sum(a == 1 for a in y[k + h:])
                    formula = (big_c**s * tilted_expectation(degree, s, j_v, 0) / 2
                               + eta**h * big_c**k / 2
                               * ((u - lam * v)
                                  * tilted_expectation(degree - 1, k - 1, j_u, 1)
                                  + (1 - u)
                                  * tilted_expectation(degree - 1, k - 1, j_u, -1)))
                    assert direct == formula, ('orientation', eta, degree, h, k, x)
                    COUNTS['residual_orientation_checks'] += 1
                    COUNTS['residual_epsilon_terms'] += len(signs)

                t_s = big_c**s * max(tilted_expectation(degree, s, j, 0)
                                    for j in range(s + 1))
                displayed_m = (t_s / 2 + eta**h * big_c**k / 2
                               * max((u - lam * v)
                                     * tilted_expectation(degree - 1, k - 1, j, 1)
                                     + (1 - u)
                                     * tilted_expectation(degree - 1, k - 1, j, -1)
                                     for j in range(k)))
                assert max(direct_values) == displayed_m, ('maximum', eta, degree, h, k)
                assert displayed_m >= 0
                COUNTS['residual_maximum_checks'] += 1
        print(f'PASS eta={eta} degree={degree}', flush=True)

print(json.dumps(dict(status='PASS', eta_values=['11/13', '3/5', '1/3'],
                      degree_range=[1, 7], h_range='0..degree-1',
                      k_range='1..degree-h', counts=COUNTS,
                      elapsed_seconds=time.monotonic() - START), indent=2))

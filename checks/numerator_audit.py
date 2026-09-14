"""Independent rational audit of Appendix A's paired numerator bound."""

from fractions import Fraction as F
from itertools import product
from math import prod
import json
import time
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")


START = time.monotonic()
COUNTS = dict(parameter_degree_cases=0, pure_coordinate_cases=0,
              raw_numerator_sign_terms=0, tilted_identity_checks=0,
              pointwise_h_checks=0, pure_jensen_checks=0,
              pure_a_bound_checks=0, mixture_cases=0,
              mixture_jensen_checks=0, mixture_a_bound_checks=0)


def weighted_numerator(r_values, w_values, q):
    count = len(r_values)
    wr = sum(w * r for w, r in zip(w_values, r_values)) / count
    inverse = sum(w / r**2 for w, r in zip(w_values, r_values)) / count
    assert inverse >= q**3 / wr**2
    return wr, inverse


for t in (F(1, 12), F(1, 4), F(1, 2)):
    eta = (1 - t) / (1 + t)
    big_c = (1 + t**2) / (1 - t**2)
    v = 2 * t / (1 + t**2)

    for degree in range(1, 8):
        COUNTS['parameter_degree_cases'] += 1
        signs = tuple((1,) + tail for tail in product((-1, 1), repeat=degree - 1))
        # Scale R by c0^degree so every quantity is rational at rational t.
        # eps0=y0=+1 fixes the global sign symmetries used in the manuscript.
        r_by_y = {}
        for y in signs:
            r_by_y[y] = tuple((prod(1 + t * e * a for e, a in zip(eps, y))
                               + prod(1 - t * e * a for e, a in zip(eps, y))) / 2
                              for eps in signs)

        for k in range(1, degree + 1):
            ell = k - 1
            q = big_c**ell
            h = ((1 + t) * (1 + t * v)**ell
                 + (1 - t) * (1 - t * v)**ell) / 2
            w_values = tuple(eta**(-sum(eps[1:k])) for eps in signs)
            assert sum(w_values) / len(signs) == q
            prefactor = t * (1 - t**2)**(degree - 1)
            n_bound = prefactor * q / h**2

            for y, r_values in r_by_y.items():
                COUNTS['pure_coordinate_cases'] += 1
                wr, inverse = weighted_numerator(r_values, w_values, q)
                COUNTS['raw_numerator_sign_terms'] += len(signs)
                COUNTS['pure_jensen_checks'] += 1
                display = ((1 + t) * prod(1 + t * v * a for a in y[1:k])
                           + (1 - t) * prod(1 - t * v * a for a in y[1:k])) / 2
                assert wr / q == display, ('tilted identity', t, degree, k, y)
                assert display <= h, ('H upper bound', t, degree, k, y)
                assert prefactor * inverse >= n_bound, ('A bound', t, degree, k, y)
                COUNTS['tilted_identity_checks'] += 1
                COUNTS['pointwise_h_checks'] += 1
                COUNTS['pure_a_bound_checks'] += 1

            first, last = signs[0], signs[-1]
            alternating = tuple(1 if i % 2 == 0 else -1 for i in range(degree))
            # R is averaged over the cavity distribution before it is inverted.
            mixtures = (
                ((F(1, 3), first), (F(2, 3), last)),
                ((F(1, 7), first), (F(2, 7), alternating), (F(4, 7), last)),
                tuple((F(1, len(signs)), y) for y in signs),
            )
            for mixture in mixtures:
                r_values = tuple(sum(weight * r_by_y[y][i] for weight, y in mixture)
                                 for i in range(len(signs)))
                wr, inverse = weighted_numerator(r_values, w_values, q)
                assert wr / q <= h
                assert prefactor * inverse >= n_bound, ('mixture A bound', t, degree, k)
                COUNTS['mixture_cases'] += 1
                COUNTS['mixture_jensen_checks'] += 1
                COUNTS['mixture_a_bound_checks'] += 1
        print(f'PASS t={t} degree={degree}', flush=True)

print(json.dumps(dict(status='PASS', t_values=['1/12', '1/4', '1/2'],
                      degree_range=[1, 7], k_range='1..degree',
                      cavity_coordinates='all y with y0=+1',
                      deleted_signs='all epsilon with epsilon0=+1',
                      mixture_distributions_per_parameter_degree_k=3,
                      counts=COUNTS, elapsed_seconds=time.monotonic() - START), indent=2))

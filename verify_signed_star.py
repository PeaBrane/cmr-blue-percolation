"""Exact signed-residual certificate for CMR percolation in dimension 22.

The five exceptional local cases use a direct binomial count of every sign
configuration modulo permutation and global-sign symmetries. Other cases use
the previously proved nonnegative-residual finite sums, with outward rounding.
No optimizer output enters any comparison. Only the standard library is used.
"""

from fractions import Fraction as Q
from itertools import product
from math import comb
import hashlib
import json
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")


M = 44
T = Q(11, 125)
TANGENT = Q(119, 100)
TARGETS = tuple(Q(x, 100000) for x in (5860, 5253, 4650, 4080, 3610))
CREDIT = Q(20311, 20000)
GROUPS = (7, 7, 8)
SIGNED_CASES = {(0, 1), (0, 2), (1, 1), (1, 2), (2, 1)}


def convolve(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def binomial_polynomial(n, a, b):
    return [comb(n, j) * a ** (n - j) * b ** j for j in range(n + 1)]


def signed_residual(h, r):
    """Return every relative-cavity-orientation upper bound, exactly.

    Fix epsilon_0=x_0=+1 using the two separate global symmetries. For r<=2,
    enumerate the remaining active signs explicitly. Inactive signs are counted
    in the x=+1 and x=-1 groups. Failed signs remain in the full cavity field.
    """
    assert r in (1, 2)
    eta = (1 - T) / (1 + T)
    s = M - h - r
    lam = 1 / TARGETS[h]
    c = TANGENT
    cosh2 = {z: (eta ** z + eta ** -z) / 2 for z in range(-M, M + 1)}
    cosh1 = {z: cosh2[z // 2] for z in range(-M, M + 1, 2)}
    inverse = {z: 1 / value ** 2 for z, value in cosh1.items()}
    denominator = 2 ** (M - 1)
    bounds = []
    for active_x in product((-1, 1), repeat=r - 1):
        for plus_x in range(s + 1):
            positive = {z: Q(0) for z in inverse}
            negative = {z: Q(0) for z in inverse}
            count = 0
            for active_eps in product((-1, 1), repeat=r - 1):
                active_sum = 1 + sum(active_eps)
                active_field = 1 + sum(x * e for x, e in zip(active_x, active_eps))
                a = eta ** h * (1 - eta ** 2) * eta ** (-active_sum) / 4
                for i in range(plus_x + 1):
                    for j in range(s - plus_x + 1):
                        inactive_sum = 2 * (i + j) - s
                        inactive_field = 2 * i - 2 * j + s - 2 * plus_x
                        w = (eta ** h * cosh2[active_sum] + cosh2[inactive_sum]) / 2
                        residual = w - lam * a
                        binom = comb(plus_x, i) * comb(s - plus_x, j)
                        bucket = positive if residual >= 0 else negative
                        for k in range(h + 1):
                            ways = binom * comb(h, k)
                            field = active_field + inactive_field + 2 * k - h
                            bucket[field] += ways * abs(residual)
                            count += ways
            assert count == denominator
            value = sum(positive[z] * inverse[z] + negative[z] *
                        (2 * cosh1[z] / c ** 3 - 3 / c ** 2) for z in inverse)
            bounds.append(value / denominator)
    assert max(bounds) < 0, (h, r, float(max(bounds)))
    return bounds


def old_residual_cases():
    """Re-evaluate all 205 nonexceptional cases from the original formulas."""
    eta = (1 - T) / (1 + T)
    C = (1 + T * T) / (1 - T * T)
    v = 2 * T / (1 + T * T)
    u = (1 + v) / 2
    scale = 2 ** 80
    kernels = []
    for j in range(M + 1):
        z = eta ** abs(2 * j - M)
        scaled = scale * 4 * z / (1 + z) ** 2
        lo = scaled.numerator // scaled.denominator
        hi = -(-scaled.numerator // scaled.denominator)
        kernels.append((lo, hi))
    tn, td = T.numerator, T.denominator
    polynomials = [binomial_polynomial(n, (td - tn) ** 2, (td + tn) ** 2)
                   for n in range(M)]
    fair = [binomial_polynomial(n, 1, 1) for n in range(M)]
    moments = []
    inactive_maxima = []
    H = td * td + tn * tn
    for ell in range(M):
        row = []
        denominator = scale * 2 ** (M - 1) * H ** ell
        for j in range(ell + 1):
            coeff = convolve(convolve(polynomials[j], polynomials[ell - j][::-1]),
                             fair[M - 1 - ell])
            assert sum(coeff) == 2 ** (M - 1) * H ** ell
            vals = [Q(sum(a * kernels[i + shift][side] for i, a in enumerate(coeff)),
                      denominator) for shift, side in ((1, 0), (1, 1), (0, 0), (0, 1))]
            row.append(vals)
        moments.append(row)
        inactive_maxima.append(C ** ell * max((x[1] + x[3]) / 2 for x in row))
    bounds = []
    by_h = []
    for h, target in enumerate(TARGETS):
        row = []
        for r in range(1, M - h + 1):
            if (h, r) in SIGNED_CASES:
                continue
            s = M - h - r
            ell = r - 1
            lam = (1 + eta ** (2 * r) + 2 * (1 if s % 2 == 0 else C) *
                   eta ** (r - h)) / (1 - eta ** 2)
            hr = ((1 + T) * (1 + T * v) ** ell +
                  (1 - T) * (1 - T * v) ** ell) / 2
            numerator = eta ** h * T * (1 - T * T) ** (M - 1) * C ** ell / hr ** 2
            coefficient = u - lam * v
            active = max(coefficient * (x[0] if coefficient < 0 else x[1]) +
                         (1 - u) * x[3] for x in moments[ell])
            residual = inactive_maxima[s] / 2 + eta ** h * C ** r * active / 2
            assert residual >= 0
            bound = numerator / (lam * numerator + residual)
            assert bound > target, (h, r, float(bound), float(target))
            bounds.append(bound)
            row.append((r, bound))
        least = min(row, key=lambda entry: entry[1])
        by_h.append({'h': h, 'cases': len(row), 'minimum_r': least[0],
                     'minimum_bound': float(least[1])})
    assert len(bounds) == 205
    return bounds, by_h


def cumulative_certificate():
    f = 1 - TARGETS[0]
    rho = (1 - TARGETS[1]) / f
    assert rho >= 1
    assert all(1 - p <= f * rho ** h for h, p in enumerate(TARGETS))
    assert CREDIT ** 5 >= rho ** 12
    for mask in range(63):
        chosen = [i for i in range(6) if mask & (1 << i)]
        n = len(chosen)
        pairs = sum(2 * j in chosen and 2 * j + 1 in chosen for j in range(3))
        penalty = n * (n - 1) // 2 + pairs
        assert rho ** penalty <= CREDIT ** n
        if n <= 4:
            for i in range(6):
                if i not in chosen:
                    q = f ** GROUPS[i // 2] * rho ** (n + int((i ^ 1) in chosen))
                    assert 0 < q < 1
    bonds = [1 - CREDIT * f ** g for g in GROUPS]
    assert all(0 < b < 1 for b in bonds)
    score = sum(bonds) - bonds[0] * bonds[1] * bonds[2]
    assert score > 1
    return {'proper_subsets': 63, 'bonds': [float(b) for b in bonds],
            'score': float(score), 'score_exact': str(score)}


def main():
    report = {'dimension': M // 2, 't': str(T), 'tangent': str(TANGENT),
              'targets': [str(p) for p in TARGETS], 'credit': str(CREDIT)}
    signed = []
    all_values = []
    for h, r in sorted(SIGNED_CASES):
        bounds = signed_residual(h, r)
        entry = {'h': h, 'r': r, 'orientations': len(bounds),
                 'residual_upper': float(max(bounds))}
        signed.append(entry)
        all_values.extend(bounds)
        print('PASS signed', entry, flush=True)
    report['signed_cases'] = signed
    bounds, by_h = old_residual_cases()
    all_values.extend(bounds)
    report['old_cases'] = by_h
    report['cumulative'] = cumulative_certificate()
    payload = '\n'.join(str(value) for value in all_values)
    report['bound_vector_sha256'] = hashlib.sha256(payload.encode()).hexdigest()
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

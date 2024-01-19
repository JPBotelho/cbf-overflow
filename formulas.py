from collections import defaultdict
from fractions import Fraction
from math import perm
from gmpy2 import mpz, bincoef


def exact_ovf_rate(n, m, k, h):
    j, nk = 2**h-1, n*k
    stirling = defaultdict(mpz)
    stirling[1,   1] = 1  # propagate stirling[0,0]
    stirling[1+j, 1] = -1

    for row in range(1, nk-j):
        comb_ = bincoef(row+j, j)
        for col in range(1, row+1):
            val = stirling[row, col]
            stirling[row+1,   col] += val * col
            stirling[row+1,   col+1] += val
            stirling[row+1+j, col+1] -= val * comb_
            del stirling[row, col]  # saving memory

    for row in range(nk-j, nk):
        for col in range(1, row+1):
            val = stirling[row, col]
            stirling[row+1,   col] += val * col
            stirling[row+1,   col+1] += val

    total = sum(stirling[nk, col] * perm(m, col)
                for col in range(1, nk+1))

    return float(1 - Fraction(total, m**nk))


def upper_bound(n, m, k, h):
    def p(val):
        return bincoef(n*k, val) * (1 - 1/m)**(n*k-val) * (1 / m)**val

    j = 2**h
    a = (n * k) / (m-1)
    val = a*(j+1) / (j * (j+1 - a))
    val *= p(j-1)
    val *= m
    return val

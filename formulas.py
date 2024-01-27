from collections import defaultdict
from fractions import Fraction
from math import perm, e, log
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

# Summary Cache: A Scalable Wide-Area
# Web Cache Sharing Protocol 
# https://pages.cs.wisc.edu/~jussara/papers/00ton.pdf 
# page 287
def worse_upper_bound(n, m, k, h):
    i = 2**h

    return m * ((e*log(2)/i)**i)

# Given a line of the binomical coefficient table
# compute the next lines until target_line
def cached_perm(line, target_line):   
    currentLine = line
    currentLineNumber = len(line) - 1

    while(currentLineNumber < target_line):
        newLine = [mpz(0) for i in range(len(currentLine)+1)]
        for col, val in enumerate(currentLine):
            newLine[col] += val 
            newLine[col+1] += val * (col+1)
        currentLineNumber += 1
        currentLine = newLine
    
    return newLine

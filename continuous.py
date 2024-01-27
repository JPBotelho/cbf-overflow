from collections import defaultdict
from fractions import Fraction
from math import log
from gmpy2 import mpz, bincoef
from time import time
from formulas import upper_bound, worse_upper_bound

# Notes:
# Naive upperbound is not always getting closer to the exact

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


def exact_ovf_rate(h):
    j = 2**h-1
    stirling = defaultdict(mpz)

    stirling[1,   1] = 1
    stirling[1+j, 1] = -1
    perm_row = [mpz(1)]
    start = time()
    row = 0
    prevExact = 0
    prevGood = 0
    prevBad = 0
    while(True):
        row += 1
        comb_ = bincoef(row+j, j)

        if(row % 10 == 0):
            lineTotal = 0
            nk = row
            k = 10
            n = int(row / k)
            m = int(round(row / log(2)))
            perm_row = cached_perm(perm_row, m)

        for col in range(1, row+1):
            val = stirling[row, col]
            stirling[row+1,   col] += val * col
            stirling[row+1,   col+1] += val
            stirling[row+1+j, col+1] -= val * comb_
            if(row % 10) == 0:
                lineTotal += (val * perm_row[col])
            del stirling[row, col]  # saving memory

        if(row % 10) == 0:
            exact =  float(1 - Fraction(lineTotal, int(m**nk)))
            approx_good = upper_bound(n, m, k, h)
            approx_bad = worse_upper_bound(n, m, k, h)

            print(f"-----------------")
            print(f"n={n}, m={m}, k={10}, h={h}")
            print(f"Exact: {exact:.5E}, var {(exact - prevExact):.5E}")

            if exact != 0:
                print(f"Tight Upperbound: {approx_good:.5E}, var {(approx_good - prevGood):.5E}, diff to exact: {(100*approx_good/exact - 100):.2}%")
                print(f"Naive Upperbound: {approx_bad:.5E}, var {(approx_bad - prevBad):.5E}, diff to tight: {round(100*approx_bad/approx_good - 100, 2)}%")
            else:
                print(f"Good Upperbound: {approx_good:.5E}, +_%")
                print(f"Bad Upperbound: {approx_bad:.5E}, +_%")

            prevExact = exact
            prevGood = approx_good
            prevBad = approx_bad
    print(f"took {time()-start}s")
exact = exact_ovf_rate(4)
            
from collections import defaultdict
from fractions import Fraction
from math import log
from gmpy2 import mpz, bincoef
from formulas import upper_bound, worse_upper_bound, cached_perm

def exact_ovf_rate_generator(h):
    j = 2**h-1

    stirling = defaultdict(mpz)
    stirling[1,   1] = 1
    stirling[1+j, 1] = -1

    # cached line of pascal's triangle
    # starts at line 1, gets expanded as needed
    perm_row = [mpz(1)]

    # Keep track of the previous values 
    prevExact = 0
    prevGood = 0
    prevBad = 0

    row = 0
    while(True):
        row += 1

        comb_ = bincoef(row+j, j)

        # for 1 in 1024 FP rate, k = 10 (https://hur.st/bloomfilter/)
        # we are interested in the (nk)th row of the stirling numbers
        # when row % 10, n = row / 10
        if(row % 10 == 0):
            lineTotal = 0
            nk = row
            k = 10
            n = int(row / k)
            m = int(round(row / log(2)))

            # generate pascal triangle line from cached line
            perm_row = cached_perm(perm_row, m)

        for col in range(1, row+1):
            val = stirling[row, col]
            stirling[row+1,   col] += val * col
            stirling[row+1,   col+1] += val
            stirling[row+1+j, col+1] -= val * comb_

            # same logic as above, when row % 10,
            # calculate row sum
            if(row % 10) == 0:
                lineTotal += (val * perm_row[col])

            del stirling[row, col]  # saving memory

        # print results
        if(row % 10) == 0:
            exact =  float(1 - Fraction(lineTotal, int(m**nk)))
            approx_good = upper_bound(n, m, k, h)
            approx_bad = worse_upper_bound(n, m, k, h)

            print(f"-----------------")
            print(f"n={n}, m={m}, k={10}, h={h}")
            print(f"Exact: {exact:.5E}, var {(exact - prevExact):.5E}")

            # in case exact = 0, the diff to exact calculations divide by 0 and break
            # the second case removes them
            if exact != 0:
                print(f"Tight Upperbound: {approx_good:.5E}, var {(approx_good - prevGood):.5E}, diff to exact: {(100*approx_good/exact - 100):.2}%")
                print(f"Naive Upperbound: {approx_bad:.5E}, var {(approx_bad - prevBad):.5E}, diff to exact: {round(100*approx_bad/exact - 100, 2)}%")
            else:
                print(f"Good Upperbound: {approx_good:.5E}, +_%")
                print(f"Bad Upperbound: {approx_bad:.5E}, +_%")

            prevExact = exact
            prevGood = approx_good
            prevBad = approx_bad

# arg: 4 is number of bits per counter
exact = exact_ovf_rate_generator(4)
            
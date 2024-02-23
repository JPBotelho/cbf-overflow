from collections import defaultdict

from fractions import Fraction

from math import log, e

from operator import delitem

from gmpy2 import mpz, bincoef

import time


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

    while (currentLineNumber < target_line):

        newLine = [mpz(0) for i in range(len(currentLine)+1)]

        for col, val in enumerate(currentLine):

            newLine[col] += val

            newLine[col+1] += val * (col+1)
        currentLineNumber += 1

        currentLine = newLine

    return newLine


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

    while (True):

        t0 = time.time()

        row += 1

        comb_ = bincoef(row+j, j)

        t1 = time.time()

        # for 1 in 1024 FP rate, k = 10 (https://hur.st/bloomfilter/)

        # we are interested in the (nk)th row of the stirling numbers

        # when row % 10, n = row / 10

        if (row % 10 == 0):

            lineTotal = 0

            nk = row

            k = 10

            n = int(row / k)

            m = int(round(row / log(2)))

            # generate pascal triangle line from cached line

            perm_row = cached_perm(perm_row, m)

        t2 = time.time()

        for col in range(1, row+1):

            val = stirling[row, col]

            stirling[row+1,   col] += val * col

            stirling[row+1,   col+1] += val

            stirling[row+1+j, col+1] -= val * comb_

            # same logic as above, when row % 10,

            # calculate row sum

            if (row % 10) == 0:

                lineTotal += (val * perm_row[col])

            del stirling[row, col]  # saving memory

        t3 = time.time()

        # print results

        if (row % 10) == 0:

            t4 = time.time()

            exact = float(1 - Fraction(int(lineTotal), int(m**nk)))

            t5 = time.time()

            approx_good = upper_bound(n, m, k, h)

            t6 = time.time()

            approx_bad = worse_upper_bound(n, m, k, h)

            t7 = time.time()

            # print(str(t0),)

            # print(f"-----------------")

            # print(f"n={n}, m={m}, k={10}, h={h}")

            # print(f"Exact: {exact:.5E}, var {(exact - prevExact):.5E}")

            # in case exact = 0, the diff to exact calculations divide by 0 and break

            # the second case removes them

            if exact != 0:

                print(t0, n, m, k, h, f"{exact:.5E}",
                      f"{(exact - prevExact):.5E}",
                      str(t5-t4), f"{approx_good:.5E}", f"{(approx_good - prevGood):.5E}",
                      f"{(100*approx_good/exact - 100):.5}%",
                      str(t6-t5), f"{approx_bad:.5E}",
                      f"{(approx_bad - prevBad):.5E}",
                      f"{round(100*approx_bad/exact -100, 2)}%", str(t7-t6),
                      sep=',')

                # print(f"Tight Upperbound: {approx_good:.5E}, var {(approx_good - prevGood):.5E}, diff to exact: {(100*approx_good/exact - 100):.2}%")

                # print(f"Naive Upperbound: {approx_bad:.5E}, var {(approx_bad - prevBad):.5E}, diff to exact: {round(100*approx_bad/exact - 100, 2)}%")

            else:

                print(t0, n, m, k, h, f"{exact:.5E}",
                      f"{(exact - prevExact):.5E}",
                      str(t5-t4), f"{approx_good:.5E}", f"{(approx_good - prevGood):.5E}",
                      "-", str(t6-t5), f"{approx_bad:.5E}",
                      f"{(approx_bad - prevBad):.5E}",
                      "-", str(t7-t6), sep=',')

                # print(f"Good Upperbound: {approx_good:.5E}, +_%")

                # print(f"Bad Upperbound: {approx_bad:.5E}, +_%")

            prevExact = exact

            prevGood = approx_good

            prevBad = approx_bad


print('t0', 'n', 'm', 'k',
      'h', 'eX', 'eVar', 'eTime',
      'tU', 'tVar', 'tDiff', 'tTime',
      'nU', 'nVar', 'nDiff', 'nTime',
      sep=',')

# arg: 4 is number of bits per counter

exact = exact_ovf_rate_generator(4)

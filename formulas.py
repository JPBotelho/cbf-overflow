from decimal import Decimal
from scipy.special import comb
from math import perm

def exact_ovf_rate(n, m, k, h, verbose=False):    
    j = 2**h-1
    stirling = {}
    stirling[(1, 1)] = 1
    stirling[(1+j, 1)] = -comb(j, j, exact=True)
    results = []    
    for row in range(1, n*k+1):
        for col in range(1, row+1):
            val = stirling[(row, col)]
            if row == n*k:
                results.append(((row, col), val)) 
            else:
                stirling[(row+1, col)] = stirling.get((row+1, col), 0) + col * val 
                stirling[(row+1, col+1)] = stirling.get((row+1, col+1), 0) + val
                stirling[(row+1+j, col+1)] = stirling.get((row+1+j, col+1), 0) - val * comb(row+j, j, exact=True)
            del stirling[(row, col)]
        if verbose and row % 100 == 0:
            print(f"Processed row {row}")

    total = 0
    for item in results:
        coords, val = item
        total += val * perm(m, coords[1])

    totalChoices = m**(n*k)

    probability = Decimal('1') - (Decimal(total) / Decimal(totalChoices)) 
    return probability

def upper_bound(n, m, k, h):
    j = 2**h 
    a = (n * k) / (m-1)

    def p(val):    
        return comb(n*k, val, exact=True) * (1 - 1/m)**(n*k-val) * (1 / m)**val

    val = a*(j+1) / (j * (j+1 - a))
    val *= p(j-1)
    val *= m
    return val

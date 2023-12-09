from decimal import Decimal as D, getcontext
from time import time
from formulas import exact_ovf_rate, upper_bound

cbfConstructions = [
    (100, 1443, 10, 4),
    (300, 4329, 10, 4),
    (500, 7214, 10, 4),
    (600, 8657, 10, 4),
    (700, 10099, 10, 4),
    (800, 11542, 10, 4),
    (900, 12985, 10, 4),
    (1000, 14426, 10, 4),
]

for cbf in cbfConstructions:
    start = time()
    getcontext().prec = 25000

    n, m, k, h = cbf
    exact = exact_ovf_rate(n, m, k, h, verbose=False)
    upper = upper_bound(n, m, k, h)
    

    print(f"-----------------")
    print(f"n={n}, m={m}, k={k}, h={h}")
    print(f"Exact: {exact:.5E}")
    print(f"Upper bound: {upper:.5E}")
    print(f"Deviation: {round(D('100')*D(upper) / exact - 100, 2)}%")
    print(f"Finished in {time()-start}s")



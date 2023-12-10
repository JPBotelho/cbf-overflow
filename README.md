# Exact overflow rate of a counting bloom filter

Some notes on the implementation:
The only tricky part is the restricted stirling numbers. 
For large CBFs doing naive recursion is too slow and they take up too much memory to be loaded at the same time. 
We construct them bottom-up in a set so that unnecessary numbers can be discarded.
Note that each number influences 3 other numbers and that S(0, 0, m) = 1. We iterate through each number in order, propagate its changes to the 3 numbers it affects, and delete it from memory.

Code uses factorials and combinations very naively, that bit can be speeded up a lot.

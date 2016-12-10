#!/usr/local/bin/python3

from scipy.special import binom
from sys import argv

def pr(N, A , k, t):
    numerator = binom(N-t * (k - 1),t)
    denominator = A ** (t * k)
    return binom(N-t * (k - 1),t) / A ** (t * k)

if len(argv) == 5:
    num_strings = 1
else:
    num_strings = argv[5]

print(pr(int(argv[1]),int(argv[2]),int(argv[3]),int(argv[4]))*int(num_strings))

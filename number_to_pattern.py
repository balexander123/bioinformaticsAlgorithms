#!/usr/local/bin/python3

# NumberToPattern(index, k)
#         if k = 1
#             return NumberToSymbol(index)
#         prefixIndex ← Quotient(index, 4)
#         r ← Remainder(index, 4)
#         symbol ← NumberToSymbol(r)
#         PrefixPattern ← NumberToPattern(prefixIndex, k − 1)
#         return concatenation of PrefixPattern with symbol

from sys import argv

number_to_symbol = {0:'A',1:'C',2:'G',3:'T'}

def quotient(index, divisor):
    return int(index / divisor)

def remainder(index, divisor):
    return index % divisor

def number_to_pattern(index, k):
    if k == 1:
        return [number_to_symbol[index]]
    prefix_index = quotient(index, 4)
    r = remainder(index, 4)
    symbol = number_to_symbol[r]
    return [symbol] + number_to_pattern(prefix_index, k - 1)

pattern = number_to_pattern(int(argv[1]),int(argv[2]))
pattern.reverse()
print(''.join(pattern))

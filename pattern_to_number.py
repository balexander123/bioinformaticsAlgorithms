#!/usr/local/bin/python3

# https://stepik.org/lesson/CS-Converting-Patterns-to-Numbers-and-Vice-Versa-3010/step/2?course=Stepic-Interactive-Text-for-Week-1&unit=8235
#
# PatternToNumber(Pattern)
#         if Pattern contains no symbols
#             return 0
#         symbol ← LastSymbol(Pattern)
#         Prefix ← Prefix(Pattern)
#         return 4 · PatternToNumber(Prefix) + SymbolToNumber(symbol)

from sys import argv

symbol_to_number = dict(A=0,C=1,G=2,T=3)

def pattern_prefix(pattern):
    if len(pattern) == 0:
        return ""
    return pattern[0:len(pattern)-1]

def pattern_to_number(pattern):
    if len(pattern) == 0:
        return 0
    symbol = pattern[len(pattern)-1]
    prefix = pattern_prefix(pattern)
    return 4 * pattern_to_number(prefix) + symbol_to_number[symbol]

print(pattern_to_number(argv[1]))
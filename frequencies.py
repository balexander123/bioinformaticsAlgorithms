#!/usr/local/bin/python3

from sys import argv

import array

import numpy

#
# utilities
#
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
#
# end utilities
#

# The pseudocode below generates a frequency array by first initializing every element in the frequency array to zero (4k operations) and then making a single pass down Text (approximately |Text| · k operations). For each k-mer Pattern that we encounter, we add 1 to the value of the frequency array corresponding to Pattern. As before, we refer to the k-mer beginning at position i of Text as Text(i, k).
#
#     ComputingFrequencies(Text, k)
#         for i ← 0 to 4k − 1
#             FrequencyArray(i) ← 0
#         for i ← 0 to |Text| − k
#             Pattern ← Text(i, k)
#             j ← PatternToNumber(Pattern)
#             FrequencyArray(j) ← FrequencyArray(j) + 1
#         return FrequencyArray



def computing_frequencies(text, k):
    frequency_array = numpy.zeros(4**k,numpy.int)
    for i in range(0,len(text) - k +1):
        pattern = text[i:i+k]
        j = pattern_to_number(pattern)
        frequency_array[j]+=1
    return ' '.join(map(str,frequency_array))

print(computing_frequencies(argv[1], int(argv[2])))
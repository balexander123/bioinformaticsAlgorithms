#!/usr/local/bin/python3

from sys import argv

import array

import numpy

import argparse

parser = argparse.ArgumentParser(description='Find frequent k-mer words with maximum allowed mismatches d.')
parser.add_argument('--file', help='the file that contains a sequence of nucleotide symbols')
parser.add_argument('--genome', help='genome string of nucleotide symbols')
parser.add_argument('-k', help='k-mer length')
parser.add_argument('-d', help='maximum hamming distance')
parser.add_argument('--reverse', help='include reverse complemnet')

#
# utilities
#
symbol_to_number = dict(A=0, C=1, G=2, T=3)


def pattern_prefix(pattern):
    if len(pattern) == 0:
        return ""
    return pattern[0:len(pattern) - 1]


def pattern_to_number(pattern):
    if len(pattern) == 0:
        return 0
    symbol = pattern[len(pattern) - 1]
    prefix = pattern_prefix(pattern)
    return 4 * pattern_to_number(prefix) + symbol_to_number[symbol]


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
    frequency_array = numpy.zeros(4 ** k, numpy.int)
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i + k]
        j = pattern_to_number(pattern)
        frequency_array[j] += 1
    return ' '.join(map(str, frequency_array))

def suffix(pattern):
    if len(pattern) == 0:
        return ""
    return pattern[1:len(pattern)]

def hamming_distance(p,q):
    if len(p) != len(q):
        return -1
    hamming_count = 0
    for i, val in enumerate(p):
        if p[i] != q[i]:
            hamming_count += 1
    return hamming_count

def first_symbol(pattern):
    if len(pattern) == 0:
        return ""
    return pattern[0:1]

def print_flat(s):
    for item in s:
        print(item)

def neighbors(pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighborhood = set()
    suffix_neighbors = neighbors(suffix(pattern), d)
    for text in suffix_neighbors:
        if hamming_distance(suffix(pattern), text) < d:
            for x in {"A", "C", "G", "T"}:
                neighborhood.add(x + text)
        else:
            neighborhood.add(first_symbol(pattern) + text)
    return neighborhood

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

def approximate_pattern_count(text, pattern, d):
    count = 0
    indexes = []
    for i in range(0, len(text) - len(pattern) + 1):
        pattern_prime = text[i:i+len(pattern)]
        if hamming_distance(pattern,pattern_prime) <= d:
            count += 1
            indexes.append(i)
    return len(indexes)

base_complement = {'A':'T','C':'G','T':'A','G':'C'}

def reverse_complement(pattern):
    reverse_comp = ""
    for c in pattern[::-1]:
        reverse_comp+=base_complement[c]
    return reverse_comp

#
# end utilities
#


# FrequentWordsWithMismatches(Text, k, d)
#     FrequentPatterns ← an empty set
#     for i ← 0 to 4k − 1
#         Close(i) ← 0
#         FrequencyArray ← 0
#     for i ← 0 to |Text| − k
#         Neighborhood ← Neighbors(Text(i, k), d)
#         for each Pattern from Neighborhood
#             index ← PatternToNumber(Pattern)
#             Close(index) ← 1
#     for i ← 0 to 4k − 1
#         if Close(i) = 1
#             Pattern ← NumberToPattern(i, k)
#             FrequencyArray(i) ← ApproximatePatternCount(Text, Pattern, d)
#     maxCount ← maximal value in FrequencyArray
#     for i ← 0 to 4k − 1
#         if FrequencyArray(i) = maxCount
#             Pattern ← NumberToPattern(i, k)
#             add Pattern to the set FrequentPatterns
#    return FrequentPatterns

def frequent_words_with_mismatches(text, k, d, reverse):
    frequent_patterns = set()
    close = numpy.zeros(4 ** k, numpy.int)
    frequency_array = numpy.zeros(4**k,numpy.int)

    for i in range(0, len(text) - k):
        neighborhood = neighbors(text[i:i+k], d)
        for pattern in neighborhood:
            index = pattern_to_number(''.join(map(str, pattern)))
            close[index] = 1
    for i in range(0, (4 ** k)-1):
        if close[i] == 1:
            pattern = ''.join(map(str, number_to_pattern(i,k)))
            frequency_array[i] = approximate_pattern_count(text, pattern, d)
    maxCount = max(frequency_array)
    for i in range(0, (4**k)-1):
        if frequency_array[i] == maxCount:
            pattern = ''.join(map(str, number_to_pattern(i, k)))
            frequent_patterns.add(pattern)
    return frequent_patterns

args = parser.parse_args()

if args.file != None:
    with open(args.file, 'r') as myfile: data=myfile.read().replace('\n', '')
else:
    data = args.genome

print(' '.join(map(str, frequent_words_with_mismatches(data,int(args.k),int(args.d),args.reverse))))
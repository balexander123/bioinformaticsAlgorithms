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

args = parser.parse_args()

# ClumpFinding(Genome, k, t, L)
#     FrequentPatterns ← an empty set
#     for i ← 0 to 4k − 1
#         Clump(i) ← 0
#     for i ← 0 to |Genome| − L
#         Text ← the string of length L starting at position i in Genome 
#         FrequencyArray ← ComputingFrequencies(Text, k)
#         for index ← 0 to 4k − 1
#             if FrequencyArray(index) ≥ t
#                 Clump(index) ← 1
#     for i ← 0 to 4k − 1
#         if Clump(i) = 1
#             Pattern ← NumberToPattern(i, k)
#             add Pattern to the set FrequentPatterns
#     return FrequentPatterns

# BetterClumpFinding(Genome, k, t, L)
#     FrequentPatterns ← an empty set
#     for i ← 0 to 4k − 1
#         Clump(i) ← 0
#     Text ← Genome(0, L)
#     FrequencyArray ← ComputingFrequencies(Text, k)
#     for i ← 0 to 4k − 1
#         if FrequencyArray(i) ≥ t
#             Clump(i) ← 1
#     for i ← 1 to |Genome| − L
#         FirstPattern ← Genome(i − 1, k)
#         index ← PatternToNumber(FirstPattern)
#         FrequencyArray(index) ← FrequencyArray(index) − 1
#         LastPattern ← Genome(i + L − k, k)
#         index ← PatternToNumber(LastPattern)
#         FrequencyArray(index) ← FrequencyArray(index) + 1
#         if FrequencyArray(index) ≥ t
#             Clump(index) ← 1
#     for i ← 0 to 4k − 1
#         if Clump(i) = 1
#             Pattern ← NumberToPattern(i, k)
#             add Pattern to the set FrequentPatterns
#     return FrequentPatterns
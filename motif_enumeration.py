#!/usr/local/bin/python3

from sys import argv

#
# utilities
#
import array

import numpy

# Neighbors(Pattern, d)
#   if d = 0
#       return {Pattern}
#   if | Pattern | = 1
#       return {A, C, G, T}
#   Neighborhood ← an empty set
#   SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
#   for each string Text from SuffixNeighbors
#       if HammingDistance(Suffix(Pattern), Text) < d
#           for each nucleotide x
#               add x • Text to Neighborhood
#       else
#           add FirstSymbol(Pattern) • Text to Neighborhood
#   return Neighborhood

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

def pattern_count(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

def frequent_words(text, k):
    frequent_patterns = []
    count=array.array('i',(0,)*len(text))
    for i in range(len(text) - k):
        pattern = text[i:i+k]
        count[i] = pattern_count(text,pattern)
    maxCount = max(count)
    for i in range(len(text) - k):
        if count[i] == maxCount:
            frequent_patterns.append(text[i:i+k])
    return set(frequent_patterns)

#
# end utilities
#

# MotifEnumeration(Dna, k, d)
#     Patterns ← an empty set
#     for each k-mer Pattern in Dna
#         for each k-mer Pattern’ differing from Pattern by at most d
#           mismatches
#             if Pattern' appears in each string from Dna with at most d
#             mismatches
#                 add Pattern' to Patterns
#     remove duplicates from Patterns
#     return Patterns
        
def motif_enumeration(dna, k, d):
    patterns= set()
    for text in dna:
        for k_mer_pattern in frequent_words(text, d):
            for pattern_prime in neighbors(k_mer_pattern, d):
                print(pattern_prime)

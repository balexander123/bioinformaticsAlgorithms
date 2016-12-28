#!/usr/local/bin/python3

from sys import argv

#
# utilities
#
import array

import numpy

import argparse

parser = argparse.ArgumentParser(description='Find frequent k-mer words with maximum allowed mismatches d.')
parser.add_argument('--file', nargs='?', help='the file that contains dna motifs')
parser.add_argument('-k', nargs='?', help='k-mer size')
parser.add_argument('--distance', nargs='?', help='allowable hamming distance')

args = parser.parse_args()

# Neighbors(Pattern, d)
#   if d = 0
#       return {Pattern}
#   if | Pattern | = 1
#       return {A, C, G, T}
#   Neighborhood <- an empty set
#   SuffixNeighbors <- Neighbors(Suffix(Pattern), d)
#   for each string Text from SuffixNeighbors
#       if HammingDistance(Suffix(Pattern), Text) < d
#           for each nucleotide x
#               add Text to Neighborhood
#       else
#           add FirstSymbol(Pattern) . Text to Neighborhood
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

#
# end utilities
#

# MotifEnumeration(Dna, k, d)
#     Patterns <- an empty set
#     for each k-mer Pattern in Dna
#         for each k-mer Pattern' differing from Pattern by at most d mismatches
#             if Pattern' appears in each string from Dna with at most d
#             mismatches
#                 add Pattern' to Patterns
#     remove duplicates from Patterns
#     return Patterns

def k_mer_patterns_from_text(text, k):
    k_mer_patterns = []
    for i in range(0, len(text) - k + 1):
        k_mer_patterns.append(text[i:i + k])
    return set(k_mer_patterns)

def k_mer_patterns_from_dna(dna,k):
    k_mer_patterns = []
    for text in dna:
        k_mer_patterns.append(k_mer_patterns_from_text(text,k))
    return k_mer_patterns

def appears_in(dna, pattern, d):
    k = len(pattern)
    num_strings = len(dna)
    num_matches = 0
    for q in k_mer_patterns_from_dna(dna, k):
        if hamming_distance(pattern, q) <= d:
            return True
    return False
        
def motif_enumeration(dna, k, d):
    patterns= set()
    k_mer_patterns = k_mer_patterns_from_dna(dna,k)
    for k_mer_pattern in k_mer_patterns:
        pattern_neighbors = neighbors(k_mer_pattern, d)
        for pattern_prime in pattern_neighbors:
            if appears_in(dna, pattern_prime, d):
                patterns.add(pattern_prime)
    return patterns

def dna_from_file(dna_file):
    with open(dna_file, 'r') as myfile:
        dna = myfile.readlines()
        for i, text in enumerate(dna):
            dna[i] = text.replace('\n', '')
    return dna

print(motif_enumeration(dna_from_file(args.file),int(args.k),int(args.distance)))
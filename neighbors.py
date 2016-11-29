#!/usr/local/bin/python3

from sys import argv

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

print_flat(neighbors(argv[1], int(argv[2])))
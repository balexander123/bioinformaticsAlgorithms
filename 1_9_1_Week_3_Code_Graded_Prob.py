#!/usr/local/bin/python3

import sys

from sys import argv


def suffix(pattern):
    if len(pattern) == 0:
        return ""
    return pattern[1:len(pattern)]


def hamming_distance(p, q):
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


def neighbors(pattern, d):
    if d == 0:
        return {pattern}
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

def k_mer_patterns_from_text(text, k):
    k_mer_patterns = []
    for i in range(0, len(text) - k + 1):
        k_mer_patterns.append(text[i:i + k])
    return k_mer_patterns


def k_mer_patterns_from_dna(dna, k):
    k_mer_patterns = []
    for text in dna:
        for k_mer in k_mer_patterns_from_text(text, k):
            k_mer_patterns.append(k_mer)
    return k_mer_patterns


def appears_in(dna, pattern, d):
    k = len(pattern)
    num_strings = len(dna)
    for i in range(0, num_strings):
        num_matches = 0
        for q in k_mer_patterns_from_text(dna[i], k):
            if hamming_distance(pattern, q) <= d:
                num_matches += 1
        if num_matches == 0:
            return False
    return True


def motif_enumeration(k, d, dna):
    patterns = set()
    k_mer_patterns = k_mer_patterns_from_dna(dna, k)
    for k_mer_pattern in k_mer_patterns:
        pattern_neighbors = neighbors(k_mer_pattern, d)
        for pattern_prime in pattern_neighbors:
            if appears_in(dna, pattern_prime, d):
                patterns.add(pattern_prime)
    return patterns

if len(argv) == 2:
    with open(argv[1], 'r') as input_file:
        first_line = input_file.readline()
        k = int(first_line.split()[0])
        d = int(first_line.split()[1])
        dna = []
        for dna_string in input_file.readlines():
            dna.append(dna_string)
else:
    first_line = sys.stdin.readline()
    k = int(first_line.split()[0])
    d = int(first_line.split()[1])
    dna = []
    for dna_string in sys.stdin.readlines():
        dna.append(dna_string)

    for i, text in enumerate(dna):
        dna[i] = text.replace('\n', '')

print(' '.join(sorted(motif_enumeration(k, d, dna))))
#!/usr/local/bin/python3

import sys

from sys import argv

def hamming_distance(p,q):
    if len(p) != len(q):
        return -1
    hamming_count = 0
    for i, val in enumerate(p):
        if p[i] != q[i]:
            hamming_count += 1
    return hamming_count

def k_mer_patterns_from(text,k):
    k_mer_patterns = []
    for i in range(0,len(text) - k + 1):
        k_mer_patterns.append(text[i:i+k])
    return set(k_mer_patterns)

def distance_between_pattern_and_strings(pattern, dna):
    k = len(pattern)
    distance = 0
    for text in dna:
        hamming_dist = 999
        k_mer_patterns = k_mer_patterns_from(text,k)
        for pattern_prime in k_mer_patterns:
            dist = hamming_distance(pattern, pattern_prime)
            if hamming_dist > dist:
                hamming_dist = dist
        distance += hamming_dist
    return distance




if len(argv) == 2:
    with open(argv[1], 'r') as input_file:
        pattern = input_file.readline()
        dna = []
        dna = input_file.readline()
        dna = dna.split()
else:
    pattern = sys.stdin.readline()
    dna = []
    dna = sys.stdin.readline()
    dna = dna.split()

for i, text in enumerate(dna):
    dna[i] = text.replace('\n', '')


print(distance_between_pattern_and_strings(pattern[0:len(pattern)-1],dna))


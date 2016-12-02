#!/usr/local/bin/python3

from sys import argv

def hamming_distance(p,q):
    if len(p) != len(q):
        return -1
    hamming_count = 0
    for i, val in enumerate(p):
        if p[i] != q[i]:
            hamming_count += 1
    return hamming_count

def approximate_pattern_matches(pattern, text, d):
    count = 0
    indexes = []
    for i in range(0, len(text) - len(pattern) + 1):
        pattern_prime = text[i:i+len(pattern)]
        if hamming_distance(pattern,pattern_prime) <= d:
            count += 1
            indexes.append(i)
    return indexes

approximate_match_indexes = approximate_pattern_matches(argv[1],argv[2],int(argv[3]))
print(len(approximate_match_indexes))
print(*approximate_match_indexes, sep=' ')


#!/usr/local/bin/python3

from sys import argv

import sys

import array

import numpy

import argparse

parser = argparse.ArgumentParser(description='Find frequent k-mer words with maximum allowed mismatches d.')
parser.add_argument('--file', nargs='?', help='the file that contains dna motifs')
parser.add_argument('-k', nargs='?', help='k-mer size')
parser.add_argument('--pattern', nargs='?', help='pattern used to determine hamming distance of dna motifs')

args = parser.parse_args()

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
        for pattern_prime in k_mer_patterns_from(text,k):
            dist = hamming_distance(pattern, pattern_prime)
            if hamming_dist > dist:
                hamming_dist = dist
        distance += hamming_dist
    return distance

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

def median_string(dna, k):
    distance = 999
    median = ''
    for i in range(0, (4 ** k) - 1):
        pattern = ''.join(number_to_pattern(i,k))
        k_mer_dist = distance_between_pattern_and_strings(pattern, dna)
        if distance > k_mer_dist:
            distance = k_mer_dist
            median = pattern
    return median

def dna_from_file():
    with open(args.file, 'r') as myfile:
        dna = myfile.readlines()
        for i, text in enumerate(dna):
            dna[i] = text.replace('\n', '')
    return dna


print(median_string(dna_from_file(),int(args.k)))



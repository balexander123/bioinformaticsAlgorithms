#!/usr/local/bin/python3

from sys import argv

#
# utilities
#
import array

import numpy

# import argparse
#
# parser = argparse.ArgumentParser(description='Find frequent k-mer words with maximum allowed mismatches d.')
# parser.add_argument('--file', nargs='?', help='the file that contains dna motifs')
# parser.add_argument('-k', nargs='?', help='k-mer size')
# parser.add_argument('--distance', nargs='?', help='allowable hamming distance')
#
# args = parser.parse_args()

class MotifEnumeration:

    def __init__(self, file):
        self.initialize_from_file(file)

    def suffix(self, pattern):
        if len(pattern) == 0:
            return ""
        return pattern[1:len(pattern)]

    def hamming_distance(self, p,q):
        if len(p) != len(q):
            return -1
        hamming_count = 0
        for i, val in enumerate(p):
            if p[i] != q[i]:
                hamming_count += 1
        return hamming_count

    def first_symbol(self, pattern):
        if len(pattern) == 0:
            return ""
        return pattern[0:1]

    def neighbors(self, pattern, d):
        if d == 0:
            return { pattern }
        if len(pattern) == 1:
            return {"A", "C", "G", "T"}
        neighborhood = set()
        suffix_neighbors = self.neighbors(self.suffix(pattern), d)
        for text in suffix_neighbors:
            if self.hamming_distance(self.suffix(pattern), text) < d:
                for x in {"A", "C", "G", "T"}:
                    neighborhood.add(x + text)
            else:
                neighborhood.add(self.first_symbol(pattern) + text)
        return neighborhood
#
# end utilities
#

    def k_mer_patterns_from_text(self, text, k):
        k_mer_patterns = []
        for i in range(0, len(text) - k + 1):
            k_mer_patterns.append(text[i:i + k])
        return k_mer_patterns

    def k_mer_patterns_from_dna(self, dna, k):
        k_mer_patterns = []
        for text in dna:
            for k_mer in self.k_mer_patterns_from_text(text,k):
                k_mer_patterns.append(k_mer)
        return k_mer_patterns

    def appears_in(self, dna, pattern, d):
        k = len(pattern)
        num_strings = len(dna)
        for i in range(0,num_strings):
            num_matches = 0
            for q in self.k_mer_patterns_from_text(dna[i], k):
                if self.hamming_distance(pattern, q) <= d:
                    num_matches += 1
            if num_matches == 0:
                return False
        return True
        
    def motif_enumeration(self):
        patterns= set()
        k_mer_patterns = self.k_mer_patterns_from_dna(self.dna,self.k)
        for k_mer_pattern in k_mer_patterns:
            pattern_neighbors = self.neighbors(k_mer_pattern, self.distance)
            for pattern_prime in pattern_neighbors:
                if self.appears_in(self.dna, pattern_prime, self.distance):
                    patterns.add(pattern_prime)
        return patterns

    def initialize_from_file(self, dna_file):
        with open(dna_file, 'r') as myfile:
            first_line = myfile.readline()
            self.k = int(first_line.split(' ')[0])
            self.distance = int(first_line.split(' ')[1].replace('\n', ''))
            self.dna = myfile.readlines()
            for i, text in enumerate(self.dna):
                self.dna[i] = text.replace('\n', '')

if __name__ == "__main__":
    motif_enum = MotifEnumeration(argv[1])
    print(' '.join(sorted(motif_enum.motif_enumeration())))
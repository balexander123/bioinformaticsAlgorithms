#!/usr/local/bin/python3

import random

import sys

import profile

from sys import argv

def k_mer_patterns_from(text,k):
    k_mer_patterns = []
    for i in range(0,len(text) - k + 1):
        k_mer_patterns.append(text[i:i+k])
    return k_mer_patterns

def first_kmers_from_dna(dna_strings,k):
    first_kmers = []
    for dna in dna_strings:
        first_kmers.append(dna[0:k])
    return first_kmers

def _random_k_mers(dna, k):
    random_k_kers = []
    for dna in dna:
        # generate random index range 0:len(dna)-self.k
        start_index = random.randrange(0, len(dna) - k)
        random_k_kers.append(dna[start_index:start_index + k])
    return random_k_kers

def most_probable(profile,dna_motifs):
    most_probables = []
    for dna in dna_motifs:
        m_probable = 0
        m_probable_k_mer = None
        # for each k-mer in dna
        for column in range(0,len(dna) - profile.k + 1):
            probability = 1
            k_mer = dna[column:column+profile.k]
            for j in range(len(k_mer)):
                # compute probability of k-mer using profile
                row = profile.nucleotideIndex[k_mer[j]]
                prob = profile.profile_matrix[profile.nucleotideIndex[k_mer[j]]][j]
                probability *= profile.profile_matrix[profile.nucleotideIndex[k_mer[j]]][j]
            if probability > m_probable:
                m_probable = probability
                m_probable_k_mer = k_mer
        most_probables.append(m_probable_k_mer)
    return most_probables

def profile_from_motifs(dna, k):
    return profile.Profile(None,k,dna,0)

# GreedyMotifSearch(Dna, k, t)
#     BestMotifs ← motif matrix formed by first k-mers in each string from Dna
#     for each k-mer Motif in the first string from Dna
#         Motif1 ← Motif
#         for i = 2 to t
#             form Profile from motifs Motif1, …, Motifi - 1
#             Motifi ← Profile-most probable k-mer in the i-th string in Dna
#         Motifs ← (Motif1, …, Motift)
#         if Score(Motifs) < Score(BestMotifs)
#             BestMotifs ← Motifs
#     return BestMotifs

def hamming_distance(p,q):
    if len(p) != len(q):
        return -1
    hamming_count = 0
    for i, val in enumerate(p):
        if p[i] != q[i]:
            hamming_count += 1
    return hamming_count

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

def median_strings(dna, k):
    medians = []
    for text in dna:
        medians.append(median_string([text],k))
    return medians

def first_kmers(dna, k):
    kmers = []
    for text in dna:
        kmers.append(text[0:k])
    return kmers

def greedy_motif_search(dna, k, t):
    first_motifs = first_kmers(dna, k)
    best_profile = profile.Profile(None, k, first_motifs, 1)
    best_score = [best_profile.motif_score(), first_motifs]
    for i in range(0,len(dna[0])-k+1):
        motifs = [dna[0][i:i+k]]
        motif_profile = profile.Profile(None, k, motifs, 1)
        for j in range(1,len(dna)):
            most_probable_kmers = most_probable(motif_profile, [dna[j]])
            for motif in most_probable_kmers:
                if motif is not None:
                    motifs.extend([motif])
                else:
                    motifs.extend([first_motifs[j]])
            motif_profile = profile.Profile(None, k, motifs, 1)
            resulting_kmers = [dna[j], most_probable_kmers, motif_profile.motif_score()]
        motif_profile = profile.Profile(None, k, motifs, 1)
        current_score = motif_profile.motif_score()
        if current_score < best_score[0]:
            best_score = [current_score, motifs]
    return best_score[1]

if len(argv) == 2:
    with open(argv[1], 'r') as input_file:
        first_line = input_file.readline()
        k = int(first_line.split()[0])
        t = int(first_line.split()[1])
        dna = []
        for dna_string in input_file.readlines():
            dna.append(dna_string)
else:
    first_line = sys.stdin.readline()
    k = int(first_line.split()[0])
    t = int(first_line.split()[1])
    dna = []
    for dna_string in sys.stdin.readlines():
        dna.append(dna_string)

for i, text in enumerate(dna):
    dna[i] = text.replace('\n', '')

print('\n'.join(greedy_motif_search(dna,k,t)))
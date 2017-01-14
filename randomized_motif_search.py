#!/usr/local/bin/python3

import random

import sys

import profile

import timeit

from sys import argv

class RandomizedMotifSearch:

    def _random_k_mers(self, k, dna_strings):
        random_k_kers = []
        for dna in dna_strings:
            # generate random index range 0:len(dna)-k
            start_index = random.randrange(0, len(dna) - k + 1)
            random_k_kers.append(dna[start_index:start_index + k])
        return random_k_kers

    def most_probable(self, profile, dna_motifs):
        most_probables = []
        for dna in dna_motifs:
            m_probable = 0
            m_probable_k_mer = None
            # for each k-mer in dna
            for column in range(0, len(dna) - profile.k + 1):
                probability = 1
                k_mer = dna[column:column + profile.k]
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

    def randomized_motif_search(self, k, dna):
        motifs = self._random_k_mers(k, dna)
        best_motifs = motifs
        while True:
            motif_profile = profile.Profile(None, k, motifs)
            motifs = self.most_probable(motif_profile, dna)
            motif_score = profile.Profile(None, k, motifs).motif_score()
            best_score = profile.Profile(None, k, best_motifs).motif_score()
            if motif_score < best_score:
                best_motifs = motifs
            else:
                return [best_score, best_motifs]

    def randomized_motif_searchX(self, k, dna, times):
        best_score = len(dna)*k
        while times > 0:
            motif_score = self.randomized_motif_search(k, dna)
            if motif_score[0] < best_score:
                best_motifs = motif_score[1]
                best_score = motif_score[0]
            times -= 1
        return best_motifs

if __name__ == '__main__':
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

    motifs = []
    randomized_motif_searcher = RandomizedMotifSearch()
    motifs = randomized_motif_searcher.randomized_motif_searchX(k,dna,1000)
    print(' '.join(motifs))
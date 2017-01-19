#!/usr/local/bin/python3

import random

import sys

import timeit

from sys import argv

from sys import argv

from enum import Enum

import numpy

import math

class Profile:


    class Nucleotides(Enum):
        A = 'A'
        C = 'C'
        G = 'G'
        T = 'T'

    nucleotideIndex = {'A': 0, 'C': 1, 'G': 2, 'T': 3 }

    def __init__(self, file=None, k=0, dna=None, psuedocount=1):
        if dna is None:
            dna = []
        if file != None:
            self.initialize_from_file(file)
        else:
            self.initialize_with(k,dna)
        self.initialize_profile(psuedocount)

    def initialize_profile(self, psuedocount):
        self.pseudocount = psuedocount
        self.score = 0
        if len(self.dna) != 0:
            self.motif_matrix = Matrix = [['' for x in range(len(self.dna[0]))] for y in range(len(self.dna))]
            for i in range(0, len(self.dna)):
                self.motif_matrix[i] = list(self.dna[i])
            self.profile_matrix = Matrix = [[0 for x in range(len(self.dna[0]))] for y in range(4)]
            self.count_matrix = Matrix = [[0 for x in range(len(self.dna[0]))] for y in range(4)]
            self.score_array = numpy.zeros(len(self.motif_matrix[0]), numpy.int)
            self.motif_count()
            self.apply_laplace_rule_of_succession(self.pseudocount)
            self.compute_profile_matrix()
            self.motif_score()

    def motif_score(self):
        self.score = 0
        if len(self.dna) > 0:
            for column in range(len(self.motif_matrix[0])):
                motif_column = [self.motif_matrix[i][column] for i in range(0,self.num_motifs)]
                self.score_array[column] = len(motif_column) - max(
                    motif_column.count(self.Nucleotides.A.value),
                    motif_column.count(self.Nucleotides.C.value),
                    motif_column.count(self.Nucleotides.G.value),
                    motif_column.count(self.Nucleotides.T.value))
            self.score = sum(self.score_array)
        return(self.score)

    def motif_count(self):
        for column in range(len(self.motif_matrix[0])):
            motif_column = [self.motif_matrix[i][column] for i in range(0,self.num_motifs)]
            self.count_matrix[0][column] = motif_column.count(self.Nucleotides.A.value)
            self.count_matrix[1][column] = motif_column.count(self.Nucleotides.C.value)
            self.count_matrix[2][column] = motif_column.count(self.Nucleotides.G.value)
            self.count_matrix[3][column] = motif_column.count(self.Nucleotides.T.value)

    def apply_laplace_rule_of_succession(self,pseudocount):
        self.pseudocount = pseudocount
        if pseudocount != 0:
            for row in range(len(self.count_matrix)):
                for col in range(len(self.count_matrix[0])):
                    self.count_matrix[row][col] += pseudocount

    def compute_profile_matrix(self):
        for row in range(len(self.count_matrix)):
            for col in range(len(self.count_matrix[0])):
                count_column = [self.count_matrix[i][col] for i in range(0, 4)]
                self.profile_matrix[row][col] = self.count_matrix[row][col] / sum(count_column)

    def compute_profile_matrix_entropy(self):
        total_entropy = 0
        for column in range(len(self.profile_matrix[0])):
            profile_column = [self.profile_matrix[i][column] for i in range(0,4)]
            entropy = 0
            for probability in profile_column:
                entropy += probability * math.log(probability,2)
            entropy *= -1
            total_entropy += entropy
        return total_entropy

    def initialize_with(self, k, dna):
        self.k = k
        self.num_motifs = len(dna)
        self.dna = dna

    def initialize_from_motifs_file(self, file):
        self.k = 0
        with open(file, 'r') as motif_file:
            self.dna = motif_file.readlines()
            self.num_motifs = len(self.dna)
            for i, text in enumerate(self.dna):
                self.dna[i] = text.replace('\n', '')

    def initialize_from_file(self, file):
        with open(file, 'r') as myfile:
            first_line = myfile.readline()
            self.k = int(first_line.split(' ')[0])
            self.num_motifs = int(first_line.split(' ')[1].replace('\n', ''))
            self.dna = myfile.readlines()
            for i, text in enumerate(self.dna):
                self.dna[i] = text.replace('\n', '')

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
            motif_profile = Profile(None, k, motifs)
            motifs = self.most_probable(motif_profile, dna)
            motif_score = Profile(None, k, motifs).motif_score()
            best_score = Profile(None, k, best_motifs).motif_score()
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
    print('\n'.join(motifs))
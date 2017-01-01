#!/usr/local/bin/python3

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

    def __init__(self, file='', k=0, dna=None, N=1):
        if dna is None:
            dna = []
        if file != '':
            self.initialize_from_file(file)
        else:
            self.initialize_with(k,dna)
        self.initialize_profile()

    # def __init__(self, dna):
    #     self.initialize_profile()

    def initialize_profile(self):
        self.motif_matrix = Matrix = [['' for x in range(len(self.dna[0]))] for y in range(len(self.dna))]
        for i in range(0, len(self.dna)):
            self.motif_matrix[i] = list(self.dna[i])
        self.profile_matrix = Matrix = [[0 for x in range(len(self.dna[0]))] for y in range(4)]
        self.count_matrix = Matrix = [[0 for x in range(len(self.dna[0]))] for y in range(4)]
        self.score_array = numpy.zeros(len(self.motif_matrix[0]), numpy.int)
        self.pseudocount = 0
        self.score = 0

        self.motif_count()
        self.apply_laplace_rule_of_succession(1)
        self.compute_profile_matrix()

        self.data = []

    def motif_score(self):
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
        self.num_motifs = len(dna[0])
        self.dna = dna

    def initialize_from_file(self, file):
        with open(file, 'r') as myfile:
            first_line = myfile.readline()
            self.k = int(first_line.split(' ')[0])
            self.num_motifs = int(first_line.split(' ')[1].replace('\n', ''))
            self.dna = myfile.readlines()
            for i, text in enumerate(self.dna):
                self.dna[i] = text.replace('\n', '')


if __name__ == "__main__":
    profile = Profile(argv[1])
    print(profile.motif_score())
    print(profile.compute_profile_matrix_entropy())
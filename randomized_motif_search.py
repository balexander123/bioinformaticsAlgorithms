#!/usr/local/bin/python3

import random

import profile

class RandomizedMotifSearch:

    def _random_k_mers(self):
        random_k_kers = []
        for dna in self.dna:
            # generate random index range 0:len(dna)-self.k
            start_index = random.randrange(0, len(dna) - self.k)
            random_k_kers.append(dna[start_index:start_index + self.k])
        return random_k_kers

    def most_probable(self,profile,dna_motifs):
        most_probables = []
        for dna in dna_motifs:
            m_probable = 0
            # for each k-mer in dna
            for column in range(0,len(dna) - profile.k + 1):
                probability = 1
                k_mer = dna[column:column+profile.k]
                for j in range(len(k_mer)):
                    # compute probability of k-mer using profile
                    prob = profile.profile_matrix[profile.nucleotideIndex[k_mer[j]]][j]
                    probability *= profile.profile_matrix[profile.nucleotideIndex[k_mer[j]]][j]
                if probability > m_probable:
                    m_probable = probability
                    m_probable_k_mer = k_mer
            most_probables.append(m_probable_k_mer)
        return most_probables
#!/usr/local/bin/python3

import random
import random_utils
import profile

def _random_k_mers(k, dna_strings):
    random_k_kers = []
    for dna in dna_strings:
        # generate random index range 0:len(dna)-k
        start_index = random.randrange(0, len(dna) - k + 1)
        random_k_kers.append(dna[start_index:start_index + k])
    return random_k_kers

def random_generated_kmer(profile,dna_string):
    distribution = []
    kmers = []
    i=0
    for column in range(0, len(dna_string) - profile.k + 1):
        probability = 1
        k_mer = dna_string[column:column + profile.k]
        kmers.append(k_mer)
        for j in range(len(k_mer)):
            # compute probability of k-mer using profile
            probability *= profile.profile_matrix[profile.nucleotideIndex[k_mer[j]]][j]
        distribution[i] = probability
        i += 1
    return kmers[random_utils.roll(distribution)]


def gibbs_sampler(dna, k, t, n):
    best_motifs = motifs = _random_k_mers(k, dna)
    for j in range(1,n):
        i = random.randrange(0, t-1)
        motif = motifs.remove(i)
        motifs.insert(i,random_generated_kmer(profile.Profile(None, k, motifs),motif))
        if profile.Profile(None, k, motifs).motif_score() < profile.Profile(None, k, best_motifs).motif_score():
            best_motifs = motifs
    return best_motifs
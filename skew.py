#!/usr/local/bin/python3

from sys import argv

import array

import numpy

skew = {'G':1, 'C':-1, 'A':0, 'T':0}

def compute_skew(genome):
    skew_array = numpy.zeros(len(genome)+1, numpy.int)
    for i in range(0, len(genome)):
        skew_array[i+1] = skew_array[i] + skew[genome[i]]
    return skew_array

genome_skew = compute_skew(argv[1])
print(genome_skew)
print('minimum skew value of',numpy.amin(genome_skew),'found at position',numpy.argmin(genome_skew),sep=' ')

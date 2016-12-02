#!/usr/local/bin/python3

from sys import argv

import numpy

import argparse

parser = argparse.ArgumentParser(description='Find genome lowest skew values.')
parser.add_argument('--file', help='the file that contains a sequence of nucleotide symbols')
parser.add_argument('--genome', help='genome string of nucleotide symbols')

skew = {'G':1, 'C':-1, 'A':0, 'T':0}

def compute_skew(genome):
    skew_array = numpy.zeros(len(genome)+1, numpy.int)
    for i in range(0, len(genome)):
        skew_array[i+1] = skew_array[i] + skew[genome[i]]
    return skew_array


args = parser.parse_args()

print(args)

if args.file != None:
    with open(args.file, 'r') as myfile: data=myfile.read().replace('\n', '')
else:
    data = args.genome

genome_skew = compute_skew(data)
print(genome_skew)
min_skew_index = numpy.argmin(genome_skew)
ii = numpy.where(genome_skew == genome_skew[min_skew_index])[0]
print(ii)
print('minimum skew value of', numpy.amin(genome_skew), 'found at position', ii ,sep=' ')

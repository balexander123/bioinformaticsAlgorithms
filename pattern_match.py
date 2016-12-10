#!/usr/local/bin/python3

from sys import argv

import argparse

parser = argparse.ArgumentParser(description='Find frequent k-mer words with maximum allowed mismatches d.')
parser.add_argument('--file', help='the file that contains a sequence of nucleotide symbols')
parser.add_argument('--genome', help='genome string of nucleotide symbols')
parser.add_argument('--pattern', help='pattern to search within file or genome')

args = parser.parse_args()

def pattern_match(text, pattern):
    offset_str=""
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            offset_str+=str(i)+" "
    return offset_str

if args.file != None:
    with open(args.file, 'r') as myfile: data=myfile.read().replace('\n', '')
else:
    data = args.genome

print(pattern_match(data, args.pattern))

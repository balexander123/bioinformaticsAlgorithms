#!/usr/local/bin/python3

import sys

from sys import argv

nucleotideIndex = {'A': 0, 'C': 1, 'G': 2, 'T': 3 }

def k_mer_patterns_from_text(text, k):
    k_mer_patterns = []
    for i in range(0, len(text) - k + 1):
        k_mer_patterns.append(text[i:i + k])
    return k_mer_patterns


def profile_most_probable_k_mer(text, k, profile_matrix):
    most_probable_k_mer = ''
    most_probable = 0
    for k_mer in k_mer_patterns_from_text(text, k):
        probability = 1
        for col in range(0,k):
            nucletide = k_mer[col:col+1]
            row = nucleotideIndex[nucletide]
            prob = profile_matrix[nucleotideIndex[nucletide]][col]
            probability *= profile_matrix[nucleotideIndex[nucletide]][col]
        if probability > most_probable:
            most_probable = probability
            most_probable_k_mer = k_mer
    return most_probable_k_mer

def profile_matrix_from_file(file):
    profile_matrix = []
    with open(file, 'r') as input_file:
        text = input_file.readline()
        k = int(input_file.readline())
        for profile_line in sys.stdin.readlines():
            profile_matrix.append([float(p) for p in profile_line.split()])
    return profile_matrix

profile_matrix = []

if len(argv) == 2:
    with open(argv[1], 'r') as input_file:
        text = input_file.readline()
        k = int(input_file.readline())
        for profile_line in input_file.readlines():
            profile_matrix.append([float(p) for p in profile_line.split()])
else:
    text = sys.stdin.readline()
    k = int(sys.stdin.readline())
    profile_matrix = []
    for profile_line in sys.stdin.readlines():
        profile_matrix.append([float(p) for p in profile_line.split()])

print(profile_most_probable_k_mer(text[0:len(text)-1],k,profile_matrix))
#!/usr/local/bin/python3

#
# https://stepik.org/lesson/Some-Hidden-Messages-are-More-Surprising-than-Others-3/step/2?course=Bioinformatics-Algorithms&unit=7
#
#
# Reverse Complement Problem: Find the reverse complement of a DNA string.
#      Input: A DNA string Pattern.
#      Output: Patternrc , the reverse complement of Pattern.

from sys import argv

base_complement = {'A':'T','C':'G','T':'A','G':'C'}

def reverse_complement(pattern):
    reverse_comp = ""
    for c in pattern[::-1]:
        reverse_comp+=base_complement[c]
    return reverse_comp

print(reverse_complement(argv[1]))
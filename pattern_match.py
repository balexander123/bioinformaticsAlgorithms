#!/usr/local/bin/python3

from sys import argv

def pattern_match(pattern, text):
    offset_str=""
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            offset_str+=str(i)+" "
    return offset_str

print(pattern_match(argv[1], argv[2]))

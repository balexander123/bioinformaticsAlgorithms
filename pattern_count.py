#!/usr/local/bin/python3

from sys import argv

def pattern_count(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

print(pattern_count(argv[1], argv[2]))

#!/usr/local/bin/python3

from sys import argv
import array

def pattern_count(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

def frequent_words(text, k):
    frequent_patterns = []
    count=array.array('i',(0,)*len(text))
    for i in range(len(text) - k):
        pattern = text[i:i+k]
        count[i] = pattern_count(text,pattern)
    maxCount = max(count)
    for i in range(len(text) - k):
        if count[i] == maxCount:
            frequent_patterns.append(text[i:i+k])
    return ' '.join(set(frequent_patterns))

print(frequent_words(argv[1], int(argv[2])))
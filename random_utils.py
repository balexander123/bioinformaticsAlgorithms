#!/usr/local/bin/python3

import numpy
import random

def random_distribution(probabilities):
    result_dist = numpy.zeros(len(probabilities),numpy.float)
    sum_prob = sum(probabilities)
    for i in range(0, len(probabilities)):
        result_dist[i] = probabilities[i] / sum_prob
    return result_dist

def roll(distribution):
    rand_roll = random.uniform(0, sum(distribution))
    dist_sum = 0
    result = 1
    for dist in distribution:
        dist_sum += dist
        if rand_roll < dist_sum:
            return result
        result += 1
    return result

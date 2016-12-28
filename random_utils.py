#!/usr/local/bin/python3

import numpy

def random_distribution(probabilities):
    result_dist = numpy.zeros(len(probabilities),numpy.float)
    sum_prob = sum(probabilities)
    for i in range(0, len(probabilities)):
        result_dist[i] = probabilities[i] / sum_prob
    return result_dist

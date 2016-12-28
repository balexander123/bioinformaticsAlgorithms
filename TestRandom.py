import unittest
import random_utils

class TestRandom(unittest.TestCase):

    def setUp(self):
        self.probabilites = [0.1, 0.2, 0.3]

    def test_result_size(self):
        self.assertEqual(len(random_utils.random_distribution(self.probabilites)),3,"resulting probability distribution length is equal to number of probabilites")

    def test_results(self):
        results = random_utils.random_distribution(self.probabilites)
        self.assertEqual(len(results),3,"resulting probability distribution length is equal to number of probabilites")
        self.assertTrue(abs(sum(results)-1) < 0.001, "resulting sum of probabilities is equal to 1")
        self.assertTrue(abs(results[0]-1/6) < 0.001,"probability distribution 0 is equal to 1/6")
        self.assertTrue(abs(results[1]-1/3) < 0.001,"probability distribution 1 is equal to 1/3")
        self.assertTrue(abs(results[2]-1/2) < 0.001,"probabality distribution 2 is equal to 1/2")

if __name__ == '__main__':
    unittest.main()
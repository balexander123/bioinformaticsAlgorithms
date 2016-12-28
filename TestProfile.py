import unittest
import profile


class TestProfile(unittest.TestCase):
    NUMBER_OF_MOTIFS = 5
    K_MER_LENGTH = 8

    def setUp(self):
        file = 'fixtures/week4_1_1_test_data_set.txt'
        quiz_file = 'fixtures/week4_quiz_test_data.txt'
        self.profile = profile.Profile(file)
        self.profile_quiz = profile.Profile(quiz_file)

    def test_profile_attributes(self):
        self.assertEqual(self.profile.num_motifs, self.NUMBER_OF_MOTIFS, "There are 5 dna motifs")
        self.assertEqual(self.profile.k, self.K_MER_LENGTH, "k-mer length is 8")

    def test_motif_matrix_entries(self):
        self.assertEqual(''.join(self.profile.motif_matrix[0]),
                         'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                         'first entry of motif matrix is CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA')
        self.assertEqual(''.join(self.profile.motif_matrix[self.NUMBER_OF_MOTIFS - 1]),
                         'AATCCACCAGCTCCACGTGCAATGTTGGCCTA',
                         'last entry of motif matrix is AATCCACCAGCTCCACGTGCAATGTTGGCCTA')

    def test_profile_matrix_dimensions(self):
        self.assertEqual(len(self.profile.profile_matrix), 4, 'profile matrix has 4 rows')
        self.assertEqual(len(self.profile.profile_matrix[0]), 32, 'profile matrix has 32 columns')

    def test_score(self):
        score = self.profile.motif_score()
        self.assertEqual(score,83,'Score(Motifs) is equal to 83')

    # CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
    # GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
    # TAGTACCGAGACCGAAAGAAGTATACAGGCGT
    # TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
    # AATCCACCAGCTCCACGTGCAATGTTGGCCTA

    def test_laplace_rule_of_succession(self):
        self.assertEqual(' '.join(str(i) for i in self.profile.count_matrix[0]),
                         '2 4 1 2 2 3 2 2 3 2 2 1 1 2 3 3 3 1 2 3 3 2 3 3 3 2 2 1 2 2 2 3','Count(Motif) for A')
        self.assertEqual(' '.join(str(i) for i in self.profile.count_matrix[1]),
                         '2 1 2 4 3 4 3 3 1 2 2 2 4 2 1 2 1 1 2 3 3 2 1 2 1 2 2 1 3 5 3 2','Count(Motif) for C')
        self.assertEqual(' '.join(str(i) for i in self.profile.count_matrix[2]),
                         '2 3 4 1 2 1 2 3 2 3 2 3 2 4 3 2 3 4 2 2 2 3 1 2 3 3 2 6 3 1 2 2','Count(Motif) for G')
        self.assertEqual(' '.join(str(i) for i in self.profile.count_matrix[3]),
                         '3 1 2 2 2 1 2 1 3 2 3 3 2 1 2 2 2 3 3 1 1 2 4 2 2 2 3 1 1 1 2 2','Count(Motif) for T')

    def test_consensus(self):
        pass

    def test_compute_profile_matrix(self):
        for column in range(len(self.profile.profile_matrix[0])):
            profile_column = [self.profile.profile_matrix[i][column] for i in range(4)]
            self.assertAlmostEqual(sum(profile_column),1.0,7,'Column probabilities sum to 1')


if __name__ == '__main__':
    unittest.main()

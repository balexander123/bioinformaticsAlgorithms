import unittest
import profile
import randomized_motif_search

class TestRandomizedMotifSearch(unittest.TestCase):

    def setUp(self):
        self.randomized_motif_search = randomized_motif_search.RandomizedMotifSearch()

    def test_quiz_most_probable(self):
        dna_strings = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
        quiz_profile = profile.Profile('fixtures/week4_quiz_starting_motifs.txt')
        most_probable_k_mers = self.randomized_motif_search.most_probable(quiz_profile,dna_strings)
        self.assertEqual(' '.join(most_probable_k_mers),"TGA GTT GAA TGT", 'Computes most probable k-mer pattern in dna')

    def test_whether_i_am_insane(self):
        alt_file = 'fixtures/bioinformatics_book_p92.txt'
        self.alt_profile = profile.Profile(alt_file)
        dna_strings = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", 'CCCTAACGAG', "CGTCAGAGGT"]
        self.assertEqual(' '.join(self.randomized_motif_search.most_probable(self.alt_profile, dna_strings)),
                         "ACCT ATGT GCGT ACGA AGGT", 'Computes most probable k-mer pattern in dna')


if __name__ == '__main__':
    unittest.main()

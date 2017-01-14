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
        self.assertEqual(' '.join(most_probable_k_mers),"TGA TAA GGA TGT", 'Computes most probable k-mer pattern in dna')

    def test_bioinformatics_book_p92(self):
        alt_file = 'fixtures/bioinformatics_book_p92.txt'
        self.alt_profile = profile.Profile(alt_file)
        dna_strings = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", 'CCCTAACGAG', "CGTCAGAGGT"]
        self.assertEqual(' '.join(self.randomized_motif_search.most_probable(self.alt_profile, dna_strings)),
                         "ACCT ATGT GCGT ACGA AGGT", 'Computes most probable k-mer pattern in dna')

    def test_there_are_t_initial_random_kmers(self):
        self.random_profile = profile.Profile('1_1_5_Randomized_Motif_Search.txt')
        random_motifs = self.randomized_motif_search._random_k_mers(self.random_profile.k,self.random_profile.dna)
        self.assertEqual(len(random_motifs),5,'The correct number of k-mers was returned')

    def test_1000_iterations(self):
        with open('1_1_5_Randomized_Motif_Search.txt', 'r') as input_file:
            first_line = input_file.readline()
            k = int(first_line.split()[0])
            t = int(first_line.split()[1])
            dna = []
            for dna_string in input_file.readlines():
                dna.append(dna_string)
        for i, text in enumerate(dna):
            dna[i] = text.replace('\n', '')
        best_motifs = self.randomized_motif_search.randomized_motif_searchX(k, dna, 1000)
        #                                       TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG
        self.assertEqual(' '.join(best_motifs),'TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG','Best motifs after 1000 iterations')
        # python3 ./randomized_motif_search.py 1_1_5_Randomized_Motif_Search.txt
        # TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG

    def test_extra_dataset(self):
        with open('ExtraDataSet_.txt', 'r') as input_file:
            first_line = input_file.readline()
            k = int(first_line.split()[0])
            t = int(first_line.split()[1])
            dna = []
            for dna_string in input_file.readlines():
                dna.append(dna_string)
        for i, text in enumerate(dna):
            dna[i] = text.replace('\n', '')
        best_motifs = self.randomized_motif_search.randomized_motif_searchX(k, dna, 1000)
        #                                       TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG
        self.assertEqual(' '.join(best_motifs),'CATGGGGAAAACTGA CCTCTCGATCACCGA CCTATAGATCACCGA CCGATTGATCACCGA CCTTGTGCAGACCGA CCTTGCCTTCACCGA CCTTGTTGCCACCGA ACTTGTGATCACCTT CCTTGTGATCAATTA CCTTGTGATCTGTGA CCTTGTGATCACTCC AACTGTGATCACCGA CCTTAGTATCACCGA CCTTGTGAAATCCGA CCTTGTCGCCACCGA TGTTGTGATCACCGC CACCGTGATCACCGA CCTTGGTTTCACCGA CCTTTGCATCACCGA CCTTGTGATTTACGA','Best motifs after 1000 iterations')
        # python3 ./randomized_motif_search.py 1_1_5_Randomized_Motif_Search.txt
        # TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG


if __name__ == '__main__':
    unittest.main()

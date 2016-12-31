import unittest
import motif_enumeration

class TestMotifEnumeration(unittest.TestCase):

    def setUp(self):
        file='fixtures/stepic_1_2_Motif_Finding.txt'
        self.motif_enum = motif_enumeration.MotifEnumeration(file)

    def flatten(self, s):
        return ' '.join(str(e) for e in sorted(list(s)))

    def test_neighbors_with_distance_zero(self):
        self.assertEqual(self.flatten(self.motif_enum.neighbors("AAA", 0)), "AAA", 'Hamming distance of zero returns same pattern')

    def test_neighbors_with_pattern_length_one(self):
        self.assertEqual(self.flatten(self.motif_enum.neighbors("A", 1)), "A C G T", 'Pattern length 1 returns single nucleotides')

    def test_neighbors_with_length_two(self):
        self.assertEqual(self.flatten(self.motif_enum.neighbors("AA", 1)), "AA AC AG AT CA GA TA", 'Pattern length 2 neighbors')

    def test_neighbors_stepic_1_16(self):
        self.assertEqual(self.flatten(self.motif_enum.neighbors("ACG", 1)), "AAG ACA ACC ACG ACT AGG ATG CCG GCG TCG", 'Pattern length 3 neighbors')

    def test_3_mers_in_dna(self):
        self.assertEqual(' '.join(self.motif_enum.k_mer_patterns_from_dna(self.motif_enum.dna,3)),"ATT TTT TTG TGG GGC TGC GCC CCT CTT TTA CGG GGT GTA TAT ATC GAA AAA AAA AAT ATT",'All 3-mer patterns found in dna')

    def test_appears_in(self):
        motifs = ['ATA', 'ATT', 'GTT', 'TTT']
        for motif in motifs:
            self.assertTrue(self.motif_enum.appears_in(self.motif_enum.dna, motif, 1),'All motifs appear in dna with at most 1 mis-match')

    def test_motif_enumeration(self):
        self.assertEqual(' '.join(sorted(self.motif_enum.motif_enumeration())),"ATA ATT GTT TTT",'Motif enumeration should work')

if __name__ == '__main__':
    unittest.main()
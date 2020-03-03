from analysis import markers_comparison as mc
import unittest
import pandas as pd


class TestMarkersComparison(unittest.TestCase):

    def test_get_matching_alleles(self):
        df = pd.DataFrame({'sample1': [1, 2, 1], 'sample2': [1, 1, 1], 'informative': [1, 1, 2]})

        # get match by allele 1
        r = mc.get_match_by_specific_allele(df['sample1'], df['sample2'], df['informative'], 1)
        self.assertEqual(r[0], True)
        self.assertEqual(r[1], False)
        self.assertEqual(r[2], False)

        # get match by allele 2
        r = mc.get_match_by_specific_allele(df['sample1'], df['sample2'], df['informative'], 2)
        self.assertEqual(r[0], False)
        self.assertEqual(r[1], False)
        self.assertEqual(r[2], False)

        # get match by both
        r = mc.get_matching_alleles(df, True)
        self.assertEqual(r[0], True)
        self.assertEqual(r[1], False)
        self.assertEqual(r[2], False)

        # if we don't filter informative, then the third marker should be a match
        r = mc.get_matching_alleles(df, False)
        self.assertEqual(r[0], True)
        self.assertEqual(r[1], False)
        self.assertEqual(r[2], True)


if __name__ == '__main__':
    unittest.main()
import refilter
from cyvcf2 import VCF
import unittest
import os


class TestFilter(unittest.TestCase):
    def setUp(self):
        # load test data
        # store each variant object into specific variables for tes
        test_directory = os.path.dirname(os.path.abspath(__file__))
        reader = VCF(os.path.join(test_directory, "test.vcf"))
        self.variants = [ variant for variant in reader ]
        self.test_filter = refilter.Filter(0.3, 0.7, 'AB', 5, -5.0, ['MISSING'], ['DB'])

    def test_filters_ok(self):
        self.assertTrue(self.test_filter.filters_ok(self.variants[0]))
        self.assertFalse(self.test_filter.filters_ok(self.variants[2]))
        self.assertTrue(self.test_filter.filters_ok(self.variants[3]))

    def test_info_fields_ok(self):
        self.assertFalse(self.test_filter.info_fields_ok(self.variants[0]))
        self.assertTrue(self.test_filter.info_fields_ok(self.variants[1]))

    def test_pass_variant(self):
        for variant in self.variants:
            self.test_filter.pass_variant(variant)
            self.assertEqual(variant.FILTER, None)

    def test_rescue(self):
        expected_filters = (None, None, 'VQSRTrancheSNP99.00to99.90;MISSING', 'VQSRTrancheSNP99.00to99.90', None)
        for index, variant in enumerate(self.variants):
            self.test_filter.rescue(variant)
            self.assertEqual(variant.FILTER, expected_filters[index])


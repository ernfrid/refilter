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
        self.test_filter = refilter.Filter(0.3, 0.7, 'AB', 'VAR_DP', 5, ['MISSING'], ['DB'])
        reader.add_filter_to_header(self.test_filter.filtered_header())
        reader.add_info_to_header(self.test_filter.rescued_header())

        self.variants = [ variant for variant in reader ]

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
        expected_filters = (
                self.test_filter.filter_tag,
                self.test_filter.filter_tag,
                ';'.join(('VQSRTrancheSNP99.00to99.90;MISSING', self.test_filter.filter_tag)),
                None,
                None
                )
        for index, variant in enumerate(self.variants):
            self.test_filter.rescue(variant)
            self.assertEqual(variant.FILTER, expected_filters[index])


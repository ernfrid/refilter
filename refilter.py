#!/usr/bin/env python

import click
from cyvcf2 import VCF, Writer
import numpy as np

class Filter(object):
    def __init__(self, min_allele_balance, max_allele_balance, allele_balance_tag, min_depth, min_vqslod, exclude_filters, exclude_fields):
        self.allele_bounds = (min_allele_balance, max_allele_balance)
        self.allele_balance_tag = allele_balance_tag
        self.min_depth = min_depth
        self.min_vqslod = min_vqslod
        self.exclude_filters = exclude_filters
        self.exclude_fields = exclude_fields
        self.filter_tag = 'AB_FILTERED{0}to{1}'.format(self.allele_bounds[0], self.allele_bounds[1])
        self.rescue_tag = 'AB_RESCUED'

    def rescued_header(self):
        return {
                'ID' : self.rescue_tag,
                'Description' : 'Filter status before rescued based on allele balance',
                'Type': 'String',
                'Number' : '1'
                }

    def _filter_description(self):
        desc = 'Failed allele balance filter. {1} <= {0} <= {2} and VQSLOD >= {3} and DP >= {4}.'.format(self.allele_balance_tag, self.allele_bounds[0], self.allele_bounds[1], self.min_vqslod, self.min_depth)
        if self.exclude_filters is not None:
            desc += ' Ignored sites with FILTER containing: {0}.'.format(','.join(self.exclude_filters))
        if self.exclude_fields is not None:
            desc += ' Ignored sites with any of the following tags in the INFO field: {0}.'.format(','.join(self.exclude_fields))
        return desc


    def filtered_header(self):
        return {
                'ID' : self.filter_tag,
                'Description' : self._filter_description()
                }

    def filters_ok(self, variant):
        if variant.FILTER is None:
            return True
        else:
            filters = set(variant.FILTER.split(';'))
            for filt_string in self.exclude_filters:
                if filt_string in filters:
                    return False
            return True

    def info_fields_ok(self, variant):
        info_tags = set([ t[0] for t in variant.INFO ])
        for field in self.exclude_fields:
            if field in info_tags:
                return False
        return True

    def __call__(self, variant):
        if self.filters_ok(variant) and self.info_fields_ok(variant):
            self.rescue(variant)

    def rescue(self, variant):
        ab = variant.INFO[self.allele_balance_tag]
        variant_filter = variant.FILTER
        if (ab >= self.allele_bounds[0] and
                ab <= self.allele_bounds[1] and
                variant.INFO['DP'] >= self.min_depth and
                variant.INFO['VQSLOD'] >= self.min_vqslod):
            if variant.FILTER is not None:
                variant.INFO[self.rescue_tag] = variant.FILTER
            self.pass_variant(variant)
        else:
            # FAILED
            self.fail_variant(variant)

    @staticmethod
    def pass_variant(variant):
        variant.FILTER = 'PASS'

    def fail_variant(self, variant):
        if variant.FILTER is not None:
            variant.FILTER = variant.FILTER.split(';') + [self.filter_tag]
        else:
            variant.FILTER = self.filter_tag

@click.command()
@click.option('--min-allele-balance', default=0.3, type=click.FLOAT,
        help='Minimum allele balance value to rescue a variant')
@click.option('--max-allele-balance', default=0.7, type=click.FLOAT,
        help='Maximum allele balance value to rescue a variant')
@click.option('--allele-balance-tag', default='AB_HOM', type=click.STRING,
        help='INFO containing the source of allele balance information')
@click.option('--min-depth', default=10, type=click.INT,
        help='Minimum depth to rescue a variant')
@click.option('--min-vqslod', default=None, type=click.FLOAT,
        help='Minimum VQSLOD score to allow rescue')
@click.option('--exclude-filters', default=['MISSING'], multiple=True, type=click.STRING,
        help='Filters that cannot be rescued')
@click.option('--exclude-fields', default=['OLD_MULTIALLELIC'], multiple=True, type=click.STRING,
        help='INFO fields that preclude a line from being rescued')
@click.argument('vcf', type=click.Path())
def main(min_allele_balance,
        max_allele_balance,
        allele_balance_tag,
        min_depth,
        min_vqslod,
        exclude_filters,
        exclude_fields,
        vcf):
    reader = VCF(vcf)
    refilter = Filter(min_allele_balance,
            max_allele_balance,
            allele_balance_tag,
            min_depth,
            min_vqslod,
            exclude_filters,
            exclude_fields)
    reader.add_filter_to_header(refilter.filtered_header())
    reader.add_info_to_header(refilter.rescued_header())
    writer = Writer('-', reader)

    for variant in reader:
        refilter(variant) # Modifies variant filter status in place
        writer.write_record(variant)

if __name__ == '__main__':
    main()


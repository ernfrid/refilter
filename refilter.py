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
        if (ab >= self.allele_bounds[0] and
                ab <= self.allele_bounds[1] and
                variant.INFO['DP'] >= self.min_depth and
                variant.INFO['VQSLOD'] >= self.min_vqslod):
            self.pass_variant(variant)

    @staticmethod
    def pass_variant(variant):
        variant.FILTER = 'PASS'

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
    writer = Writer('-', reader)
    refilter = Filter(min_allele_balance,
            max_allele_balance,
            allele_balance_tag,
            min_depth,
            min_vqslod,
            exclude_filters,
            exclude_fields)

    for variant in reader:
        refilter(variant) # Modifies variant filter status in place
        writer.write_record(variant)

if __name__ == '__main__':
    main()


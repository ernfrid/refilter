#!/usr/bin/env python

import click
from cyvcf2 import VCF, Writer
import numpy as np

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

    for variant in reader:
        filter_field = variant.FILTER
        if filter_field is None:
            continue
        filters = set(variant.FILTER.split(';'))
        for field in exclude_fields:
            if field in variant.INFO:
                writer.write_record(variant)
                continue
        for filt_string in exclude_filters:
            if filt_string in filters:
                writer.write_record(variant)
                continue # FIXME THIS IS A BUG
        ab = variant.INFO[allele_balance_tag]
        if (ab >= min_allele_balance and
                ab <= max_allele_balance and
                variant.INFO['DP'] >= min_depth and
                variant.INFO['VQSLOD'] >= min_vqslod):
            variant.FILTER = 'PASS'
            writer.write_record(variant)

if __name__ == '__main__':
    main()


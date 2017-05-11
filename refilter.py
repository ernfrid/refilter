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
def main(min_allele_balance,
        max_allele_balance,
        allele_balance_tag,
        min_depth,
        min_vqslod,
        exclude_filters,
        exclude_fields):
    print  ','.join(exclude_filters)

if __name__ == '__main__':
    main()


#!/usr/bin/env python

from __future__ import division
import sys

class AlleleBalanceAccumulator(object):
    def __init__(self):
        self.ref = 0
        self.alt = 0

    def __call__(self, AD):
        depths = AD.split(',')
        self.ref += int(depths[0])
        self.alt += int(depths[1])

def calculate_balance(ref, alt):
    try:
        ab = alt / (ref + alt)
    except ZeroDivisionError:
        ab = 0
    return '{0:.4f}'.format(ab)

if __name__ == '__main__':
    for line in sys.stdin:
        if line.startswith('##'):
            print line,
        elif line.startswith('#'):
            print '##INFO=<ID=HetAB,Number=1,Type=Float,Description="heterozygous genotype allele balance">'
            print '##INFO=<ID=HetAB_DP,Number=1,Type=Integer,Description="heterozygous genotype read depth">' 
            print '##INFO=<ID=HetHomAltAB,Number=1,Type=Float,Description="heterozygous + homozygous ALT genotype allele balance">'
            print '##INFO=<ID=HetHomAltAB_DP,Number=1,Type=Integer,Description="heterozygous + homozygous ALT genotype read depth">'
            print line,
        else:
            fields = line.rstrip().split('\t')
            format_fields = fields[8].split(':')
            gtindex = format_fields.index('GT')
            adindex = format_fields.index('AD')
            maxindex = max(gtindex, adindex)
            stats = { '0/1': AlleleBalanceAccumulator(), '1/1': AlleleBalanceAccumulator() }
            for sample in fields[9:]:
                sample_fields = sample.split(':', maxindex + 1)
                try:
                    stats[sample_fields[gtindex]](sample_fields[adindex])
                except KeyError:
                    pass

            reftotal = stats['0/1'].ref + 0.5 * stats['1/1'].alt + stats['1/1'].ref
            alttotal = stats['0/1'].alt + 0.5 * stats['1/1'].alt

            ab_info_fields = (
                    fields[7],
                    'HetAB={0}'.format(calculate_balance(stats['0/1'].ref, stats['0/1'].alt)),
                    'HetAB_DP={0}'.format(stats['0/1'].ref + stats['0/1'].alt),
                    'HetHomAltAB={0}'.format(calculate_balance(reftotal, alttotal)),
                    'HetHomAltAB_DP={0}'.format(int(reftotal + alttotal))
                    )
            fields[7] = ';'.join( ab_info_fields )
            print '\t'.join(fields)

#!/usr/bin/env python

import cyvcf2
import argparse

parser = argparse.ArgumentParser(description='Print the genotype for all individuals at a given SNP', add_help=True)
parser.add_argument('vcf', help='VCF path')
parser.add_argument('snp_id', help='SNP ID')
args = parser.parse_args()

f = cyvcf2.VCF(args.vcf, 'rb')

def format_genotype(g):
    delimiter = '|' if g[-1] else '/'
    return delimiter.join([str(i) for i in g[0:2]])

for variant in f:
    if variant.ID == args.snp_id:
        for sample, genotype in zip(f.samples, variant.genotypes):
            print('{}\t{}'.format(sample, format_genotype(genotype)))
        break

f.close()

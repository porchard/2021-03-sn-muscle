#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import pandas as pd
import numpy as np
import glob
import csv
import logging
import argparse
import gzip

parser = argparse.ArgumentParser(description='Add rsIDs to MAGIC files.', add_help = True)
parser.add_argument('magic', type = str,  help = 'MAGIC summary stats file (gzipped).')
parser.add_argument('out', type = str,  help = 'File to create')
parser.add_argument('--vcf', type = str, default='/lab/data/reference/human/hg19/annot/dbsnp150_variants/All_20170710.vcf.gz', help = 'dbSNP VCF file (default: /lab/data/reference/human/hg19/annot/dbsnp150_variants/All_20170710.vcf.gz).')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

# load in the diamante data. Determine the ref and alt alleles...
logging.info('Reading MAGIC GWAS file')
diamante = pd.read_csv(args.magic, delimiter='\t')
diamante['Chr'] = diamante.MarkerName.map(lambda x: x.split(':')[0])
diamante['Pos'] = diamante.MarkerName.map(lambda x: int(x.split(':')[1]))
diamante['variant_type'] = diamante.MarkerName.map(lambda x: x.split(':')[2])
diamante = diamante[diamante.variant_type=='SNP']
diamante.Allele1 = diamante.Allele1.map(lambda x: x.upper())
diamante.Allele2 = diamante.Allele2.map(lambda x: x.upper())

# map chrom:pos:ref:alt -> rsid
chrom_pos = set([i for i in (diamante.Chr.map(str) + ':' + diamante.Pos.map(str))])
logging.info('Looking for information for {} variants'.format(len(chrom_pos)))
chrom_pos_ref_alt_to_rsid = dict()
logging.info('Scanning through VCF to find variants')

VCF_HEADER_RE = re.compile('^#')

line_count = 0
with gzip.open(args.vcf, 'rt') as f:
    for line in f:
        line_count += 1
        if line_count % 5000000 == 0:
            logging.info('Processed {} lines of the VCF'.format(line_count))
        if VCF_HEADER_RE.match(line):
            continue
        line = line.rstrip().split()
        chrom, pos, rsid, ref, alt = line[0:5]
        if '{chrom}:{pos}'.format(**locals()) in chrom_pos:
            for i in alt.split(','):
                key = '{}:{}:{}:{}'.format(chrom, pos, ref, i)
                chrom_pos_ref_alt_to_rsid[key] = rsid
logging.info('Stored information for {} variants'.format(len(chrom_pos_ref_alt_to_rsid)))

logging.info('Determining the rsIDs')
conversions = []
count = 0
converted_count = 0
for index, row in diamante.iterrows():
    count += 1
    if count % 1000000 == 0:
        logging.info('Converted {} of {} SNPs so far'.format(converted_count, count))
    option_1 = '{Chr}:{Pos}:{Allele1}:{Allele2}'.format(**row)
    option_2 = '{Chr}:{Pos}:{Allele2}:{Allele1}'.format(**row)
    if option_1 not in chrom_pos_ref_alt_to_rsid and option_2 not in chrom_pos_ref_alt_to_rsid:
        # logging.info('No conversion for SNP at {Chr}:{Pos}'.format(**row))
        conversions.append('')
        continue
    conversion = chrom_pos_ref_alt_to_rsid[option_1] if option_1 in chrom_pos_ref_alt_to_rsid else chrom_pos_ref_alt_to_rsid[option_2]
    conversions.append(conversion)
    converted_count += 1

diamante['rsID'] = conversions

converted_count = sum([1 for i in conversions if i != ''])
logging.info('Successfully converted {} of {} variants'.format(converted_count, len(conversions)))
diamante = diamante.rename(columns={'TotalSampleSize': 'N'})

noconversion = diamante[diamante.rsID != '']
diamante = diamante[~diamante.rsID.isnull()]
diamante.to_csv(args.out, sep='\t', index=False)
noconversion.to_csv(args.out + '.dropped', sep='\t', index=False)

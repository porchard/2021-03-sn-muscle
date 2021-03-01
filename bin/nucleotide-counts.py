#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pysam
import argparse
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Given a file of positions (chrom, pos), output the nucleotide counts (nA, nC, ...) from a BAM file at those positions.', add_help=True)
parser.add_argument('--min-base-quality', dest='min_base_quality', type=int, default=0, help='Min. base quality (default: 0)')
parser.add_argument('--min-mapping-quality', dest='min_mapping_quality', type=int, default=0, help='Min. mapping quality (default: 0)')
parser.add_argument('--min-depth', dest='min_depth', type=int, default=0, help='Min. depth to output (default: 0)')
parser.add_argument('positions', help='File of positions for which to fetch nucleotide counts (chrom, pos; tsv).')
parser.add_argument('bam', help='BAM file.')
args = parser.parse_args()

NUCLEOTIDES = ['A', 'C', 'G', 'T']

logging.info('Reading positions')
positions = set()
with open(args.positions, 'r') as f_pos:
    for line in f_pos:
        chrom, pos = line.rstrip().split()
        positions.add('{}:{}'.format(chrom, pos))


logging.info('Read {} positions'.format(len(positions)))
position_count = 0
print('\t'.join(['chrom', 'pos', 'depth'] + ['n{}'.format(i) for i in NUCLEOTIDES]))
with pysam.AlignmentFile(args.bam, 'rb') as f_bam:
    for pileupcolumn in f_bam.pileup(max_depth=9999999, min_base_quality=args.min_base_quality, min_mapping_quality=args.min_mapping_quality, ignore_overlaps=True):
        chrom, pos = (pileupcolumn.reference_name, pileupcolumn.reference_pos + 1) # pos is 0-based in pysam, 1-based in VCF
        var = '{}:{}'.format(chrom, pos)
        if var in positions:
            position_count += 1
        else:
            continue
        if position_count % 10000 == 0:
            logging.info('Processed {} positions'.format(position_count))
        # NOTE: pileupcolumn.nsegments ignores the base quality filter. Therefore, not using it here.
        nsegments = 0
        nucleotide_counts = [0 for i in NUCLEOTIDES]
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                nsegments += 1
                nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                if not nuc in NUCLEOTIDES:
                    raise ValueError('At position {}:{}, encountered a nucleotide that is not in [{}]: {}'.format(chrom, pos, ', '.join(NUCLEOTIDES), nuc))
                nucleotide_counts[NUCLEOTIDES.index(nuc)] += 1
        if nsegments < args.min_depth:
            continue
        print('\t'.join([str(i) for i in [chrom, pos, nsegments] + nucleotide_counts]))

logging.info('Processed {} positions'.format(position_count))
logging.info('Done.')

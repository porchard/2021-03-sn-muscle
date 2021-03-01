#!/usr/bin/env python
# coding: utf-8

import sys
import logging
import pysam
import argparse

parser = argparse.ArgumentParser(description='Get fragment ends from an ATAC-seq bam file. Output is NOT sorted.', add_help=True)
parser.add_argument('bam', help='Bam file. Does not need to be sorted.')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

count = 0

with pysam.AlignmentFile(args.bam, 'rb') as f:
    for read in f.fetch(until_eof=True):
        count += 1
        if count % 1000000 == 0:
            logging.info('Processed {} reads'.format(count))
        fragment_end = None
        if read.is_read1:
            fragment_end = [read.reference_name, read.reference_start, read.reference_start + 1]
        elif read.is_read2:
            fragment_end = [read.reference_name, read.reference_end-1, read.reference_end]
        if fragment_end is None:
            logging.warning('Fragment end is None??')
            sys.exit()
        fragment_end = [str(i) for i in fragment_end]
        print('\t'.join(fragment_end + [read.query_name]))

logging.info('Done.')

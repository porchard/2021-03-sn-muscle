#!/usr/bin/env python

import sys
import os
import argparse
import pysam
import re
import logging

parser = argparse.ArgumentParser(description='Count number of times each barcode is observed (CB tag)')
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

counts = dict()
processed = 0

with pysam.AlignmentFile(args.bam_in, 'rb') as old_bam:
    for read in old_bam.fetch(until_eof = True):
        processed += 1
        if processed % 1000000 == 0:
            logging.info('Processed {} reads'.format(processed))
        if not read.has_tag('CB'):
            continue
        barcode = read.get_tag('CB')
        if barcode not in counts:
            counts[barcode] = 0
        counts[barcode] += 1

for barcode, count in sorted(counts.items()):
    print('{}\t{}'.format(barcode, count))

logging.info('Done')

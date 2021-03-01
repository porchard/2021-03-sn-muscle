#!/usr/bin/env python

import os
import gzip
import sys
import logging
import argparse
import re

parser = argparse.ArgumentParser(description='Make a file of TSS from a GTF.', add_help = True)
parser.add_argument('--attribute', type = str, default='gene_name', help = 'GTF attribute to use as name (default: gene_name)')
parser.add_argument('gtf', type = str,  help = 'GTF file (gzipped).')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

with gzip.open(args.gtf, 'rt') as f:
    for line in f:
        if re.match('^#', line):
            continue
        seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split('\t')
        if feature != 'transcript':
            continue
        tss = int(start) if strand == '+' else int(end)-1
        attributes = {x: y for x, y in [i.split() for i in attribute.split('; ')]}
        print('\t'.join([seqname, str(tss-1), str(tss), attributes[args.attribute].replace('"', ''), '.', strand]))

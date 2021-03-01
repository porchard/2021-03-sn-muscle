#!/usr/bin/env python

import os
import sys
import gzip

BED = sys.argv[1]

mappings = dict()
current_id = 0

with gzip.open(BED, 'rt') as f:
    for line in f:
        line = line.rstrip().split('\t')
        read_name = line[3]
        read_name_parts = read_name.split('_')
        read_name_no_barcodes = read_name_parts[0]
        if read_name_no_barcodes not in mappings:
            current_id += 1
            mappings[read_name_no_barcodes] = 'R' + str(current_id)
        read_name_parts[0] = mappings[read_name_no_barcodes]
        line[3] = '_'.join(read_name_parts)
        print('\t'.join(line))


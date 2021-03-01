#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from scipy.sparse import csc_matrix
import h5py
import argparse

parser = argparse.ArgumentParser(description='Convert a file of counts (library, barcode, feature, count) to hdf5 format, as expected for LIGER', add_help=True)
parser.add_argument('--include-features', dest='include_features', type=str, default=None, help='List of features to include')
parser.add_argument('counts', help='File of counts')
parser.add_argument('hdf5_out', help='Output file')
args = parser.parse_args()

KEEP_FEATURES = set()
if args.include_features is not None:
    with open(args.include_features, 'r') as f:
        for line in f:
            KEEP_FEATURES.add(line.rstrip())

counts = []
with open(args.counts, 'r') as f:
    for line in f:
        library, barcode, feature, count = line.rstrip().split()
        if args.include_features is not None and feature not in KEEP_FEATURES:
            continue
        count = int(count)
        counts.append([library, barcode, feature, count])

counts = pd.DataFrame(counts, columns=['library', 'barcode', 'feature', 'count'])
counts['nucleus'] = counts.library + '-' + counts.barcode

counts_wide = counts[['nucleus', 'feature', 'count']].pivot(index='feature', columns='nucleus', values='count').fillna(0)
M = csc_matrix(counts_wide.values)

f = h5py.File(args.hdf5_out, 'w')
g = f.create_group('matrix')
g.create_dataset('barcodes', data=[str(i).encode('ascii', 'ignore') for i in counts_wide.columns.values], dtype="S40")
g.create_dataset('data', data=M.data)
g.create_dataset('indptr', data=M.indptr)
g.create_dataset('indices', data=M.indices)
g.create_dataset('shape', data=M.shape)
g.attrs['shape'] = M.shape
g = f.create_group('matrix/features')
g.create_dataset('name', data=[str(i).encode('ascii', 'ignore') for i in counts_wide.index.values], dtype="S40")
f.close()

#!/usr/bin/env python

import os
import sys
import csv
import numpy as np

# QC file is the file with barcode, number_genes, number_umis, fraction_mitochondrial
# RNA nulcie is the file listing library and barrcode
# library is the library of interest
QC_FILE, RNA_NUCLEI, LIBRARY = sys.argv[1:]

barcodes = set()
umi_counts = list()

with open(RNA_NUCLEI, 'r') as f:
    for line in f:
        library, barcode, individual = line.rstrip().split()
        if library == LIBRARY:
            barcodes.add(barcode)

with open(QC_FILE, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for line in reader:
        if line['barcode'] in barcodes:
            umi_counts.append(int(line['number_umis']))

umi_counts = np.array(umi_counts)
print('Library: {}\t\tMean UMI counts: {}\t\tMedian UMI counts: {}'.format(LIBRARY, np.mean(umi_counts), np.median(umi_counts)))

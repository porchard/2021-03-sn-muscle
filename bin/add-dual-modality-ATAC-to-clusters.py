#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd

#CLUSTERS = '/lab/work/porchard/daily/2020-12-21/labs/clusters.txt'
CLUSTERS = sys.argv[1]

ATAC_BARCODE_LIST = '/home/porchard/github/snATACseq-NextFlow/737K-arc-v1.txt'
RNA_BARCODE_LIST = '/home/porchard/github/snRNAseq-NextFlow/737K-arc-v1.txt'
DUAL_MODALITY_RNA_LIBRARIES = ['1846_RNA-hg19', '1846_RNA-rn6']

# map RNA --> ATAC barcodes
atac_barcodes = pd.read_csv(ATAC_BARCODE_LIST, header=None, names=['atac_barcode'])
rna_barcodes = pd.read_csv(RNA_BARCODE_LIST, header=None, names=['rna_barcode'])
barcodes = pd.concat([atac_barcodes, rna_barcodes], axis=1)
barcode_rna_to_atac = dict(zip(barcodes.rna_barcode, barcodes.atac_barcode))


clusters = pd.read_csv(CLUSTERS, sep='\t', header=None, names=['library', 'barcode', 'cluster'])

# get dual modality libraries
dual_modality_rna_clusters = clusters[clusters.library.isin(DUAL_MODALITY_RNA_LIBRARIES)]
dual_modality_atac_clusters = dual_modality_rna_clusters.copy()
dual_modality_atac_clusters.library = dual_modality_atac_clusters.library.map(lambda x: x.replace('RNA', 'ATAC'))
dual_modality_atac_clusters.barcode = dual_modality_atac_clusters.barcode.map(barcode_rna_to_atac)

# sanity check to make sure the ATAC dual modality clusters aren't already there...
for i in dual_modality_atac_clusters.library.unique():
    assert(i not in clusters.library.unique())


clusters = pd.concat([clusters, dual_modality_atac_clusters])
clusters.to_csv(sys.stdout, sep='\t', index=False, header=None)

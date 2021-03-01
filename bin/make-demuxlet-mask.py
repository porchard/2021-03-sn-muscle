#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Create a BED file to make highly expressed genes')
parser.add_argument('--mask-top', dest='mask_top', type=float, default=0.3, help='Fraction of top genes (by counts) to mask. default: 0.3 (top 30% of genes)')
parser.add_argument('starsolo_features', help='features.tsv file from starsolo')
parser.add_argument('starsolo_mtx', help='matrix.mtx file from starsolo')
parser.add_argument('gtf', help='GTF file')
args = parser.parse_args()


def parse_gtf_attributes(a):
    """
    a should be a string, e.g.: 'gene_id "ENSG00000223972.5"; transcript_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";'
    """
    x = [i.split(' ') for i in a.rstrip(';').split('; ')]
    d = {key: val.replace('"', '') for key, val in x}
    return d

#STARSOLO_FEATURES = '/lab/work/porchard/sn-muscle-project/work/rnaseq/results/starsolo/63_40_rna-hg19/Solo.out/GeneFull/raw/features.tsv'
#STARSOLO_MAT = '/lab/work/porchard/sn-muscle-project/work/rnaseq/results/starsolo/63_40_rna-hg19/Solo.out/GeneFull/raw/matrix.mtx'
#GTF = '/lab/data/reference/human/hg19/index/STAR/current/gencode.v19.annotation.gtf'

feature_index_to_feature = dict()
feature_id_to_feature_name = dict()
feature_index = 0
with open(args.starsolo_features, 'r') as f:
    for line in f:
        feature_index += 1
        gene_id, gene_name = line.rstrip().split('\t')
        feature_id_to_feature_name[gene_id] = gene_name
        feature_index_to_feature[feature_index] = gene_id



counts = pd.read_csv(args.starsolo_mtx, skiprows=3, dtype=int, sep=' ', header=None, names=['feature_index', 'barcode_index', 'umis'])
counts['feature_id'] = counts.feature_index.map(feature_index_to_feature)

gene_counts = counts.groupby('feature_id').umis.sum().to_frame(name='umis').reset_index()
gene_counts['feature_name'] = gene_counts.feature_id.map(feature_id_to_feature_name)
gene_counts = gene_counts.sort_values('umis', ascending=False)
gene_counts['quant'] = np.arange(1, len(gene_counts)+1) / len(gene_counts)

MASK_FEATURES = set(gene_counts[gene_counts.quant<=args.mask_top].feature_id)

gtf = pd.read_csv(args.gtf, header=None, sep='\t', names=['chrom', 'source', 'type', 'start', 'end', 'x1', 'strand', 'x2', 'attributes'])
gtf = gtf[(gtf.type=='gene')]
gtf['gene_id'] = gtf.attributes.map(lambda x: parse_gtf_attributes(x)['gene_id'])

mask = gtf.loc[gtf.gene_id.isin(MASK_FEATURES),['chrom', 'start', 'end']].sort_values(['chrom', 'start', 'end'])
mask.to_csv(sys.stdout, sep='\t', index=False, header=False)

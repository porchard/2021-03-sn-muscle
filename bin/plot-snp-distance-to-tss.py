#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import pybedtools as bt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#SNP_BED = '/lab/work/porchard/sn-muscle-project/cicero-snps.bed'
#GENCODE_TSS = '/lab/work/porchard/data/gencode-tss/hg19.bed'

# SNP_BED should be BED4, with the name in column 4
# TSS should be BED6
SNP_BED, TSS = sys.argv[1:]

nearby_tss = bt.BedTool(SNP_BED).sort().closest(bt.BedTool(TSS).sort(), k=300, d=True).to_dataframe()
nearby_tss = nearby_tss[['chrom', 'start', 'end', 'name', 'thickEnd', 'blockSizes']].rename(columns={'blockSizes': 'distance', 'thickEnd': 'gene'})

# make sure for each SNP we get all TSS within 1Mb
assert(nearby_tss.groupby('name').distance.max().min() > 1e6)

close = nearby_tss[nearby_tss.distance<=1e6]
close = close.groupby(['name', 'gene']).distance.min().reset_index()

# for each SNP, plot distance to TSS for each gene
for snp in close.name.unique():
    df = close[close.name==snp].sort_values('distance')
    df['idx'] = range(len(df))
    df['label'] = [f'{i} ({j/1000} kb)' for i, j in zip(df.gene, df.distance)]
    n_genes = len(df)
    fig, ax = plt.subplots(figsize=(n_genes/3,5))
    ax.scatter(x='idx', y='distance', data=df)
    ax.set_ylim(bottom=0, top=1e6)
    ax.set_xticks(df.idx)
    ax.set_xticklabels(df.label, rotation=90)
    ax.set_ylabel(f'Distance from {snp} to closest TSS')
    fig.tight_layout()
    fig.savefig(f'{snp}.distance-to-tss.png')
    fig.clf()

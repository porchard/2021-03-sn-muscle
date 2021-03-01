#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import cmap_utils

#UMAP = '/lab/work/porchard/sn-muscle-project/work/liger-post-factorization/results/umap.txt'
#CLUSTERING = '/lab/work/porchard/sn-muscle-project/work/process-by-cluster-2021-01-20/clusters.txt'
#LIBRARY_LABELS = '/lab/work/porchard/sn-muscle-project/library-labels.txt'

UMAP, CLUSTERING, LIBRARY_LABELS = sys.argv[1:]

library_labels = pd.read_csv(LIBRARY_LABELS, sep='\t')
library_to_modality = dict(zip(library_labels.library, library_labels.modality))

umap = pd.read_csv(UMAP, header=None, sep='\t', names=['library', 'barcode', 'dim1', 'dim2'])

clusters = pd.read_csv(CLUSTERING, header=None, sep='\s+', names=['library', 'barcode', 'cluster'])

both = umap.merge(clusters, on=['library', 'barcode'])
both['is_dual_modality_atac'] = both.library.map(lambda x: '1846' in x and 'ATAC' in x)
both = both[~both.is_dual_modality_atac]
both['species'] = both.library.map(lambda x: x.split('-')[-1])
both['modality'] = both.library.map(lambda x: library_to_modality[x] if x in library_to_modality else x.split('-')[0].split('_')[1])

# main UMAP
cmap = cmap_utils.make_colormap_dict(both.cluster.sort_values().unique().tolist())

fig, ax = plt.subplots(figsize=(3,3))
g = sns.scatterplot(x='dim1', y='dim2', hue='cluster', data=both, ax=ax, palette=cmap, edgecolor=None, alpha=0.1, s=10)
g.get_legend().remove()
ax.set_xlabel('UMAP dim. 1')
ax.set_ylabel('UMAP dim. 2')

fig.savefig('umap.png', bbox_inches='tight')
fig.clf()


# UMAP by species and modality
cmap = cmap_utils.make_colormap_dict(both.cluster.sort_values().unique().tolist())

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(4,4), gridspec_kw={'hspace': 0.15, 'wspace': 0.02})
HUMAN_ALPHA = 0.15
RAT_ALPHA = 0.25
s = 5

# human RNA
ax = axs[0,0]
g = sns.scatterplot(x='dim1', y='dim2', hue='cluster', data=both[(both.species=='hg19') & (both.modality=='RNA')], ax=ax, palette=cmap, edgecolor=None, alpha=HUMAN_ALPHA, s=s)
g.get_legend().remove()
ax.set_xlabel('')
ax.set_xticks([])
ax.set_ylabel('UMAP dim. 2')
ax.set_title('Human RNA')

# rat RNA
ax = axs[0,1]
g = sns.scatterplot(x='dim1', y='dim2', hue='cluster', data=both[(both.species=='rn6') & (both.modality=='RNA')], ax=ax, palette=cmap, edgecolor=None, alpha=RAT_ALPHA, s=s)
g.get_legend().remove()
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('Rat RNA')

# human ATAC
ax = axs[1,0]
g = sns.scatterplot(x='dim1', y='dim2', hue='cluster', data=both[(both.species=='hg19') & (both.modality=='ATAC')], ax=ax, palette=cmap, edgecolor=None, alpha=HUMAN_ALPHA, s=s)
g.get_legend().remove()
ax.set_xlabel('UMAP dim. 1')
ax.set_ylabel('UMAP dim. 2')
ax.set_title('Human ATAC')

# rat ATAC
ax = axs[1,1]
g = sns.scatterplot(x='dim1', y='dim2', hue='cluster', data=both[(both.species=='rn6') & (both.modality=='ATAC')], ax=ax, palette=cmap, edgecolor=None, alpha=RAT_ALPHA, s=s)
g.get_legend().remove()
ax.set_yticks([])
ax.set_xlabel('UMAP dim. 1')
ax.set_ylabel('')
ax.set_title('Rat ATAC')

fig.savefig('umap-by-species-and-modality.png', bbox_inches='tight')
fig.clf()


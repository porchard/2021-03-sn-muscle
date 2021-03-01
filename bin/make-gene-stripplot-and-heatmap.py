#!/usr/bin/env python
# coding: utf-8

# In[92]:


import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clusters', help='Path to file of nucleis --> cluster')
parser.add_argument('--cluster-names', dest='cluster_names')
parser.add_argument('--genes', nargs='+', help='List of genes to plot')
parser.add_argument('--matrices', nargs='+', help='List of HDF5 files of gene counts')
args = parser.parse_args()

#CLUSTERS = '/lab/work/porchard/sn-muscle-project/work/process-by-cluster/clusters.txt'
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
#GENES = ['PITX2', 'MEF2A', 'MEF2B', 'MEF2C', 'MEF2D', 'FOS', 'JUN']
#MATRICES = ['/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/data/KSM2_RNA.hdf5']

CLUSTERS = args.clusters
CLUSTER_NAMES = args.cluster_names
GENES = args.genes
MATRICES = args.matrices


# In[101]:


def load_gene_and_total_counts(f, genes):
    tmp = pd.read_hdf(f)
    total_counts = tmp.sum(axis=1).reset_index().rename(columns={0: 'total'})
    gene_counts = tmp[[i for i in genes if i in tmp.columns]].reset_index()
    both = total_counts.merge(gene_counts, on='nucleus')
    return both


# In[102]:


clusters = pd.read_csv(CLUSTERS, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clusters['nucleus'] = clusters.library + '-' + clusters.barcode
nucleus_to_cluster = dict(zip(clusters.nucleus, clusters.cluster))
clusters.head()


# In[103]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_to_cluster_names = dict(zip(cluster_names.old_name, cluster_names.new_name))
cluster_to_cluster_names


# In[104]:


counts = pd.concat([load_gene_and_total_counts(f, GENES) for f in MATRICES])
counts.head()


# In[105]:


counts['cluster'] = counts.nucleus.map(nucleus_to_cluster)
counts.head()


# In[106]:


for GENE in GENES:
    if GENE not in counts.columns:
        continue
    counts_per_cluster = counts[['total', GENE, 'cluster']].groupby('cluster').sum().reset_index()
    counts_per_cluster['cpm'] = 1e6 * counts_per_cluster[GENE] / counts_per_cluster.total
    counts_per_cluster.sort_values('cluster')
    
    # plot:
    # stripplot (boxplot?)
    # heatmap (CPM)
    # colorbar
    fig, axs = plt.subplots(nrows=3, figsize=(5,3), gridspec_kw={'height_ratios':[0.3, 0.5, 5], 'hspace': 0.1})
    #fig, axs = plt.subplots(nrows=2, gridspec_kw={'height_ratios':[5, 0.5], 'hspace': 0.1})

    # make stripplot
    ax = axs[2]
    sns.stripplot(x='cluster', y=GENE, data=counts, color='black', alpha=0.1, ax=ax)
    ax.set_xlabel('')
    ax.set_xticklabels(counts_per_cluster.cluster.map(cluster_to_cluster_names), rotation=90)
    ax.set_ylabel(f'{GENE} counts')

    # heatmap
    ax = axs[1]
    fcbar = ax.imshow(counts_per_cluster[['cpm']].T, cmap='Reds', vmin=0, aspect='auto')
    ax.set_yticks([])
    ax.set_xticks([])

    # colorbar
    cbar = plt.colorbar(fcbar, cax=axs[0], orientation='horizontal')
    axs[0].set_xlabel(f'{GENE} CPM')
    axs[0].xaxis.tick_top()
    axs[0].xaxis.set_label_position('top')
   
    fig.savefig(f'{GENE}.png', bbox_inches='tight')
    fig.clf()
    


# In[ ]:





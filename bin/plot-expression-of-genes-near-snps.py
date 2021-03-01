#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import pybedtools as bt
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import numpy as np


#TSS = '/lab/work/porchard/data/gencode-tss/hg19.bed'
SNPS = ['chr5:53271420', 'chr12:26472562']
USE_SIGNED_DISTANCE = True
DISTANCE_TO_SHOW = 1e6
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
#CLUSTERS = '/lab/work/porchard/sn-muscle-project/work/process-by-cluster/clusters.txt'
#MATRICES = ['/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/data/KSM2_RNA.hdf5']
TSS = sys.argv[1]
CLUSTERS = sys.argv[2]
CLUSTER_NAMES = sys.argv[3]
MATRICES = sys.argv[4:]


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')#.astype(str)
cluster_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))

def get_tss_near_snps(snp_list, tss_file, distance):
    """
    snp_list = [chr5:10000, chr5:11111, ...]
    tss_file = bed6
    distance = bps from SNP to check
    """
    snps = pd.DataFrame([i.split(':') for i in snp_list], columns=['chrom', 'end'])
    snps.end = snps.end.astype(int)
    snps['start'] = snps.end - 1
    snps_bt = bt.BedTool().from_dataframe(snps[['chrom', 'start', 'end']]).sort()

    nearby_tss = snps_bt.window(bt.BedTool(tss_file).sort(), w=distance).to_dataframe()
    nearby_tss['distance'] = (nearby_tss.end - nearby_tss.strand).abs()
    nearby_tss['distance'] = -1*(nearby_tss.end - nearby_tss.strand)
    if not USE_SIGNED_DISTANCE:
        nearby_tss['distance'] = (nearby_tss.end - nearby_tss.strand).abs()
    nearby_tss['snp'] = nearby_tss.chrom + ':' + nearby_tss.end.astype(str)
    nearby_tss = nearby_tss[['snp', 'thickStart', 'distance']].rename(columns={'thickStart': 'gene'})
    return nearby_tss


def load_gene_and_total_counts(f, genes):
    tmp = pd.read_hdf(f)
    total_counts = tmp.sum(axis=1).reset_index().rename(columns={0: 'total'})
    gene_counts = tmp[[i for i in genes if i in tmp.columns]].reset_index()
    both = total_counts.merge(gene_counts, on='nucleus')
    return both


# In[44]:


clusters = pd.read_csv(CLUSTERS, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clusters['nucleus'] = clusters.library + '-' + clusters.barcode
nucleus_to_cluster = dict(zip(clusters.nucleus, clusters.cluster))
clusters.head()


# In[45]:


# first, get all TSS w/in 1Mb of the SNP. Note the distance to closest TSS
tss_near_snps = None
if USE_SIGNED_DISTANCE:
    tss_near_snps = get_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW).groupby(['snp', 'gene']).distance.agg(lambda x: x[x.abs()==x.abs().min()]).reset_index()
else:
    tss_near_snps = get_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW).groupby(['snp', 'gene']).distance.min().reset_index()
tss_near_snps.head()


# In[46]:


# now, get the expression of each gene in each cluster
GENES = tss_near_snps.gene.unique().tolist()
counts = pd.concat([load_gene_and_total_counts(f, GENES) for f in MATRICES])
counts['cluster'] = counts.nucleus.map(nucleus_to_cluster)
counts.head()


# In[47]:


# get CPMs
totals = counts[GENES + ['total', 'cluster']].groupby('cluster').sum().reset_index()
cpm = totals.melt(id_vars=['cluster', 'total'], var_name='gene', value_name='counts')
cpm['cpm'] = 1e6*cpm.counts / cpm.total
cpm = cpm[['cluster', 'gene', 'cpm']]
cpm.head()


# In[52]:


for SNP in tss_near_snps.snp.unique():
    df = tss_near_snps[tss_near_snps.snp==SNP].merge(cpm, on='gene').sort_values('distance')
    df['label'] = df[['gene', 'distance']].apply(lambda x: '{} ({} kb)'.format(x[0], x[1]/1000), axis
    =1)
    label_order = df[['label', 'distance']].drop_duplicates().sort_values('distance').label.to_list()
    df = df.pivot(index='label', values='cpm', columns='cluster')
    df = df.loc[label_order,:]
    df.columns = df.columns.to_series().map(cluster_to_cluster_name)
    #df_norm = df.transform(lambda x: x/max(x) if max(x) > 0 else x, axis=1)
    df_norm = df.applymap(lambda x: np.log(x+1))
    COLOR_LIMIT = max([df_norm.applymap(abs).max().max(), 0.2])
    fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(3, len(df)/3), gridspec_kw={'hspace': 0.02, 'height_ratios': [1, len(df)]})
    cbar_ax = axs[0]
    ax = axs[1]
    my_cm = cm.get_cmap('Reds', 256)
    my_cm.set_bad(color='grey')
    #norm = matplotlib.colors.Normalize(vmin=-1.,vmax=1.)
    fcbar = ax.imshow(df_norm, aspect='auto', cmap=my_cm, vmin=0, vmax=COLOR_LIMIT)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df.index.to_list())
    ax.set_ylim(top=-0.5, bottom=ax.get_ylim()[0]+0.5)

    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns.to_list(), rotation=90)

    # add annotations
    for r_index in range(len(df)):
        for c_index in range(len(df.columns)):
            if not np.isnan(df.iat[r_index,c_index]):
                ax.text(x=c_index, y=r_index, s=round(df_norm.iat[r_index,c_index], 1), va='center', ha='center')

    # add colorbar
    plt.colorbar(fcbar, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel('Log(CPM + 1)')
    cbar_ax.xaxis.tick_top()
    cbar_ax.xaxis.set_label_position('top')
    
    fig.savefig('{}.png'.format(SNP).replace(':', '_'), bbox_inches='tight')
    fig.clf()

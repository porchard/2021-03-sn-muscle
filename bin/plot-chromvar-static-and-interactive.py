#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import re

#DEVIATIONS = '/lab/work/porchard/sn-muscle-project/work/chromvar/results/chromvar/deviation-zscores.txt' #sys.argv[1]
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
#CLUSTERS = '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/liger/round-1/second-louvain-clusters.txt'
#MOTIF_NAMES = '/lab/work/porchard/sn-muscle-project/data/cisbp-motifs/motif-id-to-tf-name.txt'

DEVIATIONS, CLUSTERS, CLUSTER_NAMES, MOTIF_NAMES = sys.argv[1:]
LABEL_MOTIFS = ['5710', '6341', '4556', '4570', '4595', '3988', '5867', '3589', '5292', '6355', '3913', '5863'] #  '4541']
LABEL_MOTIFS = ['M{}_1.02'.format(i) for i in LABEL_MOTIFS]

def tissue_specificity_index(l):
    """
    Calculate tissue specificity index (as in Yanai et al 2015: https://pubmed.ncbi.nlm.nih.gov/15388519/)

    Input: list of gene expression value across tissues (one element per tissue)
    Output: tissue specificity index (float)

    Example: [0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0] --> 0.95
    """
    if not isinstance(l, list) and not isinstance(l, np.ndarray) and not isinstance(l, pd.Series):
        raise TypeError('l must be a list, numpy array, or pandas series')
    a = np.array(l)
    if a.max() == 0:
        return 0
    else:
        return sum(1 - (a / a.max())) / (len(l) - 1)


motif_names = pd.read_csv(MOTIF_NAMES, sep='\t', header=None, names=['motif', 'name'])
motif2name = {motif: name for motif, name in zip(motif_names.motif, motif_names.name)}


deviations = pd.read_csv(DEVIATIONS, sep='\t')

cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t', dtype=str).rename(columns={'old_name': 'cluster', 'new_name': 'cluster_name'})

clusters = pd.read_csv(CLUSTERS, sep='\t', header=None, names=['library', 'barcode', 'cluster'], dtype=str)
clusters['nucleus'] = clusters.library + '-' + clusters.barcode
clusters = clusters[['nucleus', 'cluster']].merge(cluster_names, on='cluster')
clusters = clusters[['nucleus', 'cluster_name']].rename(columns={'cluster_name': 'cluster'})

deviations = deviations.merge(clusters, on='nucleus')

per_cluster = deviations.groupby(['motif', 'cluster']).deviation_z.mean().reset_index()
per_cluster['motif_name'] = per_cluster.motif.map(lambda x: motif2name[x] if x in motif2name else '???')
per_cluster['tf_score'] = per_cluster.groupby('motif').deviation_z.transform(lambda x: x/max(x))
per_cluster.loc[per_cluster.deviation_z<0,'tf_score'] = 0


motif_specificity = per_cluster[['motif', 'tf_score']].groupby('motif').agg(tsi=pd.NamedAgg(column='tf_score', aggfunc=tissue_specificity_index)).reset_index()
per_cluster = per_cluster.merge(motif_specificity, on='motif')


# make one bigger interactive heatmap
SHOW_MOTIFS = per_cluster[(per_cluster.deviation_z>=1) & (per_cluster.tsi>=0.7)].motif.unique()
for_heatmap = per_cluster[per_cluster.motif.isin(SHOW_MOTIFS)]
for_heatmap.motif = for_heatmap.motif_name + ' (' + for_heatmap.motif + ')'
for_heatmap = for_heatmap[['motif', 'cluster', 'tf_score']].pivot(index='motif', columns='cluster', values='tf_score')
g = sns.clustermap(for_heatmap, cmap='viridis', figsize=(10, 30))
fig = px.imshow(g.data2d.T.loc[cluster_names.sort_values('cluster').cluster_name.to_list(),:], color_continuous_midpoint=0.5, color_continuous_scale='Reds', origin='lower', labels={'color': 'TF score'}, aspect='auto')
fig.update_layout(title='Motif scores', template='plotly_white', yaxis_title='', width=3000, height=800)
fig.write_html('tf-scores-full.html')


#SHOW_MOTIFS = per_cluster[(per_cluster.deviation_z>=1.85) & (per_cluster.tsi>=0.85)].motif.unique()
SHOW_MOTIFS = per_cluster[(per_cluster.deviation_z>=2) & (per_cluster.tsi>=0.85)].motif.unique()
for_heatmap = per_cluster[per_cluster.motif.isin(SHOW_MOTIFS)]
for_heatmap.motif = for_heatmap.motif_name + ' (' + for_heatmap.motif + ')'
for_heatmap = for_heatmap[['motif', 'cluster', 'tf_score']].pivot(index='motif', columns='cluster', values='tf_score')
g = sns.clustermap(for_heatmap, cmap='viridis', figsize=(10, 20))
fig = px.imshow(g.data2d.T.loc[cluster_names.sort_values('cluster').cluster_name.to_list(),:], color_continuous_midpoint=0.5, color_continuous_scale='Reds', origin='lower', labels={'color': 'TF score'}, aspect='auto')
fig.update_layout(title='Motif scores', template='plotly_white', yaxis_title='', width=2000, height=800)
fig.write_html('tf-scores.html')


g = sns.clustermap(for_heatmap, cmap='viridis', figsize=(10, 20), yticklabels=1)
new_labels = []
for t in g.ax_heatmap.get_yticklabels():
    motif_id = re.match('.*\((.*)\)', t.get_text()).group(1)
    if motif_id not in LABEL_MOTIFS:
        new_labels.append('')
    else:
        new_labels.append(t.get_text())
g.ax_heatmap.set_yticklabels(new_labels)


df = g.data2d.T
df = df.loc[cluster_names.sort_values('cluster').cluster_name.to_list(),:] # re-order cell types
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(8, 3), gridspec_kw={'wspace': 0.05, 'width_ratios': [1, 0.05]})
ax = axs[0]
cbar_ax = axs[1]
fcbar = ax.imshow(df, aspect='auto')
ax.set_yticks(range(len(cluster_names)))
ax.set_ylim(bottom=ax.get_ylim()[0]+0.5, top=ax.get_ylim()[1]-0.5)
ax.set_yticklabels(df.index.to_list())
ax.set_xticks([count for count, i in enumerate(new_labels) if i != ''])
ax.set_xticklabels([i for i in new_labels if i != ''], rotation=45, ha='right')
# add colorbar
plt.colorbar(fcbar, cax=cbar_ax, orientation='vertical')
cbar_ax.set_ylabel('Motif score')
fig.savefig('static-heatmap.png', bbox_inches='tight')
fig.clf()


# make per-motif stripplots
tmp = deviations.copy()
tmp['motif_name'] = tmp.motif.map(lambda x: motif2name[x] if x in motif2name else '???')
tmp = tmp[tmp.motif.isin(LABEL_MOTIFS)]
tmp.motif = tmp.motif_name + ' (' + tmp.motif + ')'
for motif, df in tmp.groupby('motif'):
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.stripplot(y='cluster', x='deviation_z', data=df, alpha=0.3, orient='h', color='black', ax=ax)
    ax.set_ylabel('')
    ax.set_xlabel('chromVAR deviation z-score')
    ax.axvline(x=0, linestyle='dashed', color='red')
    ax.set_title(motif)
    FILE_NAME = motif.replace(' ', '_').replace('(', '').replace(')', '') + '.png'
    fig.tight_layout()
    fig.savefig(FILE_NAME)
    fig.clf()

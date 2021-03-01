#!/usr/bin/env python
# coding: utf-8

# In[11]:


import os
import sys
import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#LIGER_IN_MATRIX_GLOB = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/results/long-matrices/*'
#CLUSTERING = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/results/find-promising/KSM2_RNA/interesting-clustering.15___0.05.txt'
#LIGER_IN_MATRICES = glob.glob(LIGER_IN_MATRIX_GLOB)

CLUSTERING = sys.argv[1]
OUT = sys.argv[2]
LIGER_IN_MATRICES = sys.argv[3:]

def feature_count_file_to_cluster_counts(f, nucleus_to_cluster):
    counts = {}
    SKIPPED = set()
    with open(f, 'r') as fh:
        for line in fh:
            nucleus, feature, count = line.rstrip().split('\t')
            count = int(count)
            if nucleus not in nucleus_to_cluster:
                SKIPPED.add(nucleus)
                continue
            cluster = nucleus_to_cluster[nucleus]
            if cluster not in counts:
                counts[cluster] = dict()
            if feature not in counts[cluster]:
                counts[cluster][feature] = 0
            counts[cluster][feature] += count
    if len(SKIPPED) > 0:
        for nucleus in SKIPPED:
            sys.stderr.write('WARNING: SKipped nucleus {} in file {} (cluster unknown)\n'.format(nucleus, f))
    return counts

clustering = pd.read_csv(CLUSTERING, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clustering['library'] = clustering[['library', 'barcode']].apply(lambda x: x[0].replace(x[1], ''), axis=1)
clustering['nucleus'] = clustering.library + '-' + clustering.barcode
clustering.cluster = clustering.cluster.astype(int)
nucleus_to_cluster = dict(zip(clustering.nucleus, clustering.cluster))


# In[3]:


count_matrices = dict()
for f in LIGER_IN_MATRICES:
    sys.stderr.write(f'{f}\n')
    individual, modality = os.path.basename(f).split('.')[0].split('_')
    count_matrices[(individual, modality)] = feature_count_file_to_cluster_counts(f, nucleus_to_cluster)


counts = []
for key in count_matrices:
    individual, modality = key
    for cluster in count_matrices[key]:
        for gene in count_matrices[key][cluster]:
            counts.append([individual, modality, cluster, gene, count_matrices[key][cluster][gene]])
counts = pd.DataFrame(counts, columns=['individual', 'modality', 'cluster', 'feature', 'count'])
counts['species'] = counts.individual.map(lambda x: 'rat' if 'rat' in x else 'human')


# In[6]:


counts_per_species = counts.groupby(['species', 'modality', 'cluster', 'feature'])['count'].sum().reset_index()


# In[13]:


counts_per_species['cpm'] = counts_per_species.groupby(['species', 'modality', 'cluster'])['count'].transform(lambda x: 1e6*x/sum(x))
counts_per_species.feature = counts_per_species.feature.map(lambda x: x.upper())

MARKER_GENES = ['MYH1', 'MYH7', 'PDGFRA', 'VWF', 'ACTA2', 'PTPRC', 'PAX7']
cpm = counts_per_species.loc[counts_per_species.feature.isin(MARKER_GENES),['species', 'modality', 'cluster', 'feature', 'cpm']].sort_values('cluster')


# In[53]:


fig = plt.figure(constrained_layout=True, figsize=(5,4))
widths = [1, 1, 0.05]
heights = [1, 1]
gs = fig.add_gridspec(ncols=3, nrows=2, width_ratios=widths,
                          height_ratios=heights)
cbar_ax = fig.add_subplot(gs[:,2])
human_rna_ax = fig.add_subplot(gs[0,0])
human_atac_ax = fig.add_subplot(gs[1,0])
rat_rna_ax = fig.add_subplot(gs[0,1])
rat_atac_ax = fig.add_subplot(gs[1,1])
ax = np.array([[human_rna_ax, rat_rna_ax],
     [human_atac_ax, rat_atac_ax]])

species = ['human', 'rat']
for i in species:
    rna = cpm[(cpm.species==i) & (cpm.modality=='RNA')]
    if len(rna) > 0:
        rna = rna.pivot(index='cluster', columns='feature', values='cpm').fillna(0)
        rna = rna.transform(lambda x: x/max(x)).loc[:,[i for i in MARKER_GENES if i in rna.columns]]
        rna_ax = ax[0,species.index(i)]
        rna_ax.set_title('{} RNA'.format(i.capitalize()))
        if species.index(i) == 0:
            sns.heatmap(rna, ax=rna_ax, cmap='Reds', cbar_ax=cbar_ax)
            cbar_ax.set_ylabel('Relative expression or accessibility', size=12)
        else:
            sns.heatmap(rna, ax=rna_ax, cmap='Reds', cbar=None)
        rna_ax.set_ylim(top=rna_ax.get_ylim()[1]-0.5)
        rna_ax.set_ylim(bottom=rna_ax.get_ylim()[0]+0.5)
        rna_ax.xaxis.set_ticks([])
        rna_ax.set_xlabel('')
        for tmp in ['top', 'bottom','right', 'left']:
            rna_ax.spines[tmp].set_visible(True)
        # remove y-axis labels if not left-most panel
        if species.index(i) != 0:
            rna_ax.yaxis.set_ticks([])
            rna_ax.set_ylabel('')
        else:
            rna_ax.set_ylabel('Cluster', size=12)

    atac = cpm[(cpm.species==i) & (cpm.modality=='ATAC')]
    if len(atac) > 0:
        atac = atac.pivot(index='cluster', columns='feature', values='cpm').fillna(0)
        atac = atac.transform(lambda x: x/max(x)).loc[:,[i for i in MARKER_GENES if i in atac.columns]]
        atac_ax = ax[1,species.index(i)]
        atac_ax.set_title('{} ATAC'.format(i.capitalize()))
        sns.heatmap(atac, ax=atac_ax, cmap='Reds', cbar=None)
        atac_ax.set_ylim(top=atac_ax.get_ylim()[1]-0.5)
        atac_ax.set_ylim(bottom=atac_ax.get_ylim()[0]+0.5)
        atac_ax.set_xlabel('')
        for tmp in ['top', 'bottom','right', 'left']:
            atac_ax.spines[tmp].set_visible(True)
        # remove y-axis labels if not left-most panel
        if species.index(i) != 0:
            atac_ax.yaxis.set_ticks([])
            atac_ax.set_ylabel('')
        else:
            atac_ax.set_ylabel('Cluster', size=12)
fig.savefig(OUT)


# In[ ]:





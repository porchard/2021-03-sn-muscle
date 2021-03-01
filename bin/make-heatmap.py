#!/usr/bin/env python
# coding: utf-8

# In[38]:


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

# tissue_specificity_index([0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0])


# In[39]:


def feature_count_file_to_cluster_counts(f, nucleus_to_cluster):
    counts = {}
    with open(f, 'r') as fh:
        for line in fh:
            nucleus, feature, count = line.rstrip().split('\t')
            count = int(count)
            cluster = nucleus_to_cluster[nucleus]
            if cluster not in counts:
                counts[cluster] = dict()
            if feature not in counts[cluster]:
                counts[cluster][feature] = 0
            counts[cluster][feature] += count
    return counts


# In[40]:


clustering = pd.read_csv(CLUSTERING, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clustering['nucleus'] = clustering.library + '-' + clustering.barcode
clustering.cluster = clustering.cluster.astype(int)
nucleus_to_cluster = dict(zip(clustering.nucleus, clustering.cluster))


# In[41]:


count_matrices = dict()
for f in LIGER_IN_MATRICES:
    sys.stderr.write(f'{f}\n')
    individual, modality = os.path.basename(f).split('.')[0].split('_')
    count_matrices[(individual, modality)] = feature_count_file_to_cluster_counts(f, nucleus_to_cluster)


# In[19]:


counts = []
for key in count_matrices:
    individual, modality = key
    for cluster in count_matrices[key]:
        for gene in count_matrices[key][cluster]:
            counts.append([individual, modality, cluster, gene, count_matrices[key][cluster][gene]])
counts = pd.DataFrame(counts, columns=['individual', 'modality', 'cluster', 'feature', 'count'])
counts.head()


# In[22]:


counts['cpm'] = counts.groupby(['individual', 'modality', 'cluster'])['count'].transform(lambda x: 1e6*x/sum(x))
counts.feature = counts.feature.map(lambda x: x.upper())
counts.head()


# In[25]:


MARKER_GENES = ['MYH1', 'MYH2', 'MYH4', 'MYH7', 'PDGFRA', 'VWF', 'ACTA2', 'CD163', 'PAX7']
cpm = counts.loc[counts.feature.isin(MARKER_GENES),['individual', 'modality', 'cluster', 'feature', 'cpm']].sort_values('cluster')
cpm.head()


# In[46]:


individuals = list(sorted(cpm.individual.unique().tolist()))
fig, ax = plt.subplots(ncols=len(individuals), nrows=2, figsize=(3*len(individuals), 6))
for i in individuals:
    rna = cpm[(cpm.individual==i) & (cpm.modality=='RNA')]
    if len(rna) > 0:
        rna = rna.pivot(index='cluster', columns='feature', values='cpm').fillna(0)
        rna = rna.transform(lambda x: x/max(x)).loc[:,[i for i in MARKER_GENES if i in rna.columns]]
        rna_ax = ax[0,individuals.index(i)]
        rna_ax.set_title('{}, RNA'.format(i))
        sns.heatmap(rna, ax=rna_ax, cmap='Reds', cbar=None)
        rna_ax.set_ylim(top=rna_ax.get_ylim()[1]-0.5)
        rna_ax.set_ylim(bottom=rna_ax.get_ylim()[0]+1.5)

    atac = cpm[(cpm.individual==i) & (cpm.modality=='ATAC')]
    if len(atac) > 0:
        atac = atac.pivot(index='cluster', columns='feature', values='cpm').fillna(0)
        atac = atac.transform(lambda x: x/max(x)).loc[:,[i for i in MARKER_GENES if i in atac.columns]]
        atac_ax = ax[1,individuals.index(i)]
        atac_ax.set_title('{}, ATAC'.format(i))
        sns.heatmap(atac, ax=atac_ax, cmap='Reds', cbar=None)
        atac_ax.set_ylim(top=atac_ax.get_ylim()[1]-0.5)
        atac_ax.set_ylim(bottom=atac_ax.get_ylim()[0]+1.5)
fig.tight_layout()
fig.savefig(OUT)


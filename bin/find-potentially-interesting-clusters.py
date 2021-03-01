#!/usr/bin/env python
# coding: utf-8

# In[74]:


import os
import sys
import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#LIGER_IN_MATRIX_GLOB = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat/data/*'
#CLUSTERING_GLOB = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat/results/clusters/KSM2_RNA/*'
LIGER_IN_MATRIX_GLOB = sys.argv[1]
CLUSTERING_GLOB = sys.argv[2]

LIGER_IN_MATRICES = glob.glob(LIGER_IN_MATRIX_GLOB)
CLUSTERINGS = glob.glob(CLUSTERING_GLOB)

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


# In[ ]:


# given clusterings and marker genes -- score the specificity of each marker gene in each cluster


# In[8]:


count_matrices = dict()
for f in LIGER_IN_MATRICES:
    individual, modality = os.path.basename(f).split('.')[0].split('_')
    count_matrices[(individual, modality)] = pd.read_hdf(f)


# In[114]:


def load_clustering(f):
    _, k, resolution = os.path.basename(f).replace('.txt', '').split('-')
    tmp = pd.read_csv(f, header=None, sep=' ', names=['library', 'barcode', 'cluster'])
    tmp['k'] = k
    tmp['resolution'] = resolution
    return tmp

clusterings = pd.concat([load_clustering(f) for f in CLUSTERINGS])
clusterings['library'] = clusterings[['library', 'barcode']].apply(lambda x: x[0].replace(x[1], ''), axis=1)


# In[115]:


TARGET_NUMBER_CLUSTERS = [7, 8]
number_clusters = clusterings.groupby(['k', 'resolution']).cluster.max() + 1
number_clusters = number_clusters.reset_index()


# In[116]:


interesting_clusterings = number_clusters[number_clusters.cluster.isin(TARGET_NUMBER_CLUSTERS)]


# In[ ]:


all_tsi = []


# In[117]:


for k, resolution in zip(interesting_clusterings.k, interesting_clusterings.resolution):
    # for a given clustering, get the per-individual, per-modality
    #k = '20'
    #resolution = '0.05'
    clustering = clusterings[(clusterings.k==k) & (clusterings.resolution==resolution)]
    clustering['nucleus'] = clustering.library + '-' + clustering.barcode
    nucleus_to_cluster = dict(zip(clustering.nucleus, clustering.cluster))
    per_individual_modality_and_cluster_counts = dict()
    for key in count_matrices:
        per_individual_modality_and_cluster_counts[key] = count_matrices[key].groupby(nucleus_to_cluster).sum()
    
    # convert to CPM
    cpm = {key: val.transform(lambda x: 1e6*x/sum(x), axis=1) for key, val in per_individual_modality_and_cluster_counts.items()}
    
    # subset to marker genes
    MARKER_GENES = ['MYH1', 'MYH2', 'MYH4', 'MYH7', 'PDGFRA', 'VWF', 'ACTA2', 'CD163', 'PAX7']
    markers = {}
    for key in cpm:
        tmp = cpm[key]
        tmp.columns = [x.upper() for x in tmp.columns]
        markers[key] = tmp[MARKER_GENES]
    
    
    # make heatmap
    # top row: RNA
    # bottom row: ATAC
    # each column is an individual
    individuals = list(set([i[0] for i in markers]))
    fig, ax = plt.subplots(ncols=len(individuals), nrows=2, figsize=(3*len(individuals), 6))
    for i in individuals:
        rna = markers[(i, 'RNA')]
        # column normalize
        rna = rna.transform(lambda x: x/max(x))
        rna_ax = ax[0,individuals.index(i)]
        rna_ax.set_title('{}, RNA'.format(i))
        sns.heatmap(rna, ax=rna_ax, cmap='Reds', cbar=None)
        rna_ax.set_ylim(top=rna_ax.get_ylim()[1]-0.5)
        rna_ax.set_ylim(bottom=rna_ax.get_ylim()[0]+1.5)

        atac = markers[(i, 'ATAC')]
        # column normalize
        atac = atac.transform(lambda x: x/max(x))
        atac_ax = ax[1,individuals.index(i)]
        atac_ax.set_title('{}, ATAC'.format(i))
        sns.heatmap(atac, ax=atac_ax, cmap='Reds', cbar=None)
        atac_ax.set_ylim(top=atac_ax.get_ylim()[1]-0.5)
        atac_ax.set_ylim(bottom=atac_ax.get_ylim()[0]+1.5)

    fig.tight_layout()
    fig.savefig(f'heatmap-{k}-{resolution}.pdf')
    fig.clf()
    
    tsi = pd.concat([markers[key].apply(tissue_specificity_index).to_frame(name='tsi').assign(individual=key[0], modality=key[1]).reset_index() for key in markers])
    tsi = tsi.rename(columns={'index': 'gene'}).assign(k=k, resolution=resolution)
    all_tsi.append(tsi)


# In[ ]:

if len(all_tsi) > 0:
    all_tsi = pd.concat(all_tsi)
    all_tsi.to_csv('tsi.tsv', index=False, sep='\t')


#!/usr/bin/env python
# coding: utf-8

# In[35]:


import os
import sys
import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#CLUSTERING_GLOB = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-10/results/clusters/KSM2_RNA/*'
CLUSTERING_GLOB = sys.argv[1]
CLUSTERINGS = glob.glob(CLUSTERING_GLOB)
NUCLEUS_INDIVIDUAL_ASSIGNMENTS = '/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files/data/nuclei-with-individual-assignments.txt'
LIBRARY_LABELS = '/lab/work/porchard/sn-muscle-project/library-labels.txt'
TARGET_NUMBER_CLUSTERS = 7


# In[37]:


nucleus_to_individual = pd.read_csv(NUCLEUS_INDIVIDUAL_ASSIGNMENTS, sep='\t', header=None, names=['library', 'barcode', 'individual'])
nucleus_to_individual['nucleus'] = nucleus_to_individual.library + '-' + nucleus_to_individual.barcode
nucleus_to_individual = {nucleus: individual for nucleus, individual in zip(nucleus_to_individual.nucleus, nucleus_to_individual.individual)}


# In[38]:


def load_clustering(f):
    _, k, resolution = os.path.basename(f).replace('.txt', '').split('-')
    tmp = pd.read_csv(f, header=None, sep=' ', names=['library', 'barcode', 'cluster'])
    tmp['k'] = k
    tmp['resolution'] = resolution
    return tmp

clusterings = pd.concat([load_clustering(f) for f in CLUSTERINGS])
clusterings['library'] = clusterings[['library', 'barcode']].apply(lambda x: x[0].replace(x[1], ''), axis=1)
clusterings['nucleus'] = clusterings.library + '-' + clusterings.barcode
clusterings['individual'] = clusterings.nucleus.map(nucleus_to_individual)
clusterings['id'] = clusterings['k'] + '___' + clusterings.resolution


# In[39]:





# In[40]:


library_labels = pd.read_csv(LIBRARY_LABELS, sep='\t')
library_to_modality = dict(zip(library_labels.library, library_labels.modality))

MISSING_MODALITY = [i for i in clusterings.library.unique() if i not in library_to_modality]
for i in MISSING_MODALITY:
    if 'RNA' in i:
        library_to_modality[i] = 'RNA'
    else:
        library_to_modality[i] = 'ATAC'


# In[42]:


number_clusters = (clusterings.groupby('id').cluster.max() + 1).reset_index()


# In[43]:


clusterings['modality'] = clusterings.library.map(library_to_modality)


# In[58]:


## look for clusterings with:
# target number of clusters
keep_clusterings = number_clusters[number_clusters.cluster==TARGET_NUMBER_CLUSTERS].id.to_list()
interesting_clusterings = clusterings[clusterings.id.isin(keep_clusterings)]

# each modality and individual in each cluster
N_LIBRARY_MODALITY_COMBOS = len(clusterings[['individual', 'modality']].drop_duplicates())
per_id_per_cluster_count = interesting_clusterings[['id', 'cluster', 'individual', 'modality']].drop_duplicates().groupby(['id', 'cluster']).size()
keep_clusterings = per_id_per_cluster_count.reset_index().groupby('id')[0].min().where(lambda x: x==N_LIBRARY_MODALITY_COMBOS).dropna().index.to_list()
interesting_clusterings = interesting_clusterings[interesting_clusterings.id.isin(keep_clusterings)]
# two large clusters > 20% of data each, rest small clusters (<15% of data)

interesting_clusterings.head()


# In[ ]:


interesting_clusterings.to_csv(sys.stdout, sep='\t', index=False)

for i in interesting_clusterings.id.unique():
    FILENAME = f'interesting-clustering.{i}.txt'
    interesting_clusterings.loc[interesting_clusterings.id==i,['library', 'barcode', 'cluster']].to_csv(FILENAME, sep='\t', index=False, header=False)

# In[20]:


#counts = interesting_clusterings.groupby(['id', 'cluster']).size().reset_index().rename(columns={0: 'n'})
#counts['fraction'] = counts.groupby('id').n.transform(lambda x: x/sum(x))
#counts.head()


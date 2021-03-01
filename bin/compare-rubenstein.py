#!/usr/bin/env python
# coding: utf-8

# In[52]:


import os
import sys
import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

#LIGER_IN_MATRIX_GLOB = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/results/long-matrices/*RNA*'
#CLUSTERING = '/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/results/find-promising/KSM2_RNA/interesting-clustering.15___0.05.txt'
#LIGER_IN_MATRICES = glob.glob(LIGER_IN_MATRIX_GLOB)
INCLUDE_LIBRARIES = ['63_20_rna-hg19', '63_40_rna-hg19', '133155-hg19', '133156-hg19', '133157-hg19', '133158-hg19']
#RUBENSTEIN = '/lab/work/porchard/sn-muscle-project/data/rubenstein-fiber-type-differential-genes/fiber-type-differential-genes.txt'

CLUSTERING = sys.argv[1]
RUBENSTEIN = sys.argv[2]
OUT = sys.argv[3]
LIGER_IN_MATRICES = sys.argv[4:]

def nucleus_to_library(n):
    return '-'.join(n.split('-')[:-1])


# In[32]:


rubenstein = pd.read_csv(RUBENSTEIN, sep='\t')
rubenstein = rubenstein[rubenstein.padj<=0.05]
rubenstein.loc[rubenstein.fiber_type=='Type I','log2FC'] = -1*rubenstein.loc[rubenstein.fiber_type=='Type I','log2FC']


# In[33]:


def feature_count_file_to_cluster_counts(f, nucleus_to_cluster):
    counts = {}
    with open(f, 'r') as fh:
        for line in fh:
            nucleus, feature, count = line.rstrip().split('\t')
            if nucleus_to_library(nucleus) not in INCLUDE_LIBRARIES:
                continue
            count = int(count)
            cluster = nucleus_to_cluster[nucleus]
            if cluster not in counts:
                counts[cluster] = dict()
            if feature not in counts[cluster]:
                counts[cluster][feature] = 0
            counts[cluster][feature] += count
    return counts

clustering = pd.read_csv(CLUSTERING, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clustering['nucleus'] = clustering.library + '-' + clustering.barcode
clustering.cluster = clustering.cluster.astype(int)
nucleus_to_cluster = dict(zip(clustering.nucleus, clustering.cluster))


# In[8]:


count_matrices = dict()
for f in LIGER_IN_MATRICES:
    sys.stderr.write(f'{f}\n')
    individual, modality = os.path.basename(f).split('.')[0].split('_')
    if modality != 'RNA':
        continue
    count_matrices[(individual, modality)] = feature_count_file_to_cluster_counts(f, nucleus_to_cluster)


counts = []
for key in count_matrices:
    individual, modality = key
    for cluster in count_matrices[key]:
        for gene in count_matrices[key][cluster]:
            counts.append([individual, modality, cluster, gene, count_matrices[key][cluster][gene]])
counts = pd.DataFrame(counts, columns=['individual', 'modality', 'cluster', 'feature', 'count'])


# In[71]:


#cpm = counts[(counts.individual=='KSM1') & (counts.cluster.isin([0, 1]))].pivot(index='feature', columns='cluster', values='count').fillna(0)
cpm = counts[counts.cluster.isin([0, 1])].groupby(['cluster', 'feature'])['count'].sum().reset_index().pivot(index='feature', columns='cluster', values='count').fillna(0)
cpm = cpm.transform(lambda x: 1e6*x/sum(x))
cpm = cpm.loc[rubenstein.gene,:]
cpm['our_log2FC'] = np.log2(cpm[0]/cpm[1])
both = cpm[['our_log2FC']].join(rubenstein[['gene', 'log2FC']].set_index('gene'))


# In[77]:


fig, ax = plt.subplots(figsize=(9,9))
ax.scatter(x=both.log2FC, y=both.our_log2FC)
ax.axhline(0, color='grey', linestyle='dashed')
ax.axvline(0, color='grey', linestyle='dashed')
ax.set_xlabel('Rubenstein et al log2(Type IIA / Type I fibers)', size=15)
ax.set_ylabel('Our log2(Type II / Type I fibers)', size=15)
MIN = min(ax.get_ylim()[0], ax.get_xlim()[0])
MAX = max(ax.get_ylim()[1], ax.get_xlim()[1])
ax.set_xlim(left=MIN, right=MAX)
ax.set_ylim(bottom=MIN, top=MAX)
ax.plot([MIN, MAX], [MIN, MAX], linestyle='dotted', color='grey')
text = [ax.text(x=r['log2FC'], y=r['our_log2FC'], s=n) for n, r in both.iterrows()]
adjust_text(text, arrowprops={'color': 'red', 'arrowstyle':'-'})
fig.tight_layout()
fig.savefig(OUT)
#'rubenstein-vs-our-fiber-type-lfcs.pdf')


# In[ ]:





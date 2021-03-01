#!/usr/bin/env python
# coding: utf-8

# In[68]:


import os
import sys
import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import cmap_utils

#CLUSTERING = '/lab/work/porchard/sn-muscle-project/work/process-by-cluster-promising-2021-01-12/clusters.txt'
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
#LIBRARY_LABELS = '/lab/work/porchard/sn-muscle-project/library-labels.txt'

CLUSTERING = sys.argv[1]
CLUSTER_NAMES = sys.argv[2]
LIBRARY_LABELS = sys.argv[3]


# In[60]:


def is_dual_modality_atac(library):
    return ('1846' in library and 'ATAC' in library)


# In[69]:


library_labels = pd.read_csv(LIBRARY_LABELS, sep='\t')
library_labels.head()


# In[61]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))


# In[62]:


clusters = pd.read_csv(CLUSTERING, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clusters['library'] = clusters[['library', 'barcode']].apply(lambda x: x[0].replace(x[1], ''), axis=1)
clusters['nucleus'] = clusters.library + '-' + clusters.barcode
nucleus_to_cluster = dict(zip(clusters.nucleus, clusters.cluster))


# In[74]:


library_labels = pd.read_csv(LIBRARY_LABELS, sep='\t')
human_atac_libraries = [i for i in clusters.library.unique() if 'hg19' in i and 'ATAC' in i] + library_labels[(library_labels.modality=='ATAC') & (library_labels.species=='Human')].library.to_list()
rat_atac_libraries = [i for i in clusters.library.unique() if 'rn6' in i and 'ATAC' in i] + library_labels[(library_labels.modality=='ATAC') & (library_labels.species=='Rat')].library.to_list()


# In[65]:


# Fig 2E
total_nucleus_counts = clusters[~clusters.library.map(is_dual_modality_atac)].groupby('cluster').size().reset_index().rename(columns={0: 'nuclei'})
total_nucleus_counts['frac'] = total_nucleus_counts.nuclei / total_nucleus_counts.nuclei.sum()
total_nucleus_counts['cluster_name'] = total_nucleus_counts.cluster.map(cluster_to_cluster_name)


# In[66]:


cluster_name_to_color = cmap_utils.make_colormap_dict(total_nucleus_counts.cluster_name.to_list())


# In[51]:


fig, ax = plt.subplots(figsize=(4,3))
sns.barplot(y='cluster_name', x='frac', data=total_nucleus_counts, palette=cluster_name_to_color, ax=ax)
ax.set_xlabel('Fraction of nuclei')
ax.set_ylabel('')
for n, r in total_nucleus_counts.iterrows():
    counts = f'{r["nuclei"]:,}'
    if r['cluster'] == 0:
        label = '{} ({}%)'.format(counts, round(100*r['frac'],1))
        color = 'white'
        ax.text(x=0.05, y=r['cluster'], s=label, color=color, ha='left', va='center')
    elif r['cluster'] == 1:
        label_1 = '{}'.format(counts)
        label_2 = '({}%)'.format(round(100*r['frac'],1))
        color_1 = 'white'
        color_2 = 'black'
        ax.text(x=0.05, y=r['cluster'], s=label_1, color=color_1, ha='left', va='center')
        ax.text(x=r['frac'], y=r['cluster'], s=label_2, color=color_2, ha='left', va='center')
    else:
        label = '{} ({}%)'.format(counts, round(100*r['frac'],1))
        color = 'black'
        ax.text(x=0.05, y=r['cluster'], s=label, color=color, ha='left', va='center')
fig.tight_layout()
fig.savefig('all-nuclei-cell-type-proportions.png')
fig.clf()


# In[121]:


# Fig 2H
fig, axs = plt.subplots(nrows=2, figsize=(3.5,4))

# plot for human
human = clusters[clusters.library.isin(human_atac_libraries)].groupby('cluster').size().reset_index().rename(columns={0: 'nuclei'})
human['cluster_name'] = human.cluster.map(cluster_to_cluster_name)
human_ax = axs[0]
sns.barplot(x='nuclei', y='cluster_name', data=human, palette=cluster_name_to_color, ax=human_ax)
human_ax.set_ylabel('')
human_ax.set_xlabel('')
human_ax.yaxis.set_label_position('right')
human_ax.yaxis.tick_right()
for n, r in human.iterrows():
    counts = f'{r["nuclei"]:,}'
    color = 'white' if r['cluster'] == 0 else 'black'
    align =  'right' if r['cluster'] == 0 else 'left'
    human_ax.text(x=r['nuclei'], y=r['cluster'], s=counts, color=color, ha=align, va='center')

# plot for rat
rat = clusters[clusters.library.isin(rat_atac_libraries)].groupby('cluster').size().reset_index().rename(columns={0: 'nuclei'})
rat['cluster_name'] = rat.cluster.map(cluster_to_cluster_name)
rat_ax = axs[1]
sns.barplot(x='nuclei', y='cluster_name', data=rat, palette=cluster_name_to_color, ax=rat_ax)
rat_ax.set_ylabel('')
rat_ax.set_xlabel('# nuclei with\nchromatin accessibility')
rat_ax.yaxis.set_label_position('right')
rat_ax.yaxis.tick_right()
for n, r in rat.iterrows():
    counts = f'{r["nuclei"]:,}'
    color = 'white' if r['cluster'] == 0 else 'black'
    align =  'right' if r['cluster'] == 0 else 'left'
    rat_ax.text(x=r['nuclei'], y=r['cluster'], s=counts, color=color, ha=align, va='center')

fig.tight_layout()
fig.savefig('snATAC-nuclei-cell-type-proportions-by-species.png')
fig.clf()


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[13]:


import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


# In[4]:


#MODEL_RESULTS_HUMAN = '/lab/work/porchard/sn-muscle-project/work/logistic-regression-multiple-states-hg19/results/enhancer-regression/model_results.txt'
#MODEL_RESULTS_RAT = '/lab/work/porchard/sn-muscle-project/work/logistic-regression-multiple-states-rn6/results/enhancer-regression/model_results.txt'
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/work/manuscript-figures-v2/work/59/3d57abecd67d873791b89ec0b58bc9/cluster-names.txt'
#CELL_TYPE_DECODING = '/lab/work/porchard/sn-muscle-project/data/roadmap-posteriors/roadmap_cell_types.txt'

MODEL_RESULTS_HUMAN, MODEL_RESULTS_RAT, CLUSTER_NAMES, CELL_TYPE_DECODING = sys.argv[1:]


# In[70]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_to_cluster_name = {'cluster_{}'.format(i): j for i, j in zip(cluster_names.old_name, cluster_names.new_name)}


# In[6]:


human = pd.read_csv(MODEL_RESULTS_HUMAN, sep='\t').assign(species='human')
rat = pd.read_csv(MODEL_RESULTS_RAT, sep='\t').assign(species='rat')
models = pd.concat([human, rat])


# In[10]:


# normalize coefficient
models.coef = models.coef.map(lambda x: 0 if x < 0 else x)
models.coef = models.groupby(['cluster' , 'species']).coef.transform(lambda x: x/max(x))


# In[8]:


tissue_mappings = pd.read_csv(CELL_TYPE_DECODING, sep='\t')[['color', 'eid', 'standardized_name']]


# In[18]:


cell_type_to_label = {eid: '{} ({})'.format(standardized_name, eid) for eid, standardized_name in zip(tissue_mappings.eid, tissue_mappings.standardized_name)}


# In[24]:


# plot out the full cluster map
# human first
tmp = models.loc[models.species=='human',['cluster', 'cell_type', 'coef']].pivot(index='cluster', columns='cell_type', values='coef')
tmp.columns = tmp.columns.to_series().map(cell_type_to_label)
g = sns.clustermap(tmp)
fig = px.imshow(g.data2d, color_continuous_midpoint=0.5, color_continuous_scale='Reds', origin='lower', labels={'color': 'Enhancer score'}, aspect='auto')
fig.update_layout(title='', template='plotly_white', yaxis_title='', width=2000, height=700)
fig.write_html('human.html')

# rat second
tmp = models.loc[models.species=='rat',['cluster', 'cell_type', 'coef']].pivot(index='cluster', columns='cell_type', values='coef')
tmp.columns = tmp.columns.to_series().map(cell_type_to_label)
g = sns.clustermap(tmp)
fig = px.imshow(g.data2d, color_continuous_midpoint=0.5, color_continuous_scale='Reds', origin='lower', labels={'color': 'Enhancer score'}, aspect='auto')
fig.update_layout(title='', template='plotly_white', yaxis_title='', width=2000, height=700)
fig.write_html('rat.html')


# In[82]:


# plot the figure
KEEP = ['E100', 'E023', 'E122', 'E111', 'E029', 'E089']
LABELS = {'Psoas muscle' : 'E100','MSC derived adipocytes' : 'E023', 'HUVEC cells (Endothelial)' : 'E122', 'Smooth muscle' : 'E111', 'Monocytes' : 'E029', 'Fetal trunk muscle' : 'E089'}
LABELS = {val: key for key, val in LABELS.items()}

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(3,4), gridspec_kw={'height_ratios': [0.1, 1, 1], 'hspace': 0.08})

cbar_ax = axs[0]

# human first
ax = axs[1]
tmp = models.loc[(models.species=='human') & (models.cell_type.isin(KEEP)),['cluster', 'cell_type', 'coef']].pivot(index='cluster', columns='cell_type', values='coef')
tmp = tmp[KEEP].sort_index()
fcbar = ax.imshow(tmp, cmap='bwr', vmin=0, vmax=1, aspect='auto')
ax.set_xlabel('')
ax.set_xticks([])
ax.set_yticks(range(len(tmp)))
ax.set_yticklabels([cluster_to_cluster_name[i] for i in tmp.index])
ax.set_ylim(bottom=ax.get_ylim()[0]+0.5, top=ax.get_ylim()[1]-0.5)

plt.colorbar(fcbar, cax=cbar_ax, orientation='horizontal')
cbar_ax.set_xlabel('Enhancer similarity score')
cbar_ax.xaxis.tick_top()
cbar_ax.xaxis.set_label_position('top') 

# rat
ax = axs[2]
tmp = models.loc[(models.species=='rat') & (models.cell_type.isin(KEEP)),['cluster', 'cell_type', 'coef']].pivot(index='cluster', columns='cell_type', values='coef')
tmp = tmp[KEEP].sort_index()
ax.imshow(tmp, cmap='bwr', vmin=0, vmax=1, aspect='auto')
ax.set_xticks(range(len(LABELS)))
ax.set_xticklabels([LABELS[i] for i in KEEP], rotation=30, ha='right')
ax.set_yticks(range(len(tmp)))
ax.set_yticklabels([cluster_to_cluster_name[i] for i in tmp.index])
ax.set_ylim(bottom=ax.get_ylim()[0]+0.5, top=ax.get_ylim()[1]-0.5)
ax.set_xlabel('Roadmap cell type')
fig.savefig('enhancer-similarity-heatmap.pdf', bbox_inches='tight')
fig.clf()


#!/usr/bin/env python
# coding: utf-8

# In[209]:


import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import re
import glob
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

SOURCE_SN_ATAC_OURS = ['Type I muscle fibers', 'Type II muscle fibers', 'Mesenchymal stem cells', 'Endothelial cells', 'Smooth muscle', 'Immune cells', 'Muscle satellite cells']
SOURCE_SN_ATAC_OTHER = ['Beta cells']
SOURCE_BULK_ATAC = ['Adipose']


# In[210]:


#LDSC_RESULT_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/ldsc/UKB/hg19/joint-model/results/partitioned-heritability/*.results')
#PHENOTYPE_TSV = '/lab/work/porchard/sn-muscle-project/data/ukb-summary-stats/ukb31063_h2_all.02Oct2019.tsv.gz'
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
CLUSTER_NAMES = sys.argv[1]
PHENOTYPE_TSV = sys.argv[2]
LDSC_RESULT_FILES = sys.argv[3:]


# In[211]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_id_to_cluster_name = {f'cluster_{i}': j for i, j in zip(cluster_names.old_name, cluster_names.new_name)}


# In[212]:


phenotypes = pd.read_csv(PHENOTYPE_TSV, sep='\t')[['phenotype', 'description']].astype(str)
UPDATE_TRAIT_DESCRIPTIONS = {'Diseases of veins, lymphatic vessels and lymph nodes, not elsewhere classified': 'Diseases of veins, lymphatic vessels and lymph nodes', 'Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor: Hayfever, allergic rhinitis or eczema': 'Diagnosed by doctor: Hayfever, allergic rhinitis or eczema'}               
phenotypes.description = phenotypes.description.map(lambda x: UPDATE_TRAIT_DESCRIPTIONS[x] if x in UPDATE_TRAIT_DESCRIPTIONS else x)
phenotype_to_description = dict(zip(phenotypes.phenotype, phenotypes.description))


# In[213]:


def load_ldsc_result_file(f):
    tmp = pd.read_csv(f, sep='\t')
    tmp = tmp[tmp.Category.map(lambda x: 'L2_1' in x)]
    tmp.Category = tmp.Category.map(lambda x: x.replace('L2_1', '')).map(lambda x: 'H9 hESC' if x == 'Germ' else x)
    tmp['phenotype'] = os.path.basename(f).split('.')[1]
    return tmp


# In[214]:


ldsc = pd.concat([load_ldsc_result_file(f) for f in LDSC_RESULT_FILES])
ldsc['Category'] = ldsc.Category.map(lambda x: cluster_id_to_cluster_name[x] if x in cluster_id_to_cluster_name else x)
ldsc['one_tailed_pvalue'] = ldsc['Coefficient_z-score'].map(lambda x: 1-norm.cdf(x))
ldsc['two_tailed_pvalue'] = ldsc['Coefficient_z-score'].map(lambda x: 2*(1-norm.cdf(abs(x))))
ldsc['pvalue'] = ldsc['two_tailed_pvalue']


# In[215]:


# Examine tissue of top coefficient z-score per trait
top_tissue_per_trait = ldsc.sort_values('Coefficient_z-score', ascending=False).groupby('phenotype').head(1).assign(description=lambda df: df.phenotype.map(phenotype_to_description))
top_tissue_per_trait.to_csv('top-z-score-tissue-per-trait.tsv', sep='\t', index=False)


# Examine tissue of top enrichment per trait
top_tissue_per_trait = ldsc.sort_values('Enrichment', ascending=False).groupby('phenotype').head(1).assign(description=lambda df: df.phenotype.map(phenotype_to_description))
top_tissue_per_trait.to_csv('top-enrichment-tissue-per-trait.tsv', sep='\t', index=False)

# plot all things that are sign. in one of our cell types after BY correction.
our_cell_types = ldsc[ldsc.Category.isin(cluster_id_to_cluster_name.values())]
our_cell_types['BY'] = multipletests(our_cell_types.pvalue, method='fdr_by')[1]
SHOW_PHENOTYPES = our_cell_types[our_cell_types.BY<=0.05].phenotype

# to plot:
# pass the matrix to plot
# also, pass the elements to asterisk
def make_heatmap(m, asterisk, cluster_rows=True):
    g = sns.clustermap(m, cmap='vlag', figsize=(10, 15))
    df = g.data2d.T.loc[m.columns.to_list(),:] if not cluster_rows else g.data2d.T
    LIMITS = df.abs().max().max()
    fig, ax = plt.subplots(figsize=(len(m)/4.5, 3+len(m.columns)/5))
    fcbar = ax.imshow(df, aspect='auto', cmap='vlag', vmin=-1*LIMITS, vmax=LIMITS)
    ax.set_yticks(range(len(df)))
    ax.set_xticks(range(len(df.columns)))
    ax.set_ylim(bottom=ax.get_ylim()[0]+0.5, top=ax.get_ylim()[1]-0.5)
    ax.set_yticklabels(df.index.to_list())
    ax.set_xticklabels(df.columns.to_list(), rotation=90)
    cbar = plt.colorbar(fcbar)
    cbar.set_label('LDSC Z-score')
    for (category, phenotype) in asterisk:
        col = df.columns.to_list().index(phenotype)
        row = df.index.to_list().index(category)
        ax.text(col, row, s='*', va='center', ha='center')
    fig.tight_layout()
    return (fig, ax)



# plot all things that are sign. in one of our cell types after BY correction.
# correct within our cell types only, but can show data for all other cell types too
our_cell_types = ldsc[ldsc.Category.isin(cluster_id_to_cluster_name.values())]
our_cell_types['BY'] = multipletests(our_cell_types.pvalue, method='fdr_by')[1]
SHOW_PHENOTYPES = our_cell_types[our_cell_types.BY<=0.05].phenotype

SIGNIFICANT_SUBSET = our_cell_types[(our_cell_types.BY<=0.05) & (our_cell_types['Coefficient_z-score']>0)]
SIGNIFICANT_PAIRS = SIGNIFICANT_SUBSET[['Category', 'phenotype']]
SIGNIFICANT_PAIRS.phenotype = SIGNIFICANT_PAIRS.phenotype.map(lambda x: '{} ({})'.format(phenotype_to_description[x], x))
SHOW_PHENOTYPES = SIGNIFICANT_SUBSET.phenotype.unique()
df = ldsc[(ldsc.Category.isin(cluster_id_to_cluster_name.values())) & (ldsc.phenotype.isin(SHOW_PHENOTYPES))].loc[:,['Category', 'phenotype', 'Coefficient_z-score']].pivot(index='phenotype', columns='Category', values='Coefficient_z-score')
df.index = df.index.to_series().map(lambda x: '{} ({})'.format(phenotype_to_description[x], x))
df = df.loc[:,cluster_names.sort_values('old_name').new_name.to_list()] # re-order cell types
asterisk = [[r['Category'], r['phenotype']] for n, r in SIGNIFICANT_PAIRS.iterrows()]
fig, ax = make_heatmap(df, asterisk)
fig.tight_layout()
fig.savefig('LDSC-UKB-across-muscle-cell-types.pdf')
fig.clf()


# BY correction across ALL cell types
all_cell_types = ldsc.copy()
all_cell_types['BY'] = multipletests(all_cell_types.pvalue, method='fdr_by')[1]
all_cell_types.Category = all_cell_types.Category.map(lambda x:  'Beta cells' if x == 'beta_ATAC' else x).map(lambda x: x.replace('_', ' ').capitalize().replace(' i ', ' I ').replace(' ii ', ' II '))

SIGNIFICANT_SUBSET = all_cell_types[(all_cell_types.BY<=0.05) & (all_cell_types['Coefficient_z-score']>0)]
SIGNIFICANT_PAIRS = SIGNIFICANT_SUBSET[['Category', 'phenotype']]
SIGNIFICANT_PAIRS.phenotype = SIGNIFICANT_PAIRS.phenotype.map(lambda x: '{} ({})'.format(phenotype_to_description[x], x))
SHOW_PHENOTYPES = SIGNIFICANT_SUBSET.phenotype.unique()
df = all_cell_types[all_cell_types.phenotype.isin(SHOW_PHENOTYPES)].loc[:,['Category', 'phenotype', 'Coefficient_z-score']].pivot(index='phenotype', columns='Category', values='Coefficient_z-score')
df.index = df.index.to_series().map(lambda x: '{} ({})'.format(phenotype_to_description[x], x))
df = df.loc[:,df.columns.to_series().sort_values()]
asterisk = [[r['Category'], r['phenotype']] for n, r in SIGNIFICANT_PAIRS.iterrows()]
fig, ax = make_heatmap(df, asterisk)
for i in ax.get_yticklabels():
    if i.get_text() in SOURCE_SN_ATAC_OURS:
        i.set(color='#a6611a')
    elif i.get_text() in SOURCE_SN_ATAC_OTHER:
        i.set(color='#dfc27d')
    elif i.get_text() in SOURCE_BULK_ATAC:
        i.set(color='#018571')
    else:
        i.set(color='#80cdc1')
fig.tight_layout()
fig.savefig('LDSC-UKB-across-all-cell-types.pdf')
fig.clf()


# In[ ]:





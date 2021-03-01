#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import glob


# In[90]:


CLUSTERS = '/lab/work/porchard/sn-muscle-project/work/process-by-cluster/clusters.txt'
CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
INDIVIDUAL_ASSIGNMENTS = '/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files/data/nuclei-with-individual-assignments.txt'
ATAC_METRICS = ['/lab/work/porchard/sn-muscle-project/work/qc/results/call-nuclei-atac/metrics.txt'] + glob.glob('/lab/work/porchard/Nova-303/work/qc/results/ataqv/*')
RNA_METRICS = glob.glob('/lab/work/porchard/sn-muscle-project/work/rnaseq-qc/results/rnaseq-qc/*') + glob.glob('/lab/work/porchard/Nova-315/work/starsolo-qc/results/qc/*')


# In[48]:


clusters = pd.read_csv(CLUSTERS, header=None, sep='\t', names=['library', 'barcode', 'cluster'])
clusters['nucleus'] = clusters.library + '-' + clusters.barcode


# In[49]:


individuals = pd.read_csv(INDIVIDUAL_ASSIGNMENTS, sep='\t', header=None, names=['library', 'barcode', 'individual'])
# drop ATAC FANS libraries
individuals = individuals[~individuals.library.isin(['133152-hg19', '133154-hg19'])]
individuals['nucleus'] = individuals.library + '-' + individuals.barcode
individuals.head()


# In[50]:


individuals = individuals.merge(clusters, on=['library', 'barcode', 'nucleus'], how='left')


# In[53]:


# remap library names

LIBRARY_FOR_TABLE = {'125589': 'Human-rat mix (ATAC)',
    '133151': 'no FANS rep1 (ATAC)',
    '133153': 'no FANS rep2 (ATAC)',
    '133155': 'no FANS rep1 (RNA)',
    '133156': 'FANS rep1 (RNA)',
    '133157': 'no FANS rep2 (RNA)',
    '133158': 'FANS rep2 (RNA)',
    '63_20': '20k nuclei (ATAC)',
    '63_40': '40k nuclei (ATAC)',
    '63_20_rna': '20k nuclei (RNA)',
    '63_40_rna': '40k nuclei (RNA)',
    '1846_ATAC': 'Human-rat mix (MULTIOME)',
    '1846_RNA': 'Human-rat mix (MULTIOME)'}

individuals['new_library'] = individuals.library.map(lambda x: LIBRARY_FOR_TABLE[x.split('-')[0]])
individuals.head()


# In[54]:


# get fragment counts for ATAC
atac_metrics = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['nucleus', 'metric', 'value']) for f in ATAC_METRICS])
atac_metrics = atac_metrics[atac_metrics.nucleus.isin(individuals.nucleus)]
atac_metrics = atac_metrics[atac_metrics.metric=='hqaa']
atac_metrics = atac_metrics[['nucleus', 'value']].rename(columns={'value': 'atac_reads'})
atac_metrics['atac_fragments'] = atac_metrics.atac_reads.astype(int).map(lambda x: round(x/2))
atac_metrics = atac_metrics[['nucleus', 'atac_fragments']]
atac_metrics.head()


# In[55]:


# get UMI counts for RNA
rna_metrics = pd.concat([pd.read_csv(f, sep='\t').assign(library=os.path.basename(f).split('.')[0]) for f in RNA_METRICS])
rna_metrics['nucleus'] = rna_metrics.library + '-' + rna_metrics.barcode
rna_metrics = rna_metrics[rna_metrics.nucleus.isin(individuals.nucleus)]
rna_metrics = rna_metrics[['nucleus', 'umis']].rename(columns = {'umis': 'rna_umis'})
rna_metrics.rna_umis = rna_metrics.rna_umis.astype(int)


# In[56]:


individuals = individuals.merge(rna_metrics, on='nucleus', how='left').merge(atac_metrics, on='nucleus', how='left')


# In[73]:


#individuals.groupby('new_library').loc[:,['rna_umis', 'atac_fragments']].agg([np.mean, np.median])
median_umis = individuals[~individuals.rna_umis.isnull()].groupby('new_library').rna_umis.median().map(round).reset_index().rename(columns={'rna_umis': 'Median RNA UMIs'})
mean_umis = individuals[~individuals.rna_umis.isnull()].groupby('new_library').rna_umis.mean().map(round).reset_index().rename(columns={'rna_umis': 'Mean RNA UMIs'})
median_atac = individuals[~individuals.atac_fragments.isnull()].groupby('new_library').atac_fragments.median().map(round).reset_index().rename(columns={'atac_fragments': 'Median ATAC fragments'})
mean_atac = individuals[~individuals.atac_fragments.isnull()].groupby('new_library').atac_fragments.mean().map(round).reset_index().rename(columns={'atac_fragments': 'Mean ATAC fragments'})


# In[82]:


# make table S3
# for counts per individual, don't double-count multi-ome nuclei
counts_per_sample = individuals[~individuals.library.map(lambda x: '1846' in x and 'ATAC' in x)].groupby(['new_library', 'individual']).size().reset_index(name='nuclei').pivot(index='new_library', values='nuclei', columns='individual').fillna(0).astype(int)
counts_per_sample.columns = counts_per_sample.columns.map(lambda x: 'nuclei from {}'.format(x.replace('rat1', 'rat').replace('KSM', 'HSM')))
counts_per_sample = counts_per_sample.join(median_umis.set_index('new_library'))
counts_per_sample = counts_per_sample.join(mean_umis.set_index('new_library'))
counts_per_sample = counts_per_sample.join(median_atac.set_index('new_library'))
counts_per_sample = counts_per_sample.join(mean_atac.set_index('new_library'))
table_s3 = counts_per_sample.fillna(0).astype(int).reset_index().rename(columns={'new_library': 'library'})
table_s3.to_csv('per-library-summary-stats.tsv', sep='\t', index=False)


# In[100]:


# make table S5
table_s5 = individuals[~individuals.library.map(lambda x: '1846' in x and 'ATAC' in x)]
table_s5['modality'] = table_s5.new_library.map(lambda x: re.match('.*\((.*)\)', x).group(1))
table_s5['species'] = table_s5.library.map(lambda x: x.split('-')[-1]).map({'hg19': 'human', 'rn6': 'rat'})
table_s5 = table_s5[~table_s5.cluster.isnull()]
table_s5.cluster = table_s5.cluster.astype(int)
table_s5 = table_s5.groupby(['modality', 'species', 'cluster']).size().reset_index(name='n')
cluster_to_cluster_name = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_to_cluster_name = {int(old_name): new_name for old_name, new_name in zip(cluster_to_cluster_name.old_name, cluster_to_cluster_name.new_name)}
table_s5['col'] = table_s5.species + ' ' + table_s5.modality
table_s5 = table_s5[['cluster', 'col', 'n']].pivot(index='cluster', columns='col', values='n').reset_index()
table_s5.cluster = table_s5.cluster.map(cluster_to_cluster_name)
table_s5.to_csv('per-species-per-modality-nucleus-counts.tsv', sep='\t', index=False)


# In[ ]:





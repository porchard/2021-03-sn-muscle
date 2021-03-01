#!/usr/bin/env python
# coding: utf-8

# In[75]:


import os
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import matplotlib.ticker as ticker
import argparse

# RNA
#DEMUXLET_UNMASKED = ['/lab/work/porchard/sn-muscle-project/work/qc/results/demuxlet/unmasked/63_20_rna-hg19.best.txt', '/lab/work/porchard/sn-muscle-project/work/qc/results/demuxlet/unmasked/63_40_rna-hg19.best.txt']
#DEMUXLET_MASKED = ['/lab/work/porchard/sn-muscle-project/work/qc/results/demuxlet/masked/63_20_rna-hg19.best.txt', '/lab/work/porchard/sn-muscle-project/work/qc/results/demuxlet/masked/63_40_rna-hg19.best.txt']
#THRESHOLDS = '/lab/work/porchard/sn-muscle-project/initial-thresholds-rna.txt'
#METRICS = glob.glob('/lab/work/porchard/sn-muscle-project/work/rnaseq-qc/results/rnaseq-qc/*')
#LIBRARY_LABELS = '/lab/work/porchard/sn-muscle-project/library-labels.txt'

LIBRARY_TO_INDIVIDUALS = {
    '133155-hg19': ['KSM1'],
    '133156-hg19': ['KSM1'],
    '133157-hg19': ['KSM1'],
    '133158-hg19': ['KSM1'],
    '63_20_rna-hg19': ['KSM1', 'KSM2'],
    '63_40_rna-hg19': ['KSM1', 'KSM2']
}


# In[ ]:


parser = argparse.ArgumentParser('Enforce QC thresholds to call RNA nuclei.')
parser.add_argument('--demuxlet-masked', dest='demuxlet_masked', nargs='+', help='Demuxlet masked best files')
parser.add_argument('--demuxlet-unmasked', dest='demuxlet_unmasked', nargs='+', help='Demuxlet masked best files')
parser.add_argument('--metrics', nargs='+', help='RNA-seq metric files')
parser.add_argument('--thresholds', help='Fle of initial thresholds')
#parser.add_argument('--library-labels', dest='library_labels', help='Ataqv metric file')
args = parser.parse_args()

DEMUXLET_MASKED = args.demuxlet_masked
DEMUXLET_UNMASKED = args.demuxlet_unmasked
THRESHOLDS = args.thresholds
METRICS = args.metrics
#LIBRARY_LABELS = args.library_labels


# In[76]:


thresholds = pd.read_csv(THRESHOLDS, sep='\t')


# In[77]:


def load_metrics(f):
    tmp = pd.read_csv(f, sep='\t')
    tmp['library'] = os.path.basename(f).replace('.rnaseq-qc.txt', '')
    tmp['nucleus'] = tmp.library + '-' + tmp.barcode
    return tmp


metrics = pd.concat([load_metrics(f) for f in METRICS])


# In[79]:


def load_demuxlet(f):
    tmp = pd.read_csv(f, sep='\t').loc[:,['BARCODE', 'BEST', 'SNG.1ST']].rename(columns={'BARCODE': 'barcode', 'BEST': 'best', 'SNG.1ST': 'sng'})
    tmp['assignment'] = tmp.best.map(lambda x: x.split('-')[0])
    tmp['library'] = os.path.basename(f).replace('.best.txt', '')
    assert(all(tmp.assignment.isin(['SNG', 'DBL', 'AMB'])))
    return tmp


# In[80]:


demuxlet_masked = pd.concat([load_demuxlet(f).assign(masked='masked') for f in DEMUXLET_MASKED])
demuxlet_unmasked = pd.concat([load_demuxlet(f).assign(masked='unmasked') for f in DEMUXLET_UNMASKED])
demuxlet = pd.concat([demuxlet_masked, demuxlet_unmasked])
demuxlet['nucleus'] = demuxlet.library + '-' + demuxlet.barcode
demuxlet_tested = demuxlet.nucleus.unique() # record all that were tested in demuxlet
demuxlet.head()


# In[93]:


number_individuals_assigned_to = demuxlet[demuxlet.assignment=='SNG'].loc[:,['nucleus', 'sng']].drop_duplicates().groupby('nucleus').size()
nuclei_with_clear_assignment = number_individuals_assigned_to[number_individuals_assigned_to==1].index.to_list()


# In[96]:


singlet_assignments = {r['nucleus']: r['sng'] for i, r in demuxlet[(demuxlet.assignment=='SNG') & (demuxlet.nucleus.isin(nuclei_with_clear_assignment))].iterrows()}


# In[82]:

tmp = metrics.merge(thresholds, on='library')
tmp = tmp[tmp.umis!='None']
tmp.umis = tmp.umis.astype(int)
pass_nondemuxlet_qc_thresholds = tmp[(tmp.umis>=tmp.min_hqaa) & (tmp.umis<=tmp.max_hqaa) & (tmp.fraction_mitochondrial<=tmp.max_mitochondrial)].nucleus.to_list()


demuxlet['for_comparison'] = demuxlet.best.map(lambda x: 'AMB/DBL' if 'SNG' not in x else x)
masked_vs_unmasked = demuxlet.loc[demuxlet.nucleus.isin(pass_nondemuxlet_qc_thresholds),['for_comparison', 'masked', 'nucleus']].pivot(index='nucleus', columns='masked', values='for_comparison').groupby(['masked', 'unmasked']).size().reset_index().rename(columns={0: 'n'})

confusion_matrix = masked_vs_unmasked.pivot(index='unmasked', columns='masked', values='n').fillna(0)
ax = sns.heatmap(confusion_matrix, annot=True)
b, t = ax.get_ylim()
b += 0.5
t -= 0.5
ax.set_ylim(b, t)
plt.savefig('masked-vs-unmasked-confusion-matrix.png')
plt.clf()


metrics = metrics.merge(thresholds, on='library')
metrics = metrics[metrics.umis!='None']
metrics.umis = metrics.umis.astype(int)
metrics = metrics[metrics.umis>=100]
metrics.head()


# make sure we ran demuxlet for all barcodes meeting the min UMI threshold
metrics_from_demuxlet_libraries = metrics[metrics.library.isin(demuxlet.library.unique())]
metrics_from_demuxlet_libraries = metrics_from_demuxlet_libraries[metrics_from_demuxlet_libraries.umis>=metrics_from_demuxlet_libraries.min_hqaa]
assert(all(metrics_from_demuxlet_libraries.nucleus.isin(demuxlet_tested)))


# In[87]:


metrics['individual'] = metrics.library.map(lambda x: LIBRARY_TO_INDIVIDUALS[x][0] if len(LIBRARY_TO_INDIVIDUALS[x]) == 1 else '')
metrics.loc[metrics.individual=='','individual'] = metrics.loc[metrics.individual=='',:].nucleus.map(lambda x: singlet_assignments[x] if x in singlet_assignments else 'DBL')


# enforce thresholds
pass_qc = metrics[(metrics.umis>=metrics.min_hqaa) & (metrics.umis<=metrics.max_hqaa) & (metrics.fraction_mitochondrial<=metrics.max_mitochondrial) & (metrics.individual!='DBL')]
pass_qc = pass_qc[['library', 'barcode', 'individual']]
pass_qc.to_csv('pass-qc.tsv', sep='\t', index=False, header=False)


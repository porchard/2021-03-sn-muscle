#!/usr/bin/env python
# coding: utf-8

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

# ATAC
#DEMUXLET = ['/lab/work/porchard/sn-muscle-project/work/qc/results/demuxlet/unmasked/63_20-hg19.best.txt', '/lab/work/porchard/sn-muscle-project/work/qc/results/demuxlet/unmasked/63_40-hg19.best.txt']
#THRESHOLDS = '/lab/work/porchard/sn-muscle-project/initial-thresholds-atac.txt'
#METRICS = '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/nucleus-qc/metrics.txt'
#LIBRARY_LABELS = '/lab/work/porchard/sn-muscle-project/library-labels.txt'
SPECIES_THRESHOLD = 0.87
REMOVE_HQAA_ABOVE_QUANTILE_SINGLE_SPECIES = 0.75 # after all other QC filtering -- if a barcode, is in the 75th percentile or above in terms of HQAA for it's library, remove it (doublets often have high HQAA)
REMOVE_HQAA_ABOVE_QUANTILE_DUAL_SPECIES = 0.9

LIBRARY_TO_INDIVIDUALS = {
    '125589-rn6': ['rat1'],
    '125589-hg19': ['KSM1'],
    '133151-hg19': ['KSM1'],
    '133152-hg19': ['KSM1'],
    '133153-hg19': ['KSM1'],
    '133154-hg19': ['KSM1'],
    '63_20-hg19': ['KSM1', 'KSM2'],
    '63_40-hg19': ['KSM1', 'KSM2']
}


# In[ ]:


parser = argparse.ArgumentParser('Enforce QC thresholds to call ATAC nuclei.')
parser.add_argument('--demuxlet', nargs='+', help='Demuxlet best files (unmaked)')
parser.add_argument('--metrics', help='Ataqv metric file')
parser.add_argument('--thresholds', help='Fle of initial thresholds')
parser.add_argument('--library-labels', dest='library_labels', help='Ataqv metric file')
args = parser.parse_args()

DEMUXLET = args.demuxlet
THRESHOLDS = args.thresholds
METRICS = args.metrics
LIBRARY_LABELS = args.library_labels


# In[239]:


@ticker.FuncFormatter
def read_count_formatter(x, pos):
    if x >= 1e9:
        return '{}B'.format(x/1e9)
    if x >= 1e6:
        return '{}M'.format(x/1e6)
    if x >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


# In[215]:


def load_demuxlet(f):
    tmp = pd.read_csv(f, sep='\t').loc[:,['BARCODE', 'BEST', 'SNG.1ST']].rename(columns={'BARCODE': 'barcode', 'BEST': 'best', 'SNG.1ST': 'sng'})
    tmp['assignment'] = tmp.best.map(lambda x: x.split('-')[0])
    tmp['library'] = os.path.basename(f).replace('.best.txt', '')
    assert(all(tmp.assignment.isin(['SNG', 'DBL', 'AMB'])))
    return tmp


# In[ ]:


library_labels = pd.read_csv(LIBRARY_LABELS, sep='\t')
FANS_LIBRARIES = library_labels[library_labels.fans_status=='FANS'].library.unique().tolist()
library_labels = dict(zip(library_labels.library, library_labels.name))


# In[228]:


demuxlet = pd.concat([load_demuxlet(f) for f in DEMUXLET])
demuxlet['nucleus'] = demuxlet.library + '-' + demuxlet.barcode
nucleus_to_doublet_singlet_call = dict(zip(demuxlet.nucleus, demuxlet.assignment))
nucleus_to_individual = dict(zip(demuxlet.nucleus, demuxlet.sng))
nonsinglets = demuxlet[demuxlet.assignment!='SNG'].nucleus.to_list()
demuxlet.head()


# In[217]:


thresholds = pd.read_csv(THRESHOLDS, sep='\t')


# In[218]:


metrics = pd.read_csv(METRICS, sep='\t', header=None, names=['nucleus', 'metric', 'value'])
metrics['base_library'] = metrics.nucleus.map(lambda x: x.split('-')[0])
metrics['genome'] = metrics.nucleus.map(lambda x: x.split('-')[1])
metrics['barcode'] = metrics.nucleus.map(lambda x: x.split('-')[2])


# In[219]:


# assign species
dual_species_libraries = metrics[['base_library', 'genome']].drop_duplicates().base_library.value_counts().where(lambda x: x>1).dropna().index.to_list()
species_comparison = metrics[(metrics.metric=='hqaa') & (metrics.base_library.isin(dual_species_libraries))].pivot(index=['base_library', 'barcode'], columns='genome', values='value').astype(int)
species_comparison['assignment'] = species_comparison.apply(lambda x: x.idxmax() if (x.max() / x.sum()) >= SPECIES_THRESHOLD else 'None', axis=1)
species_comparison = species_comparison.reset_index()
keep_nuclei = species_comparison.loc[species_comparison.assignment!='None',['base_library', 'assignment', 'barcode']].apply(lambda x: '-'.join(x), axis=1).to_list()
# drop those nuclei from dual species libraries that were assigned to the other species
metrics = metrics[(~metrics.base_library.isin(dual_species_libraries)) | (metrics.nucleus.isin(keep_nuclei))]
# finished assigning species


# In[220]:


metrics = metrics[['nucleus', 'metric', 'value']].pivot(index='nucleus', columns='metric', values='value')
metrics = metrics[metrics.tss_enrichment!='None'].astype(float).reset_index()
metrics['library'] = metrics.nucleus.map(lambda x: '-'.join(x.split('-')[:2]))
# enfore initial thresholds
metrics = metrics.merge(thresholds, on='library')
metrics['pass_initial_qc'] = ((metrics.hqaa>=metrics.min_hqaa) & (metrics.hqaa<=metrics.max_hqaa) & (metrics.tss_enrichment>=metrics.min_tss_enrichment) & (metrics.max_fraction_reads_from_single_autosome <= metrics.max_max_fraction_reads_from_single_autosome))


# In[222]:


# compare to demuxlet
metrics_with_demuxlet = metrics[metrics.nucleus.isin(nucleus_to_doublet_singlet_call)]
metrics_with_demuxlet = metrics_with_demuxlet[metrics_with_demuxlet.pass_initial_qc]
metrics_with_demuxlet['demuxlet'] = metrics_with_demuxlet.nucleus.map(nucleus_to_doublet_singlet_call)
metrics_with_demuxlet.head()


# In[124]:


fig, ax = plt.subplots()
sns.boxplot(x='library', y='hqaa', hue='demuxlet', data=metrics_with_demuxlet, ax=ax)
ax.set_yscale('log')
fig.tight_layout()
fig.savefig('doublet-vs-singlet-hqaa_boxplot.png')
fig.clf()


# In[184]:


# for each of the demuxlet libraries, plot out cumulative fraction doublets...
for l in metrics_with_demuxlet.library.unique():
    tmp = metrics_with_demuxlet[metrics_with_demuxlet.library==l].sort_values('hqaa', ascending=False).loc[:,['demuxlet', 'hqaa']]
    tmp['is_doublet'] = (tmp.demuxlet=='DBL').astype(int)
    tmp['cumulative_doublet_proportion'] = tmp.is_doublet.expanding().mean()
    tmp['window_doublet_proportion'] = tmp.is_doublet.rolling(window=100, min_periods=100, center=True).mean()
    tmp = tmp.sort_values('hqaa', ascending=True)
    tmp['rank'] = range(1, len(tmp)+1)
    tmp = tmp.sort_values('rank', ascending=False)

    # identify the point beyond which the proportion is always > 0.25
    rank_threshold = tmp.set_index('rank').window_doublet_proportion.dropna().map(lambda x: int(x > 0.25)).expanding().mean().where(lambda x: x>=1.0).dropna().index.to_series().min()
    quantile_threshold = 1 - (tmp['rank']>=rank_threshold).mean()

    fig, ax = plt.subplots(ncols=2, figsize=(10,4))
    ax[0].plot(tmp['rank'], tmp.cumulative_doublet_proportion)
    ax[1].plot(tmp['rank'], tmp.window_doublet_proportion)
    ax[0].set_ylim(bottom=0, top=1)
    ax[1].set_ylim(bottom=0, top=1)
    ax[1].axhline(0.25, color='red', linestyle='dashed')
    ax[1].axvline(rank_threshold, color='red', linestyle='dashed')
    ax[1].text(x=rank_threshold, y=0.8, s='{}%'.format(round(100*quantile_threshold, 1)))
    ax[0].set_xlabel('HQAA rank (larger = higher HQAA)')
    ax[0].set_ylabel('Prob(nucleus with greater HQAA is doublet)')
    ax[1].set_xlabel('HQAA rank (larger = higher HQAA)')
    ax[1].set_ylabel('Prob(nucleus with similar HQAA rank is doublet)\n(Rolling windows of 100 nuclei)')
    fig.tight_layout()
    fig.savefig(f'{l}.doublet-vs-hqaa-rank.png')
    fig.clf()


# In[223]:


# make plots...
library_to_hqaa_quantile = {library: REMOVE_HQAA_ABOVE_QUANTILE_DUAL_SPECIES if library.split('-')[0] in dual_species_libraries else REMOVE_HQAA_ABOVE_QUANTILE_SINGLE_SPECIES for library in thresholds.library}

new_max_hqaa_thresholds = dict()
for library, grp in metrics[metrics.pass_initial_qc].groupby('library'):
    new_max_hqaa_thresholds[library] = grp.hqaa.quantile(library_to_hqaa_quantile[library], interpolation='nearest')
new_max_hqaa_thresholds


# In[224]:


thresholds['new_max_hqaa'] = thresholds.library.map(new_max_hqaa_thresholds)


# In[230]:


metrics['new_max_hqaa'] = metrics.library.map(new_max_hqaa_thresholds)
metrics['pass_final_qc'] = ((metrics.hqaa>=metrics.min_hqaa) & (metrics.hqaa<=metrics.max_hqaa) & (metrics.tss_enrichment>=metrics.min_tss_enrichment) & (metrics.max_fraction_reads_from_single_autosome <= metrics.max_max_fraction_reads_from_single_autosome) & (~metrics.nucleus.isin(nonsinglets)))


# In[308]:


# make plots
def plot_hqaa_vs_tss_enrichment(df, thresholds, library_labels):
    number_libraries = df.library.nunique()
    COL_WRAP = 3
    NROWS = int(number_libraries/COL_WRAP) + int(divmod(number_libraries, COL_WRAP)[1] > 1)
    NCOLS = min(number_libraries, COL_WRAP)
    fig, axs = plt.subplots(nrows=NROWS, ncols=NCOLS, figsize=(4.5*NCOLS, 4.5*NROWS))
    for library, ax in zip(df.library.unique(), axs.flatten()):
        ax.scatter(x='hqaa', y='tss_enrichment', data=df[df.library==library], alpha=0.02, color='black')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(bottom=0.1, top=1000)
        ax.axhline(y=thresholds.set_index('library').at[library, 'min_tss_enrichment'], color='red', linestyle='dashed')
        ax.axvline(x=thresholds.set_index('library').at[library, 'min_hqaa'], color='red', linestyle='dashed')
        #ax.axvline(x=thresholds.set_index('library').at[library, 'new_max_hqaa'], color='red', linestyle='dashed')
        ax.xaxis.set_major_formatter(read_count_formatter)
        ax.set_ylabel('TSS enrichment')
        ax.set_xlabel('Final high-quality autosomal reads')
        ax.set_title(library_labels[library])
    fig.tight_layout()
    return (fig, axs)


def plot_hqaa_vs_max_frac_autosomal(df, thresholds, library_labels):
    number_libraries = df.library.nunique()
    COL_WRAP = 3
    NROWS = int(number_libraries/COL_WRAP) + int(divmod(number_libraries, COL_WRAP)[1] > 1)
    NCOLS = min(number_libraries, COL_WRAP)
    fig, axs = plt.subplots(nrows=NROWS, ncols=NCOLS, figsize=(4.5*NCOLS, 4.5*NROWS))
    for library, ax in zip(df.library.unique(), axs.flatten()):
        ax.scatter(x='hqaa', y='max_fraction_reads_from_single_autosome', data=df[df.library==library], alpha=0.02, color='black')
        ax.set_xscale('log')
        ax.set_ylim(bottom=0, top=0.5)
        ax.axhline(y=thresholds.set_index('library').at[library, 'max_max_fraction_reads_from_single_autosome'], color='red', linestyle='dashed')
        ax.axvline(x=thresholds.set_index('library').at[library, 'min_hqaa'], color='red', linestyle='dashed')
        #ax.axvline(x=thresholds.set_index('library').at[library, 'new_max_hqaa'], color='red', linestyle='dashed')
        ax.xaxis.set_major_formatter(read_count_formatter)
        ax.set_ylabel('Max. fraction reads from single autosome')
        ax.set_xlabel('Final high-quality autosomal reads')
        ax.set_title(library_labels[library])
    fig.tight_layout()
    return (fig, axs)

# plot HQAA vs TSS enrichment
fig, ax = plot_hqaa_vs_tss_enrichment(metrics[metrics.library.isin(['63_20-hg19', '63_40-hg19'])], thresholds, library_labels)
fig.savefig('hqaa-vs-tss-enrich-loading.png')
fig.clf()
fig, ax = plot_hqaa_vs_tss_enrichment(metrics[~metrics.library.isin(FANS_LIBRARIES)], thresholds, library_labels)
fig.savefig('hqaa-vs-tss-enrich-all-used-downstream.png')
fig.clf()
fig, ax = plot_hqaa_vs_tss_enrichment(metrics, thresholds, library_labels)
fig.savefig('hqaa-vs-tss-enrich-all.png')
fig.clf()

# plot HQAA vs max frac autosomal enrichment
fig, ax = plot_hqaa_vs_max_frac_autosomal(metrics[metrics.library.isin(['63_20-hg19', '63_40-hg19'])], thresholds, library_labels)
fig.savefig('hqaa-vs-max-frac-autosomal-loading.png')
fig.clf()
fig, ax = plot_hqaa_vs_max_frac_autosomal(metrics[~metrics.library.isin(FANS_LIBRARIES)], thresholds, library_labels)
fig.savefig('hqaa-vs-max-frac-autosomal-all-used-downstream.png')
fig.clf()
fig, ax = plot_hqaa_vs_max_frac_autosomal(metrics, thresholds, library_labels)
fig.savefig('hqaa-vs-max-frac-autosomal-all.png')
fig.clf()


# In[188]:


#species_comparison[species_comparison[['hg19', 'rn6']].max(axis=1)>20000].assignment.value_counts()
#1-455/(1410+1224+455)


# In[234]:

pass_qc = metrics[metrics.pass_final_qc]
pass_qc['individual'] = pass_qc[['library', 'nucleus']].apply(lambda x: LIBRARY_TO_INDIVIDUALS[x[0]][0] if len(LIBRARY_TO_INDIVIDUALS[x[0]]) == 1 else nucleus_to_individual[x[1]], axis=1)
pass_qc = pass_qc[~pass_qc.nucleus.isin(nonsinglets)]
pass_qc['barcode'] = pass_qc.nucleus.map(lambda x: x.split('-')[-1])
pass_qc['library'] = pass_qc.nucleus.map(lambda x: '-'.join(x.split('-')[0:-1]))
pass_qc = pass_qc[['library', 'barcode', 'individual']]
pass_qc.to_csv('pass-qc.tsv', sep='\t', index=False, header=False)


# In[ ]:
thresholds[['library', 'min_hqaa', 'max_max_fraction_reads_from_single_autosome', 'min_tss_enrichment']].to_csv('atac-qc-thresholds.tsv', sep='\t', index=False, header=True)

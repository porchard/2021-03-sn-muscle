#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import glob
import re


# In[2]:


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

    
def load_ataqv_metric_file(f):
    tmp = pd.read_csv(f, sep='\t', header=None, names=['nucleus', 'metric', 'value'])
    tmp['base_library'] = tmp.nucleus.map(lambda x: x.split('-')[0])
    tmp['genome'] = tmp.nucleus.map(lambda x: x.split('-')[1])
    tmp['barcode'] = tmp.nucleus.map(lambda x: x.split('-')[2])
    return tmp


def load_rnaqc_metric_file(f):
    FILENAME_RE = '(.*)-(.*).qc.txt'
    library, genome = re.match(FILENAME_RE, os.path.basename(f)).groups()
    tmp = pd.read_csv(f, sep='\t')
    tmp = tmp[tmp.barcode!='no_barcode'].set_index('barcode')
    for col in tmp.columns:
        if 'alignment' in col or 'reads' in col or col == 'umis':
            tmp[col] = tmp[col].astype(int)
        elif 'fraction' in col:
            tmp[col] = tmp[col].astype(float)
    tmp['base_library'] = library
    tmp['genome'] = genome
    tmp = tmp.reset_index()
    return tmp


# In[3]:


ATAC_BARCODE_LIST = '/home/porchard/github/snATACseq-NextFlow/737K-arc-v1.txt'
RNA_BARCODE_LIST = '/home/porchard/github/snRNAseq-NextFlow/737K-arc-v1.txt'
ATAQV_METRIC_FILES = glob.glob('/lab/work/porchard/Nova-303/work/qc/results/ataqv/*')
RNAQC_METRIC_FILES = glob.glob('/lab/work/porchard/Nova-315/work/starsolo-qc/results/qc/*')


# In[4]:


atac_barcodes = pd.read_csv(ATAC_BARCODE_LIST, header=None, names=['atac_barcode'])
rna_barcodes = pd.read_csv(RNA_BARCODE_LIST, header=None, names=['rna_barcode'])
barcodes = pd.concat([atac_barcodes, rna_barcodes], axis=1)
barcodes['joint'] = barcodes.rna_barcode + '_' + barcodes.atac_barcode
barcode_rna_to_joint = dict(zip(barcodes.rna_barcode, barcodes.joint))
barcode_atac_to_joint = dict(zip(barcodes.atac_barcode, barcodes.joint))
barcodes.head()


# In[5]:


# read in QC metrics
ataqv_metrics = pd.concat([load_ataqv_metric_file(f) for f in ATAQV_METRIC_FILES])
ataqv_metrics.head()


# In[6]:


rnaqc_metrics = pd.concat([load_rnaqc_metric_file(f) for f in RNAQC_METRIC_FILES])
rnaqc_metrics.head()


# In[7]:


human_vs_rat_atac = ataqv_metrics.loc[ataqv_metrics.metric=='hqaa',['value', 'genome', 'barcode']].pivot(index='barcode', columns='genome', values='value').astype(int)
human_vs_rat_atac.head()


# In[8]:


human_vs_rat_rna = rnaqc_metrics.loc[:,['umis', 'genome', 'barcode']].pivot(index='barcode', columns='genome', values='umis').fillna(0).astype(int)
human_vs_rat_rna.head()


# In[9]:


# plot hg19 vs rn6 for RNA and for ATAC separately
# left col: rna
# right col: atac
# top row: rn6 hqaa vs hg19 hqaa (or UMIs, for RNA)
# bottom row: hg19 / (rn6 + hg19)
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))

# top left
ax = axs[0,0]
sns.scatterplot(x='hg19', y='rn6', data=human_vs_rat_rna, alpha=0.1, ax=ax)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_ylabel('Barcode UMIs (rn6)')
ax.set_xlabel('Barcode UMIs (hg19)')
ax.set_title('RNA')

# top right
rn6_vs_hg19 = axs[0,1]
sns.scatterplot(x='hg19', y='rn6', data=human_vs_rat_atac, alpha=0.1, ax=rn6_vs_hg19)
rn6_vs_hg19.xaxis.set_major_formatter(read_count_formatter)
rn6_vs_hg19.yaxis.set_major_formatter(read_count_formatter)
rn6_vs_hg19.set_ylabel('Barcode HQAA (rn6)')
rn6_vs_hg19.set_xlabel('Barcode HQAA (hg19)')
rn6_vs_hg19.set_title('ATAC')


# bottom left
tmp = human_vs_rat_rna[human_vs_rat_rna.sum(axis=1)>0]
tmp['x'] = tmp.hg19
tmp['y'] = tmp.hg19 / tmp[['hg19', 'rn6']].sum(axis=1)
ax = axs[1,0]
sns.scatterplot(x='x', y='y', data=tmp, alpha=0.1, ax=ax)
ax.set_xlabel('UMI when mapping to hg19')
ax.set_ylabel('(UMI when mapping to hg19) /\n(UMI when mapping to hg19 + UMI when mapping to rn6)')
ax.set_xscale('log')
ax.set_xlim(left=1)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_title('RNA')

# bottom right
tmp = human_vs_rat_atac
tmp['x'] = tmp.hg19
tmp['y'] = tmp.hg19 / tmp[['hg19', 'rn6']].sum(axis=1)
fraction_hg19_ax = axs[1,1]
sns.scatterplot(x='x', y='y', data=tmp, alpha=0.1, ax=fraction_hg19_ax)
fraction_hg19_ax.set_xlabel('HQAA when mapping to hg19')
fraction_hg19_ax.set_ylabel('(HQAA when mapping to hg19) /\n(HQAA when mapping to hg19 + HQAA when mapping to rn6)')
fraction_hg19_ax.set_xscale('log')
fraction_hg19_ax.xaxis.set_major_formatter(read_count_formatter)
fraction_hg19_ax.set_title('ATAC')
fig.tight_layout()
fig.savefig('species-pairplot.png')
fig.clf()


# In[12]:


# read in ataqv metrics
# for rat and human separately, plot counts in RNA vs counts in ATAC
rna_vs_atac_counts = pd.DataFrame({'barcode': barcodes.joint})
rna_vs_atac_counts['rna_barcode'] = rna_vs_atac_counts.barcode.map(lambda x: x.split('_')[0])
rna_vs_atac_counts['atac_barcode'] = rna_vs_atac_counts.barcode.map(lambda x: x.split('_')[1])
rna_vs_atac_counts['atac_hg19'] = rna_vs_atac_counts.atac_barcode.map(lambda x: human_vs_rat_atac.at[x,'hg19'] if x in human_vs_rat_atac.index else 0)
rna_vs_atac_counts['atac_rn6'] = rna_vs_atac_counts.atac_barcode.map(lambda x: human_vs_rat_atac.at[x,'rn6'] if x in human_vs_rat_atac.index else 0)
rna_vs_atac_counts['rna_hg19'] = rna_vs_atac_counts.rna_barcode.map(lambda x: human_vs_rat_rna.at[x,'hg19'] if x in human_vs_rat_rna.index else 0)
rna_vs_atac_counts['rna_rn6'] = rna_vs_atac_counts.rna_barcode.map(lambda x: human_vs_rat_rna.at[x,'rn6'] if x in human_vs_rat_rna.index else 0)
#maxes = rna_vs_atac_counts[['atac_hg19', 'atac_rn6', 'rna_hg19', 'rna_rn6']].max(axis=1)
#rna_vs_atac_counts = rna_vs_atac_counts[maxes>10]
PASS_RNA_THRESHOLD = rna_vs_atac_counts[['rna_hg19', 'rna_rn6']].max(axis=1)>=100
PASS_ATAC_THRESHOLD = rna_vs_atac_counts[['atac_hg19', 'atac_rn6']].max(axis=1)>=1000
rna_vs_atac_counts = rna_vs_atac_counts[(PASS_RNA_THRESHOLD) | (PASS_ATAC_THRESHOLD)]
rna_vs_atac_counts.head()


# In[35]:


HUMAN_MIN_ATAC = 5000
HUMAN_MIN_RNA = 100
RAT_MIN_ATAC = 5000
RAT_MIN_RNA = 100


# In[31]:


fig, axs = plt.subplots(ncols=2, figsize=(10, 5))
# plot hg19 on left
ax = axs[0]
ax.scatter(x='rna_hg19', y='atac_hg19', data=rna_vs_atac_counts, alpha=0.05)
#sns.kdeplot(x='rna_hg19', y='atac_hg19', data=rna_vs_atac_counts.iloc[0:1000,:], ax=ax, color='red')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=1)
ax.set_ylim(bottom=1)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('ATAC HQAA')
ax.axvline(HUMAN_MIN_RNA, color='red', linestyle='dashed')
ax.axhline(HUMAN_MIN_ATAC, color='red', linestyle='dashed')
ax.set_title('Human')

# plot rn6 on the right
ax = axs[1]
ax.scatter(x='rna_rn6', y='atac_rn6', data=rna_vs_atac_counts, alpha=0.05)
#sns.kdeplot(x='rna_hg19', y='atac_hg19', data=rna_vs_atac_counts.iloc[0:1000,:], ax=ax, color='red')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=1)
ax.set_ylim(bottom=1)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('ATAC HQAA')
ax.axvline(RAT_MIN_RNA, color='red', linestyle='dashed')
ax.axhline(RAT_MIN_ATAC, color='red', linestyle='dashed')
ax.set_title('Rat')

fig.tight_layout()
fig.savefig('read-filter-before-species-calling.png')
fig.clf()


# In[36]:


# assign species
rna_vs_atac_counts['atac_ratio'] = rna_vs_atac_counts.atac_hg19 / rna_vs_atac_counts[['atac_hg19', 'atac_rn6']].sum(axis=1)
rna_vs_atac_counts['rna_ratio'] = rna_vs_atac_counts.rna_hg19 / rna_vs_atac_counts[['rna_hg19', 'rna_rn6']].sum(axis=1)
rna_vs_atac_counts.head()
fig, ax = plt.subplots()
ax.scatter(x='rna_ratio', y='atac_ratio', data=rna_vs_atac_counts, alpha=0.05)
ax.set_xlabel('hg19 UMI / (hg19 UMI + rn6 UMI)')
ax.set_ylabel('hg19 HQAA / (hg19 HQAA + rn6 HQAA)')
ax.axhline(0.85, color='red', linestyle='dashed')
ax.axhline(0.15, color='red', linestyle='dashed')
fig.tight_layout()
fig.savefig('human-rna-vs-atac-fractions.png')
fig.clf()


fig, ax = plt.subplots()
meets_counts_threshold_hg19 = ((rna_vs_atac_counts.atac_hg19 >= HUMAN_MIN_ATAC) & (rna_vs_atac_counts.rna_hg19 >= HUMAN_MIN_RNA))
meets_counts_threshold_rn6 = ((rna_vs_atac_counts.atac_rn6 >= RAT_MIN_ATAC) & (rna_vs_atac_counts.rna_rn6 >= RAT_MIN_RNA))
ax.scatter(x='rna_ratio', y='atac_ratio', data=rna_vs_atac_counts[(meets_counts_threshold_hg19) | (meets_counts_threshold_rn6)], alpha=0.05)
ax.set_xlabel('hg19 UMI / (hg19 UMI + rn6 UMI)')
ax.set_ylabel('hg19 HQAA / (hg19 HQAA + rn6 HQAA)')
ax.axhline(0.85, color='red', linestyle='dashed')
ax.axhline(0.15, color='red', linestyle='dashed')
fig.tight_layout()
fig.savefig('species-calling.human-rna-vs-atac-fractions-meeting-count-thresholds.png')
fig.clf()


# In[37]:


rna_vs_atac_counts['is_nucleus_hg19'] = ((rna_vs_atac_counts.atac_hg19 >= HUMAN_MIN_ATAC) & (rna_vs_atac_counts.rna_hg19 >= HUMAN_MIN_RNA) & (rna_vs_atac_counts.atac_ratio>=0.85))
rna_vs_atac_counts['is_nucleus_rn6'] = ((rna_vs_atac_counts.atac_rn6 >= RAT_MIN_ATAC) & (rna_vs_atac_counts.rna_rn6 >= RAT_MIN_RNA) & (rna_vs_atac_counts.atac_ratio<=0.15))
rna_vs_atac_counts.head()


# In[38]:


# plot human vs rat ATAC HQAA for supplement
def spec(x):
    if x > 0.85:
        return 'human'
    elif x < 0.15:
        return 'rat'
    else:
        return 'mixed'
tmp = rna_vs_atac_counts.copy()
tmp['species'] = tmp.atac_ratio.map(spec)
fig, ax = plt.subplots()
rn6_vs_hg19 = ax
sns.scatterplot(x='atac_hg19', y='atac_rn6', hue='species', data=tmp, alpha=0.1, ax=rn6_vs_hg19)
rn6_vs_hg19.xaxis.set_major_formatter(read_count_formatter)
rn6_vs_hg19.yaxis.set_major_formatter(read_count_formatter)
rn6_vs_hg19.set_ylabel('Pass filter ATAC reads\n(rat genome)')
rn6_vs_hg19.set_xlabel('Pass filter ATAC reads\n(human genome)')
fig.tight_layout()
fig.savefig('human-vs-rat-atac-hqaa.png')
fig.clf()


# In[39]:


human_nuclei = rna_vs_atac_counts.loc[rna_vs_atac_counts.is_nucleus_hg19,['barcode', 'rna_barcode', 'atac_barcode']].reset_index(drop=True)
rat_nuclei = rna_vs_atac_counts.loc[rna_vs_atac_counts.is_nucleus_rn6,['barcode', 'rna_barcode', 'atac_barcode']].reset_index(drop=True)

ataqv_metrics_human = ataqv_metrics.loc[ataqv_metrics.genome=='hg19',['barcode', 'metric', 'value']].pivot(index='barcode', columns='metric', values='value').reset_index().rename(columns={'barcode': 'atac_barcode'})
ataqv_metrics_human = ataqv_metrics_human[['atac_barcode', 'hqaa', 'tss_enrichment', 'max_fraction_reads_from_single_autosome']]
ataqv_metrics_human = ataqv_metrics_human[ataqv_metrics_human.atac_barcode.isin(human_nuclei.atac_barcode)]
ataqv_metrics_human.hqaa = ataqv_metrics_human.hqaa.astype(int)
ataqv_metrics_human.tss_enrichment = ataqv_metrics_human.tss_enrichment.astype(float)
ataqv_metrics_human.max_fraction_reads_from_single_autosome = ataqv_metrics_human.max_fraction_reads_from_single_autosome.astype(float)
human_nuclei = human_nuclei.merge(ataqv_metrics_human, on='atac_barcode')

ataqv_metrics_rat = ataqv_metrics.loc[ataqv_metrics.genome=='rn6',['barcode', 'metric', 'value']].pivot(index='barcode', columns='metric', values='value').reset_index().rename(columns={'barcode': 'atac_barcode'})
ataqv_metrics_rat = ataqv_metrics_rat[['atac_barcode', 'hqaa', 'tss_enrichment', 'max_fraction_reads_from_single_autosome']]
ataqv_metrics_rat = ataqv_metrics_rat[ataqv_metrics_rat.atac_barcode.isin(rat_nuclei.atac_barcode)]
ataqv_metrics_rat.hqaa = ataqv_metrics_rat.hqaa.astype(int)
ataqv_metrics_rat.tss_enrichment = ataqv_metrics_rat.tss_enrichment.astype(float)
ataqv_metrics_rat.max_fraction_reads_from_single_autosome = ataqv_metrics_rat.max_fraction_reads_from_single_autosome.astype(float)
rat_nuclei = rat_nuclei.merge(ataqv_metrics_rat, on='atac_barcode')

rnaqc_metrics_tmp = rnaqc_metrics.rename(columns={'barcode': 'rna_barcode', 'fraction_mitochondrial': 'rna_fraction_mitochondrial'}).loc[:,['genome', 'rna_barcode', 'umis', 'rna_fraction_mitochondrial']]
human_nuclei = human_nuclei.merge(rnaqc_metrics_tmp[rnaqc_metrics_tmp.genome=='hg19'], on='rna_barcode')
rat_nuclei = rat_nuclei.merge(rnaqc_metrics_tmp[rnaqc_metrics_tmp.genome=='rn6'], on='rna_barcode')


# In[44]:


# plot HQAA vs UMIs
# plot HQAA vs TSS enrichment
# plot HQAA vs max fraction autosomal
# plot UMIs vs fraction mitochondrial
# human on left, rat on right

HUMAN_HQAA_MIN = 20000
RAT_HQAA_MIN = 20000
HUMAN_UMIS_MIN = 500
RAT_UMIS_MIN = 500
HUMAN_TSS_ENRICHMENT_MIN = 2
RAT_TSS_ENRICHMENT_MIN = 2
HUMAN_MAX_FRAC_FROM_AUTOSOME_MAX = 0.15
RAT_MAX_FRAC_FROM_AUTOSOME_MAX = 0.15
HUMAN_RNA_FRAC_MITO_MAX = 0.02
RAT_RNA_FRAC_MITO_MAX = 0.02

human_nuclei['pass_qc'] = (human_nuclei.hqaa >= HUMAN_HQAA_MIN) & (human_nuclei.umis >= HUMAN_UMIS_MIN) & (human_nuclei.tss_enrichment >= HUMAN_TSS_ENRICHMENT_MIN) & (human_nuclei.max_fraction_reads_from_single_autosome<=HUMAN_MAX_FRAC_FROM_AUTOSOME_MAX) & (human_nuclei.rna_fraction_mitochondrial<=HUMAN_RNA_FRAC_MITO_MAX)
rat_nuclei['pass_qc'] = (rat_nuclei.hqaa >= RAT_HQAA_MIN) & (rat_nuclei.umis >= RAT_UMIS_MIN) & (rat_nuclei.tss_enrichment >= RAT_TSS_ENRICHMENT_MIN) & (rat_nuclei.max_fraction_reads_from_single_autosome<=RAT_MAX_FRAC_FROM_AUTOSOME_MAX) & (rat_nuclei.rna_fraction_mitochondrial<=RAT_RNA_FRAC_MITO_MAX)
human_nuclei.pass_qc = human_nuclei.pass_qc.map(lambda x: '{} (n={})'.format(x, sum(human_nuclei.pass_qc==x)))
rat_nuclei.pass_qc = rat_nuclei.pass_qc.map(lambda x: '{} (n={})'.format(x, sum(rat_nuclei.pass_qc==x)))

fig, axs = plt.subplots(ncols=2, nrows=4, figsize=(2*5, 4*5))

# HQAA vs UMIs, human
ax = axs[0,0]
#sns.scatterplot(x='umis', y='hqaa', hue='pass_qc', data=human_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='umis', y='hqaa', color='black', data=human_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=HUMAN_MIN_RNA)
ax.set_ylim(bottom=HUMAN_MIN_ATAC)
#ax.set_xlim(left=1)
#ax.set_ylim(bottom=1)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('ATAC pass filter reads')
ax.set_title('Human')
ax.axvline(HUMAN_UMIS_MIN, color='red', linestyle='dashed')
ax.axhline(HUMAN_HQAA_MIN, color='red', linestyle='dashed')

# HQAA vs UMIs, rat
ax = axs[0,1]
#sns.scatterplot(x='umis', y='hqaa', hue='pass_qc', data=rat_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='umis', y='hqaa', color='black', data=rat_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('ATAC pass filter reads')
ax.set_title('Rat')
ax.axvline(RAT_UMIS_MIN, color='red', linestyle='dashed')
ax.axhline(RAT_HQAA_MIN, color='red', linestyle='dashed')

# HQAA vs TSS enrichment, human
ax = axs[1,0]
#sns.scatterplot(x='hqaa', y='tss_enrichment', hue='pass_qc', data=human_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='hqaa', y='tss_enrichment', color='black', data=human_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('ATAC pass filter reads')
ax.set_ylabel('ATAC TSS enrichment')
ax.axvline(HUMAN_HQAA_MIN, color='red', linestyle='dashed')
ax.axhline(HUMAN_TSS_ENRICHMENT_MIN, color='red', linestyle='dashed')

# HQAA vs TSS enrichment, rat
ax = axs[1,1]
#sns.scatterplot(x='hqaa', y='tss_enrichment', hue='pass_qc', data=rat_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='hqaa', y='tss_enrichment', color='black', data=rat_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('ATAC pass filter reads')
ax.set_ylabel('ATAC TSS enrichment')
ax.axvline(RAT_HQAA_MIN, color='red', linestyle='dashed')
ax.axhline(RAT_TSS_ENRICHMENT_MIN, color='red', linestyle='dashed')


# HQAA vs max fraction autosomal, human
ax = axs[2,0]
#sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', hue='pass_qc', data=human_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', color='black', data=human_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('ATAC pass filter reads')
ax.set_ylabel('Max. fraction reads from single autosome')
ax.axvline(HUMAN_HQAA_MIN, color='red', linestyle='dashed')
ax.axhline(HUMAN_MAX_FRAC_FROM_AUTOSOME_MAX, color='red', linestyle='dashed')


# HQAA vs max fraction autosomal, rat
ax = axs[2,1]
#sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', hue='pass_qc', data=rat_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', color='black', data=rat_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('ATAC pass filter reads')
ax.set_ylabel('Max. fraction reads from single autosome')
ax.axvline(RAT_HQAA_MIN, color='red', linestyle='dashed')
ax.axhline(RAT_MAX_FRAC_FROM_AUTOSOME_MAX, color='red', linestyle='dashed')

# UMIs vs fraction mitochondrial, human
ax = axs[3,0]
#sns.scatterplot(x='umis', y='rna_fraction_mitochondrial', hue='pass_qc', data=human_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='umis', y='rna_fraction_mitochondrial', color='black', data=human_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('RNA frac. mitochondrial')
ax.axvline(HUMAN_UMIS_MIN, color='red', linestyle='dashed')
ax.axhline(HUMAN_RNA_FRAC_MITO_MAX, color='red', linestyle='dashed')

# UMIs vs fraction mitochondrial, rat
ax = axs[3,1]
#sns.scatterplot(x='umis', y='rna_fraction_mitochondrial', hue='pass_qc', data=rat_nuclei, alpha=0.3, ax=ax)
sns.scatterplot(x='umis', y='rna_fraction_mitochondrial', color='black', data=rat_nuclei, alpha=0.3, ax=ax)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('RNA UMIs')
ax.set_ylabel('RNA frac. mitochondrial')
ax.axvline(RAT_UMIS_MIN, color='red', linestyle='dashed')
ax.axhline(RAT_RNA_FRAC_MITO_MAX, color='red', linestyle='dashed')

fig.tight_layout()
fig.savefig('post-species-calling-qc.png')
fig.clf()


# In[17]:


HUMAN_PASS_QC = human_nuclei.loc[human_nuclei.pass_qc.map(lambda x: 'True' in x),['barcode', 'rna_barcode', 'atac_barcode']]
RAT_PASS_QC = rat_nuclei.loc[rat_nuclei.pass_qc.map(lambda x: 'True' in x),['barcode', 'rna_barcode', 'atac_barcode']]
HUMAN_PASS_QC.to_csv('human-nuclei.tsv', sep='\t', index=False)
RAT_PASS_QC.to_csv('rat-nuclei.tsv', sep='\t', index=False)


# In[ ]:





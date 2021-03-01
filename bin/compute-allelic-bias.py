#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
import glob
from scipy.stats import binom_test
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import re


# In[2]:


def make_quadrant_labels(x, y):
    # x and y must be pandas series/numpy arrays
    POSSIBLE_SIGNS = [-1, 0, 1]
    combos = pd.DataFrame([[i, j, sum(np.logical_and(np.sign(x)==i, np.sign(y)==j))] for i in POSSIBLE_SIGNS for j in POSSIBLE_SIGNS], columns=['x_sign', 'y_sign', 'count'])
    combos['fraction'] = combos['count'] / combos['count'].sum()
    combos['label'] = combos[['count', 'fraction']].apply(lambda x: '{} ({}%)'.format(int(x[0]), round(100*x[1], 1)), axis=1)
    return combos


# In[3]:


COUNT_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/atac-allelic-bias-by-cell-type/results/nucleotide-counts/*hg19.KSM*counts.txt')
CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
HET_POSITION_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/atac-allelic-bias-by-cell-type/results/nucleotide-counts/KSM*.het-positions.bed')

cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_id_to_cluster_name = dict(zip(cluster_names.old_name.astype(str), cluster_names.new_name))


def parse_file_name(f):
    library, sample, cluster, _, _2 = os.path.basename(f).split('.')
    return {'library': library, 'sample': sample, 'cluster': cluster}


def load_count_file(f):
    tmp = pd.read_csv(f, sep='\t')
    for k, v in parse_file_name(f).items():
        tmp[k] = v
    return tmp


def load_het_positions_file(f):
    tmp = pd.read_csv(f, sep='\t', header=None, names=['chrom', 'start', 'end', 'allele_1', 'allele_2', 'ref', 'sample'])
    tmp['alt'] = tmp.allele_1
    tmp.loc[tmp.alt==tmp.ref,'alt'] = tmp.allele_2
    assert(all((tmp.allele_1==tmp.ref) | (tmp.allele_2==tmp.ref)))
    assert(all(tmp.alt!=tmp.ref))
    return tmp[['chrom', 'end', 'ref', 'alt', 'sample']].rename(columns={'end': 'pos'})


het_positions = pd.concat([load_het_positions_file(f) for f in HET_POSITION_FILES])

sizes = het_positions.groupby(['chrom', 'pos', 'sample']).size()
sizes.sort_values(ascending=False).head() # TODO: there are a couple of duplicates. Remove these, or handle them properly.
multi_allelic_positions = sizes[sizes>1]
drop = multi_allelic_positions.reset_index().loc[:,['chrom', 'pos']].astype(str).apply(lambda x: ':'.join(x), axis=1).to_list()


counts = pd.concat([load_count_file(f) for f in COUNT_FILES])
counts['variant'] = counts[['chrom', 'pos']].astype(str).apply(lambda x: ':'.join(x), axis=1)
counts['cluster_name'] = counts.cluster.map(cluster_id_to_cluster_name)
counts = counts[~counts.variant.isin(drop)]
counts = counts.merge(het_positions, on=['chrom', 'pos', 'sample'])

# combine libraries
counts_per_cluster_per_sample = counts.loc[:,['variant', 'sample', 'cluster', 'ref', 'alt', 'depth', 'nA', 'nC', 'nG', 'nT']].groupby(['variant', 'sample', 'cluster', 'ref', 'alt']).sum().reset_index()
counts_per_cluster_per_sample.to_csv('counts-per-cluster-per-sample.txt', sep='\t', index=False)

# combine samples
counts_per_cluster = counts_per_cluster_per_sample.groupby(['variant', 'cluster', 'ref', 'alt']).sum().reset_index()
counts_per_cluster.to_csv('counts-per-cluster.txt', sep='\t', index=False)


# perform the testing
MIN_DEPTH = 10
allelic_bias = counts_per_cluster[counts_per_cluster.depth>=MIN_DEPTH]
allelic_bias['ref_counts'] = allelic_bias[['ref', 'alt', 'nA', 'nC', 'nG', 'nT']].apply(lambda x: x['n' + x['ref']], axis=1)
allelic_bias['alt_counts'] = allelic_bias[['ref', 'alt', 'nA', 'nC', 'nG', 'nT']].apply(lambda x: x['n' + x['alt']], axis=1)
allelic_bias = allelic_bias[(allelic_bias.ref_counts + allelic_bias.alt_counts)>=MIN_DEPTH]
allelic_bias['fraction_ref'] = allelic_bias.ref_counts/(allelic_bias.ref_counts + allelic_bias.alt_counts)
allelic_bias['pvalue'] = allelic_bias[['ref_counts', 'alt_counts']].apply(lambda x: binom_test(x=x, p=0.5), axis=1)
allelic_bias['qvalue'] = allelic_bias.groupby('cluster').pvalue.transform(lambda x: multipletests(x, method='fdr_bh')[1])

allelic_bias.to_csv('allelic-bias.txt', sep='\t', index=False)

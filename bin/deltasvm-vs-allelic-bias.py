#!/usr/bin/env python
# coding: utf-8

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


def make_quadrant_labels(x, y):
    # x and y must be pandas series/numpy arrays
    POSSIBLE_SIGNS = [-1, 0, 1]
    combos = pd.DataFrame([[i, j, sum(np.logical_and(np.sign(x)==i, np.sign(y)==j))] for i in POSSIBLE_SIGNS for j in POSSIBLE_SIGNS], columns=['x_sign', 'y_sign', 'count'])
    combos['fraction'] = combos['count'] / combos['count'].sum()
    combos['label'] = combos[['count', 'fraction']].apply(lambda x: '{} ({}%)'.format(int(x[0]), round(100*x[1], 1)), axis=1)
    return combos

COUNT_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/atac-allelic-bias-by-cell-type/results/nucleotide-counts/*hg19.KSM*counts.txt')
CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
HET_POSITION_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/atac-allelic-bias-by-cell-type/results/nucleotide-counts/KSM*.het-positions.bed')
DELTA_SVM_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/delta-svm/all-keep-40k/results/cluster_*/snp_scores_cluster_*.txt')

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


def load_deltasvm_scores(f, variants, zscore=False):
    # get deltaSVM scores for these variants
    d = {variant: None for variant in variants}

    count = 0
    all_values = []
    with open(f, 'r') as fh:
        for line in fh:
            count += 1
            if count % 1000000 == 0:
                sys.stderr.write('Processed {} snps\n'.format(count))
            variant, score = line.rstrip().split('\t')
            variant = 'chr' + variant
            if zscore:
                all_values.append(float(score))
            if variant in d:
                d[variant] = float(score)
    
    if zscore:
        avg = np.mean(all_values)
        sd = np.std(all_values)
        z = dict()
        for variant, score in d.items():
            new_score = None if score is None else (score - avg) / sd
            z[variant] = new_score
        return z
    else:
        return d


het_positions = pd.concat([load_het_positions_file(f) for f in HET_POSITION_FILES])

het_in_n_samples = het_positions[['chrom', 'pos', 'alt', 'sample']].astype(str).drop_duplicates().groupby(['chrom', 'pos', 'alt']).size().reset_index()
het_in_n_samples['snp'] = het_in_n_samples.chrom + ':' + het_in_n_samples.pos.astype(str) + ':' + het_in_n_samples.alt
het_in_n_samples = het_in_n_samples.set_index('snp')[0]

# get deltaSVM zscores for these variants in each cell type
fetch_variants = het_in_n_samples.index.to_list()
delta_svm_z_dict = dict()

for f in DELTA_SVM_FILES:
    RE = 'snp_scores_(.*).txt'
    cluster = re.match(RE, os.path.basename(f)).group(1)
    delta_svm_z_dict[cluster] = load_deltasvm_scores(f, fetch_variants, zscore=True)


# get deltaSVM scores for these variants in each cell type
delta_svm_dict = dict()

for f in DELTA_SVM_FILES:
    RE = 'snp_scores_(.*).txt'
    cluster = re.match(RE, os.path.basename(f)).group(1)
    delta_svm_dict[cluster] = load_deltasvm_scores(f, fetch_variants)

het_in_both = [i.replace('chr', '') for i in het_in_n_samples[het_in_n_samples>1].index.to_list()]
sizes = het_positions.groupby(['chrom', 'pos', 'sample']).size()
sizes.sort_values(ascending=False).head() # TODO: there are a couple of duplicates. Remove these, or handle them properly.
multi_allelic_positions = sizes[sizes>1]
drop = multi_allelic_positions.reset_index().loc[:,['chrom', 'pos']].astype(str).apply(lambda x: ':'.join(x), axis=1).to_list()

counts = pd.concat([load_count_file(f) for f in COUNT_FILES])
counts['variant'] = counts[['chrom', 'pos']].astype(str).apply(lambda x: ':'.join(x), axis=1)
counts['cluster_name'] = counts.cluster.map(cluster_id_to_cluster_name)
counts = counts[~counts.variant.isin(drop)]
counts = counts.merge(het_positions, on=['chrom', 'pos', 'sample'])

counts_per_cluster_per_sample = counts.loc[counts.cluster.isin(['0', '1']),['variant', 'sample', 'cluster', 'ref', 'alt', 'depth', 'nA', 'nC', 'nG', 'nT']].groupby(['variant', 'sample', 'cluster', 'ref', 'alt']).sum().reset_index()

# using TII fibers:
# combine across individuals, and use strict FDR threshold
counts_per_cluster = counts_per_cluster_per_sample[counts_per_cluster_per_sample.cluster=='0'].groupby(['variant', 'ref', 'alt']).sum().reset_index()
counts_per_cluster.to_csv('counts-per-cluster.txt', sep='\t', index=False)

MIN_DEPTH = 15
allelic_bias = counts_per_cluster[counts_per_cluster.depth>=MIN_DEPTH]
allelic_bias['ref_counts'] = allelic_bias[['ref', 'alt', 'nA', 'nC', 'nG', 'nT']].apply(lambda x: x['n' + x['ref']], axis=1)
allelic_bias['alt_counts'] = allelic_bias[['ref', 'alt', 'nA', 'nC', 'nG', 'nT']].apply(lambda x: x['n' + x['alt']], axis=1)
allelic_bias = allelic_bias[(allelic_bias.ref_counts + allelic_bias.alt_counts)>=MIN_DEPTH]
allelic_bias['fraction_ref'] = allelic_bias.ref_counts/(allelic_bias.ref_counts + allelic_bias.alt_counts)
allelic_bias = allelic_bias[allelic_bias.fraction_ref>0]
allelic_bias = allelic_bias[allelic_bias.fraction_ref<1]
allelic_bias['pvalue'] = allelic_bias[['ref_counts', 'alt_counts']].apply(lambda x: binom_test(x=x, p=0.5), axis=1)
allelic_bias['qvalue'] = multipletests(allelic_bias.pvalue, method='fdr_bh')[1]
allelic_bias['significant'] = (allelic_bias.qvalue<=0.001)
allelic_bias['dsvm_score'] = (allelic_bias.variant + ':' +  allelic_bias.alt).map(delta_svm_dict['cluster_0'])
allelic_bias['dsvm_zscore'] = (allelic_bias.variant + ':' +  allelic_bias.alt).map(delta_svm_z_dict['cluster_0'])

fig, ax = plt.subplots()
compare = allelic_bias[((allelic_bias.ref_counts+allelic_bias.alt_counts)>=15) & (allelic_bias.qvalue<=0.01) & (allelic_bias.fraction_ref>0) & (allelic_bias.fraction_ref<1) & (allelic_bias.dsvm_zscore.abs() > 2)]
labels = make_quadrant_labels(compare.fraction_ref-0.5, compare.dsvm_zscore)
sns.scatterplot(x='fraction_ref', y='dsvm_zscore', data=compare, ax=ax, alpha=0.3)
ax.set_xlabel('Fraction ref. (Type II fibers)')
ax.set_ylabel('deltaSVM Z-score for alt allele (Type II fibers)')
labels['xloc'] = 0.5 + labels.x_sign * 0.4
labels['yloc'] = min([abs(i) for i in ax.get_ylim()]) * 0.8 * labels.y_sign
labels = labels[labels['count']>0]
ax.axhline(0, color='red', linestyle='dashed')
ax.axvline(0.5, color='red', linestyle='dashed')
[ax.text(x=r['xloc'], y=r['yloc'], s=r['label'], ha='center') for i, r in labels.iterrows()]
fig.tight_layout()
fig.savefig('dsvm-vs-allelic-bias-type-2-fibers-fdr-1-deltasvmZ-2-min-coverage-15.png')
fig.clf()

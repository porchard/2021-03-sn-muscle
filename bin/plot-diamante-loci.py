#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import glob
import cmap_utils
import pybedtools as bt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--peaks', nargs='+', help='Glob for cluster peak files')
parser.add_argument('--cs', nargs='+', help='Glob for DIAMANTE c.s. files')
parser.add_argument('--dsvm', nargs='+', help='Glob for DSVM files')
args = parser.parse_args()

#PEAK_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/process-by-cluster/results/peaks/broad/*hg19_peaks.broadPeak.noblacklist')
#PEAK_FILES = [f for f in PEAK_FILES if 'all-' not in f]

#DSVM_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/diamante-deltaSVM/results/dsvm/*')

#DIAMANTE_CREDIBLE_SET_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/*')

PEAK_FILES = args.peaks
DSVM_FILES = args.dsvm
DIAMANTE_CREDIBLE_SET_FILES = args.cs
CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'

cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_names.old_name = 'cluster_' + cluster_names.old_name.astype(str)
cluster_id_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))


# need the summary stats to infer the alt allele of interest
#DIAMANTE_SUMMARY_STATS = '/lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt'
#s = pd.read_csv(DIAMANTE_SUMMARY_STATS, sep='\t')
#s.head()


def load_credible_set_file(f):
    tmp = pd.read_csv(f, sep='\t')
    return tmp

cs = pd.concat([load_credible_set_file(f) for f in DIAMANTE_CREDIBLE_SET_FILES])
cs['chrom'] = 'chr' + cs.Chr.astype(str)
cs['pos'] = cs.Pos.astype(int)


dsvm = pd.concat([pd.read_csv(f, sep='\t').assign(cluster=os.path.basename(f).replace('.dsvm.txt', '')) for f in DSVM_FILES])
dsvm['chrom'] = 'chr' + dsvm.snp.map(lambda x: x.split(':')[0])
dsvm['pos'] = dsvm.snp.map(lambda x: x.split(':')[1]).astype(int)
dsvm['cluster_name'] = dsvm.cluster.map(cluster_id_to_cluster_name)


# determine whether each SNP overlaps a peak in each cluster
snp_overlaps_peak = dict() # snp --> cluster --> True/False
for f in PEAK_FILES:
    cluster = 'cluster_' + os.path.basename(f).split('-')[0]
    snp_bed = dsvm[['chrom', 'pos']].assign(start=lambda df: df.pos-1).rename(columns={'pos': 'end'}).loc[:,['chrom', 'start', 'end']].drop_duplicates()
    overlaps = bt.BedTool().from_dataframe(snp_bed).sort().intersect(bt.BedTool(f).sort(), loj=True).to_dataframe()
    overlaps_peak = overlaps.loc[overlaps.name!='.',['chrom','start', 'end']].drop_duplicates()
    overlaps_no_peak = overlaps.loc[overlaps.name=='.',['chrom', 'start', 'end']].drop_duplicates()
    assert(len(snp_bed) == (len(overlaps_no_peak) + len(overlaps_peak)))
    snp_overlaps_peak[cluster] = dict()
    for snp in overlaps_peak.chrom + ':' + overlaps_peak.end.astype(str):
        snp_overlaps_peak[cluster][snp] = True
    for snp in overlaps_no_peak.chrom + ':' + overlaps_no_peak.end.astype(str):
        snp_overlaps_peak[cluster][snp] = False


cluster_to_color = cmap_utils.make_colormap_dict(cluster_names.sort_values('old_name').new_name.to_list())

for INDEX_SNP in cs.IndexSNP.sort_values().unique():
    sys.stderr.write(f'Processing {INDEX_SNP}\n')
    sys.stderr.flush()
    INDEX_SNP_CHROM = INDEX_SNP.split('_')[0]
    INDEX_SNP_POS = int(INDEX_SNP.split('_')[1])
    # get the credible set
    ppas = cs.loc[cs.IndexSNP==INDEX_SNP,['chrom', 'pos', 'PPAg']]
    pos_to_index = {pos: index for index, pos in enumerate(ppas.pos.sort_values().unique())}
    ppas = ppas.merge(dsvm, how='left', on=['chrom', 'pos'])
    ppas['ind'] = ppas.pos.map(pos_to_index)
    INDEX_SNP_INDEX = ppas[ppas.pos==INDEX_SNP_POS]['ind'].unique()[0] if INDEX_SNP_POS in ppas.pos.to_list() else None
    ppas['snp_no_allele'] = ppas.chrom + ':' + ppas.pos.astype(str)
    ppas['snp_overlaps_peak_in_cluster'] = [snp_overlaps_peak[cluster][snp] if cluster in snp_overlaps_peak else True for cluster, snp in zip(ppas.cluster, ppas.snp_no_allele)]
    ppas['z-score-overlap'] = ppas['z-score'] * ppas['snp_overlaps_peak_in_cluster'].astype(int)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(max(ppas['ind'].max()/4, 2), 3), gridspec_kw={'hspace': 0})

    # first row is deltaSVM without peak overlap
    # second row is deltaSVM with peak overlap
    # third row is SNP PPA

    LIMIT = max(ppas['z-score'].abs().max()*1.1, 2)
    if ppas.ind.max() > 100:
        continue

    ax = axs[0]
    #g = sns.scatterplot(x='ind', y='z-score', hue='cluster', palette=cluster_to_color, data=ppas[~ppas.cluster.isnull()].sort_values('cluster'), ax=ax)
    g = sns.stripplot(x='ind', y='z-score', hue='cluster_name', palette=cluster_to_color, data=ppas[~ppas.cluster.isnull()].sort_values('cluster'), ax=ax, alpha=0.8)
#    if INDEX_SNP_INDEX is not None:
#        ax.axvline(INDEX_SNP_INDEX, color='red', linestyle='dashed', alpha=0.2)
    ax.axhline(0, color='black', linestyle='dashed', alpha=0.2)
    #g.get_legend().remove()
    lgd = g.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)
    ax.set_ylabel('DeltaSVM\nz-score')
    ax.set_xlabel('')
    ax.set_ylim(bottom=-1*LIMIT, top=LIMIT)

#    ax = axs[1]
#    g = sns.scatterplot(x='ind', y='z-score-overlap', hue='cluster', palette=cluster_to_color, data=ppas[~ppas.cluster.isnull()].sort_values('cluster'), ax=ax)
#    if INDEX_SNP_INDEX is not None:
#        ax.axvline(INDEX_SNP_INDEX, color='red', linestyle='dashed', alpha=0.2)
#    ax.axhline(2, color='red', linestyle='dashed', alpha=0.2)
#    ax.axhline(-2, color='red', linestyle='dashed', alpha=0.2)
#    g.get_legend().remove()
#    ax.set_ylabel('DeltaSVM\nz-score')
#    ax.set_xlabel('')
#    ax.set_ylim(bottom=-1*LIMIT, top=LIMIT)


    ax = axs[1]
    ax.plot('ind', 'PPAg', data=ppas[['ind', 'PPAg']].sort_values('ind').drop_duplicates())
    ax.set_ylabel('PPA')
    #if INDEX_SNP_INDEX is not None:
    #    ax.axvline(INDEX_SNP_INDEX, color='red', linestyle='dashed', alpha=0.2)
    ax.set_xticks(ppas['ind'].sort_values().unique())
    ax.set_xticklabels(ppas.sort_values('ind').pos.astype(str).unique(), rotation=90)
    ax.set_ylim(bottom=0, top=ppas.PPAg.max()*1.1)
    ax.set_xlabel(f'SNP position (chr{INDEX_SNP_CHROM})')
    axs[1].set_xlim(left=axs[0].get_xlim()[0], right=axs[0].get_xlim()[1])
    
    fig.savefig(f'{INDEX_SNP}.png', bbox_inches='tight', bbox_extra_artists=(lgd,))
    fig.clf()

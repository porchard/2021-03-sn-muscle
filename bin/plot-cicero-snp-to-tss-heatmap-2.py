#!/usr/bin/env python
# coding: utf-8

# In[49]:


import os
import sys
import pandas as pd
import pybedtools as bt
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import numpy as np


# In[88]:


#TSS = '/lab/work/porchard/data/gencode-tss/hg19.bed'
SNPS = ['chr5:53271420', 'chr12:26472562']
USE_SIGNED_DISTANCE = True
DISTANCE_TO_SHOW = 1e6
#CHROM_SIZES = '/lab/work/porchard/data/chrom_sizes/hg19.chrom_sizes'
CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
#CICERO_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/cicero/results/cicero/*.1000000.cicero.txt')
TSS = sys.argv[1]
CHROM_SIZES = sys.argv[2]
CICERO_FILES = sys.argv[3:]


# In[90]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t').astype(str)
cluster_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))


# In[19]:


def get_tss_near_snps(snp_list, tss_file, distance):
    """
    snp_list = [chr5:10000, chr5:11111, ...]
    tss_file = bed6
    distance = bps from SNP to check
    """
    snps = pd.DataFrame([i.split(':') for i in snp_list], columns=['chrom', 'end'])
    snps.end = snps.end.astype(int)
    snps['start'] = snps.end - 1
    snps_bt = bt.BedTool().from_dataframe(snps[['chrom', 'start', 'end']]).sort()
    
    nearby_tss = snps_bt.window(bt.BedTool(tss_file).sort(), w=distance).to_dataframe()
    nearby_tss['distance'] = -1*(nearby_tss.end - nearby_tss.strand)
    if not USE_SIGNED_DISTANCE:
        nearby_tss['distance'] = (nearby_tss.end - nearby_tss.strand).abs()
    nearby_tss['snp'] = nearby_tss.chrom + ':' + nearby_tss.end.astype(str)
    nearby_tss = nearby_tss[['snp', 'thickStart', 'distance']].rename(columns={'thickStart': 'gene'})
    return nearby_tss
    

def make_tss_file_of_tss_near_snps(snp_list, tss_file, distance):
    snps = pd.DataFrame([i.split(':') for i in snp_list], columns=['chrom', 'end'])
    snps.end = snps.end.astype(int)
    snps['start'] = snps.end - 1
    snps_bt = bt.BedTool().from_dataframe(snps[['chrom', 'start', 'end']]).sort()
    nearby_tss = snps_bt.window(bt.BedTool(tss_file).sort(), w=distance).to_dataframe()
    nearby_tss['distance'] = (nearby_tss.end - nearby_tss.strand).abs()
    nearby_tss = nearby_tss.loc[nearby_tss.distance<=distance,['name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']]
    nearby_tss.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return nearby_tss

# In[20]:


def cicero_assign_peaks_to_tss(peak_list, tss_file, distance, chrom_sizes, upstream_only=True):
    """
    Peak list should be format chrom_start_end
    tss_file should be bed6
    distance = distance from TSS to extend
    upstream_only = only link peak to TSS if it's upstream or overlapping the TSS
    
    Returns: Dict ( chr_start_end --> set(tss_1, tss_2, ...) )
    """
    peaks = pd.DataFrame([i.split('_') for i in peak_list], columns=['chrom', 'start', 'end'])
    peaks.start = peaks.start.astype(int)
    peaks.end = peaks.end.astype(int)
    peaks_bt = bt.BedTool().from_dataframe(peaks).sort()
    tss_bt = bt.BedTool(tss_file).sort()
    
    r = 0 if upstream_only else distance
    intersection = peaks_bt.intersect(tss_bt.slop(l=distance, r=r, s=True, g=chrom_sizes), wa=True, wb=True)
    
    peak2tss = {peak: set() for peak in peak_list} # chr_start_end --> set(tss_1, tss_2, ...)
    for i in intersection:
        peak_chrom, peak_start, peak_end, tss_chrom, tss_start, tss_end, gene, score, tss_strand = i
        peak = f'{peak_chrom}_{peak_start}_{peak_end}'
        peak2tss[peak].add(gene)
    
    return peak2tss


def cicero_assign_snps_to_peaks(snp_list, peak_list):
    """
    snp = chr5:10000
    Peak list should be format chrom_start_end
    
    Returns: Dict ( chr_pos --> set(peak_1, peak_2, ...) )
    """
    peaks = pd.DataFrame([i.split('_') for i in peak_list], columns=['chrom', 'start', 'end'])
    peaks.start = peaks.start.astype(int)
    peaks.end = peaks.end.astype(int)
    peaks_bt = bt.BedTool().from_dataframe(peaks).sort()
    
    snps = pd.DataFrame([i.split(':') for i in snp_list], columns=['chrom', 'end'])
    snps.end = snps.end.astype(int)
    snps['start'] = snps.end - 1
    snps_bt = bt.BedTool().from_dataframe(snps[['chrom', 'start', 'end']]).sort()
    
    intersection = peaks_bt.intersect(snps_bt, wa=True, wb=True)
    snp2peak = {snp: set() for snp in snp_list} # chr5:10000 --> set(peak1, peak2)
    
    for i in intersection:
        peak_chrom, peak_start, peak_end, snp_chrom, snp_start, snp_end = i
        peak = f'{peak_chrom}_{peak_start}_{peak_end}'
        snp = f'{snp_chrom}:{snp_end}'
        snp2peak[snp].add(peak)
    
    return snp2peak
    


# In[21]:


# first, get all TSS w/in 1Mb of the SNP. Note the distance to closest TSS
tss_near_snps = None
if USE_SIGNED_DISTANCE:
    tss_near_snps = get_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW).groupby(['snp', 'gene']).distance.agg(lambda x: x[x.abs()==x.abs().min()]).reset_index()
else:
    tss_near_snps = get_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW).groupby(['snp', 'gene']).distance.min().reset_index()

import random
tmp_file_name = 'tmp.{}.bed'.format(random.randrange(1, 100000000))
tss_near_snps_for_new_file = make_tss_file_of_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW)
tss_near_snps_for_new_file.to_csv(tmp_file_name, sep='\t', header=False, index=False)


results = []

for f in CICERO_FILES:
    cicero = pd.read_csv(f, sep='\t', header=None, names=['peak_1', 'peak_2', 'coaccess'])
    all_peaks = list(set(cicero.peak_1.unique().tolist() + cicero.peak_2.unique().tolist()))
    
    # find the peak that overlaps each SNP
    snps_to_peaks = cicero_assign_snps_to_peaks(SNPS, list(all_peaks))

    # keep all connections involving any of these peaks
    all_snp_containing_peaks = set()
    for v in snps_to_peaks.values():
        all_snp_containing_peaks = all_snp_containing_peaks.union(v)
    all_snp_containing_peaks = list(all_snp_containing_peaks)

    cicero = cicero[(cicero.peak_1.isin(all_snp_containing_peaks)) | (cicero.peak_2.isin(all_snp_containing_peaks))]

    # invert snp --> peaks, so we know which SNPs overlap each peak
    peak_to_snps = {peak: set() for peak in all_snp_containing_peaks}
    for snp, peaks in snps_to_peaks.items():
        for peak in peaks:
            peak_to_snps[peak].add(snp)
            
    # assign peaks to TSS
    all_peaks = list(set(cicero.peak_1.unique().tolist() + cicero.peak_2.unique().tolist()))
    peaks_to_tss = cicero_assign_peaks_to_tss(all_peaks, tmp_file_name, 5000, CHROM_SIZES, upstream_only=True)
    
    # now, for each SNP and each gene, get the max co-accessibility score
    df = []
    for i, r in cicero.iterrows():
        peak_1 = r['peak_1']
        peak_2 = r['peak_2']
        coaccess = r['coaccess']
        if coaccess == 'NA':
            continue
        if peak_1 in peak_to_snps:
            for snp in peak_to_snps[peak_1]:
                for tss in peaks_to_tss[peak_2]:
                    df.append([snp, tss, coaccess])
        if peak_2 in peak_to_snps:
            for snp in peak_to_snps[peak_2]:
                for tss in peaks_to_tss[peak_1]:
                    df.append([snp, tss, coaccess])
    df = pd.DataFrame(df, columns=['snp', 'gene', 'coaccess'])
    df.coaccess = df.coaccess.astype(float)
    df = df.groupby(['snp', 'gene']).coaccess.max().reset_index()

    df = tss_near_snps.merge(df, on=['snp', 'gene'], how='left')
    df['file'] = f
    results.append(df)
    
results = pd.concat(results)
os.remove(tmp_file_name)

# In[25]:


# plot
results.file = results.file.map(os.path.basename)
results['cluster'] = results.file.map(lambda x: x.split('.')[0])
results.head()


# In[93]:


for SNP in results.snp.unique():
    df = results[results.snp==SNP]
    df['label'] = df[['gene', 'distance']].apply(lambda x: '{} ({} kb)'.format(x[0], x[1]/1000), axis=1)
    label_order = df[['label', 'distance']].drop_duplicates().sort_values('distance').label.to_list()
    df = df.pivot(index='label', values='coaccess', columns='cluster')
    df = df.loc[label_order,:]
    df.columns = df.columns.to_series().map(cluster_to_cluster_name)
    COLOR_LIMIT = max([df.applymap(abs).max().max(), 0.2])
    fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(3, len(df)/3), gridspec_kw={'hspace': 0.02, 'height_ratios': [1, len(df)]})
    cbar_ax = axs[0]
    ax = axs[1]
    my_cm = cm.get_cmap('bwr', 256)
    my_cm.set_bad(color='grey')
    #norm = matplotlib.colors.Normalize(vmin=-1.,vmax=1.)
    fcbar = ax.imshow(df, aspect='auto', cmap=my_cm, vmin=-1*COLOR_LIMIT, vmax=COLOR_LIMIT)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df.index.to_list())
    ax.set_ylim(top=-0.5, bottom=ax.get_ylim()[0]+0.5)

    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns.to_list(), rotation=90)

    # add annotations
    for r_index in range(len(df)):
        for c_index in range(len(df.columns)):
            if not np.isnan(df.iat[r_index,c_index]):
                ax.text(x=c_index, y=r_index, s=round(df.iat[r_index,c_index], 2), va='center', ha='center')

    # add colorbar
    plt.colorbar(fcbar, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel('Cicero co-access.')
    cbar_ax.xaxis.tick_top()
    cbar_ax.xaxis.set_label_position('top')

    fig.savefig('{}.png'.format(SNP).replace(':', '_'), bbox_inches='tight')
    fig.clf()


#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


#TSS = '/lab/work/porchard/data/gencode-tss/hg19.bed'
#METHOD = 'spearman' # spearman or pearson
#LINK_FILES = glob.glob(f'/lab/work/porchard/sn-muscle-project/work/find-peak-gene-links-across-clusters-2/results/link-{METHOD}/*')
SNPS = ['chr5:53271420', 'chr12:26472562']
DISTANCE_TO_SHOW = 1e6
USE_SIGNED_DISTANCE = True

TSS = sys.argv[1]
LINK_FILES = sys.argv[2:]


# In[7]:


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


# In[4]:


# first, get all TSS w/in 1Mb of the SNP. Note the distance to closest TSS
tss_near_snps = None
if USE_SIGNED_DISTANCE:
    tss_near_snps = get_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW).groupby(['snp', 'gene']).distance.agg(lambda x: x[x.abs()==x.abs().min()]).reset_index()
else:
    tss_near_snps = get_tss_near_snps(SNPS, TSS, DISTANCE_TO_SHOW).groupby(['snp', 'gene']).distance.min().reset_index()

# In[5]:


links = pd.concat([pd.read_csv(f, sep='\t') for f in LINK_FILES])


# In[8]:


snps_to_peaks = cicero_assign_snps_to_peaks(SNPS, links.peak.map(lambda x: x.replace(':', '_')).to_list())


# In[15]:


# get the peaks that the SNPs lie in
for SNP in SNPS:
    if len(snps_to_peaks[SNP]) == 0:
        sys.stderr.write(f'Skipping SNP {SNP} (no peak found)\n')
        continue
    if len(snps_to_peaks[SNP]) > 1:
        sys.stderr.write(f'Skipping SNP {SNP} (multiple peaks found??)\n')
        continue
    PEAK = list(snps_to_peaks[SNP])[0].replace('_', ':')
    peak_links = links[links.peak==PEAK]
    # line plot
    for_plot = tss_near_snps.loc[tss_near_snps.snp==SNP,['gene', 'distance']].merge(peak_links, on='gene', how='left')
    for_plot = for_plot.sort_values('distance')
    for_plot['label'] = for_plot[['gene', 'distance']].apply(lambda x: '{} ({} kb)'.format(x[0], x[1]/1e3), axis=1)
    for_plot['label'] = pd.Categorical(for_plot['label'], categories=for_plot['label'], ordered=True)
    for_plot['idx'] = list(range(len(for_plot)))

    fig, ax = plt.subplots(figsize=(len(for_plot)/5, 5))

    YLIM = max([for_plot.zscore.abs().max()*1.1, 2])

    ax = sns.scatterplot(x='idx', y='zscore', data=for_plot)
    ax.set_xticks(for_plot.idx)
    ax.set_xticklabels(for_plot.label, rotation=90)
    ax.set_xlabel('Gene (TSS distance to SNP, kb)')
    ax.set_ylim(bottom=-1*YLIM, top=YLIM)
    ax.set_ylabel('Z-score for gene - peak correlation')
    ax.axhline(0, color='black', alpha=0.2, linestyle='dashed')
    
    fig.tight_layout()
    fig.savefig(f'{SNP}-peak-gene-links.png')
    fig.clf()

    df = for_plot.set_index('label').loc[:,['zscore']]
    fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(1, len(df)/3), gridspec_kw={'hspace': 0.02, 'height_ratios': [1, len(df)]})
    cbar_ax = axs[0]
    ax = axs[1]
    my_cm = cm.get_cmap('bwr', 256)
    my_cm.set_bad(color='grey')
    COLOR_LIM = max([for_plot.zscore.abs().max()*1, 2])
    fcbar = ax.imshow(df, aspect='auto', cmap=my_cm, vmin=-1*COLOR_LIM, vmax=COLOR_LIM)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df.index.to_list())
    ax.set_ylim(top=-0.5, bottom=ax.get_ylim()[0]+0.5)
    for i, r in for_plot.iterrows():
        if np.isnan(r['zscore']):
            continue
        ax.text(x=0, y=r['idx'], s=round(r['zscore'], 2), ha='center', va='center')
    ax.set_xticks([])
    ax.set_xlim(left=-0.5, right=0.5)
    
    # add colorbar
    plt.colorbar(fcbar, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel('Z-score for\npeak-gene link')
    cbar_ax.xaxis.tick_top()
    cbar_ax.xaxis.set_label_position('top')
    fig.savefig('{}-peak-gene-links-heatmap.png'.format(SNP).replace(':', '_'), bbox_inches='tight')
    fig.clf()


#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pygenometracks.tracks as pygtk
import argparse
import glob
import pandas as pd
import cmap_utils
import gzip
import pybedtools as bt


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument('--highlight', type=int, help='Highlight this position (int)')
parser.add_argument('--out', help='Name of plot.')
parser.add_argument('--gtf_bed', help='Path to bed file of genes')
parser.add_argument('--credible_set_file', help='Path to diamante credible set file')
parser.add_argument('--region', help='Region to plot, in format chr5:53,269,334-53,278,795')
parser.add_argument('--chromhmm', nargs='*', default=[], help='ChromHMM files')
parser.add_argument('--muscle_bws', nargs='*', default=[],help='Muscle BW files')
parser.add_argument('--other_bws', nargs='*', default=[], help='Other BW files')
args = parser.parse_args()

GENCODE_BED = args.gtf_bed
CREDIBLE_SET_FILE = args.credible_set_file
REGION = args.region
HIGHLIGHT_POS = args.highlight
CHROMHMM_FILES = args.chromhmm
BW_FILES = args.muscle_bws
OTHER_BW_FILES = args.other_bws
OUT = args.out


# In[ ]:


#GENCODE_BED = '/lab/work/porchard/sn-muscle-project/data/gencode-bed/gencode.hg19.bed'
#CREDIBLE_SET_FILE = '/lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ARL15_5_53271420.txt'
#REGION = 'chr5:53,269,334-53,278,795'
#HIGHLIGHT_POS = 53271420
#CHROMHMM_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/data/chromhmm/*.dense.bed')
#BW_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/process-by-cluster/results/bigwigs/atac/*hg19*')
#OTHER_BW_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/browser-sessions/2020-sn-muscle-GWAS-overlap/results/bigwigs/*')


BW_FILES = sorted(BW_FILES, key=lambda x: os.path.basename(x))
OTHER_BW_FILES = sorted(OTHER_BW_FILES, key=lambda x: os.path.basename(x))


# In[2]:


track_config = dict(file=GENCODE_BED, arrow_length=1000, arrow_interval=10, file_type='bed', height=3, color='black', merge_transcripts=False, labels=True, section_name='', prefered_name='gene_name', global_max_row=False, font_size=8, style="UCSC", all_labels_inside=False, labels_in_margin=True)
tk = pygtk.BedTrack(track_config)


# In[3]:


#track_config = dict(file='/lab/data/reference/human/hg19/annot/gencode.v19.annotation.gtf.gz', file_type='bed', height=3, color='black', line_width=0.1, merge_transcripts=False, labels=True, section_name='', prefered_name='gene_name', global_max_row=False, font_size=8, style="UCSC", all_labels_inside=False, labels_in_margin=True)
#tk = pygtk.BedTrack(track_config)


# In[30]:


CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t').astype(str)
cluster_colors = cmap_utils.make_colormap_dict(cluster_names.new_name.to_list())
cluster_names = dict(zip(cluster_names.old_name, cluster_names.new_name))


# In[51]:


# ARL15
#CREDIBLE_SET_FILE = '/lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ARL15_5_53271420.txt'
#REGION = 'chr5:53,269,334-53,278,795'
#HIGHLIGHT_POS = 53271420

# ITPR2
#CREDIBLE_SET_FILE = '/lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ITPR2_12_26453283.txt'
#REGION = 'chr12:26,469,532-26,475,063'
#REGION = 'chr12:26,436,640-26,491,955'
#HIGHLIGHT_POS = 26472562

region = REGION.replace(',', '')
chrom, start_end = region.split(':')
start, end = start_end.split('-')
start = int(start)
end = int(end)

credible_set = pd.read_csv(CREDIBLE_SET_FILE, sep='\t')
assert(HIGHLIGHT_POS in credible_set.Pos.to_list())


# In[53]:


def get_chromhmm_in_region(f, chrom, start, end):
    tmp = pd.read_csv(f, sep='\t', header=None).iloc[:,[0, 1, 2, 8]]
    tmp.columns = ['chrom', 'start', 'end', 'name']
    chromhmm_bt = bt.BedTool().from_dataframe(tmp)
    region_bt = bt.BedTool(f'{chrom} {start} {end}', from_string=True)
    return chromhmm_bt.intersect(region_bt, wa=True).to_dataframe()


# In[54]:


bw_files = pd.DataFrame({'file': OTHER_BW_FILES + BW_FILES})
bw_files['track_name'] = bw_files.file.map(lambda x: os.path.basename(x).replace('-hg19', '').replace('.bw', ''))
bw_files['track_name'] = bw_files['track_name'].map(lambda x: cluster_names[x] if x in cluster_names else x)
bw_files['track_color'] = bw_files['track_name'].map(lambda x: cluster_colors[x] if x in cluster_colors else (1,0.5,0,1))
bw_files = bw_files.set_index('file')


# In[60]:


rows = 3 + len(bw_files) + len(CHROMHMM_FILES)
height_ratios = [2, 1] + [2 for i in range(len(bw_files))] + [1 for i in CHROMHMM_FILES] + [4]

fig, axs = plt.subplots(nrows=rows, ncols=2, figsize=(10,3+rows/3), gridspec_kw={'width_ratios': [1, 8], 'height_ratios': height_ratios, 'hspace': 0.02, 'wspace': 0})
# plot PPA
ax = axs[0,0]
ax.text(1, 0.5, s='SNP PPA                 ', ha='right', va='center')
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, 1)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)


ax = axs[0,1]
ax.set_title(f'{REGION} (human)')
ax.set_xlim(start, end)
ax.set_ylim(0, credible_set.PPAg.max()*1.2)
ax.set_xticks([])
for x, y in zip(credible_set.Pos, credible_set.PPAg):
    ax.scatter(int(x), float(y), color='black', s=2)


# plot SNPs
ax = axs[1,0]
ax.patch.set_alpha(0)
ax.text(1, 0, s='T2D credible set SNPs   ', ha='right', va='bottom')
ax.set_xlim(0, 1)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)


ax = axs[1,1]
ax.set_xlim(start, end)
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.set_xticks([])
for pos in credible_set.Pos:
    ax.axvline(int(pos))

    
# fill in the bigwig tracks
for row_index, bw_file in enumerate(bw_files.index, 2):
    properties_dict = {'file': bw_file, 'height': 3, 'color': bw_files.at[bw_file,'track_color']}
    bw = pygtk.BigWigTrack(properties_dict)
    ax = axs[row_index,0]
    ax.text(1, 0.5, s=bw_files.at[bw_file,'track_name'] + ' ', ha='right', va='center')
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 1)
    ax = axs[row_index,1]
    bw.plot(ax, chrom, int(start), int(end))
    ax.xaxis.set_ticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylim(0, 1)
    ax.set_xlim(start, end)
    ax.yaxis.set_ticklabels([])
    ax.axvline(HIGHLIGHT_POS, color='red', linestyle='dotted')
    

# fill in the chromhmm tracks
for row_index, chromhmm_file in enumerate(CHROMHMM_FILES, 2+len(bw_files)):
    df = get_chromhmm_in_region(chromhmm_file, chrom, start, end)
    ax = axs[row_index,0]
    ax.text(1, 0.5, s=os.path.basename(chromhmm_file).replace('.dense.bed', ' chromHMM'), ha='right', va='center')
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 1)
    ax = axs[row_index,1]
    ax.set_ylim(0, 1)
    ax.set_xlim(start, end)
    ax.set_yticks([])
    ax.set_xticks([])
    for i, r in df.iterrows():
        c = [int(i)/255 for i in r['name'].split(',')]
        ax.fill_between(x=[r['start'], r['end']], y1=[1, 1], color=c)


# plot gene model
ax = axs[-1,0]
ax.remove()
ax = axs[-1,1]
ax.set_xlim(start, end)
ax.set_yticks([])
tk.plot(ax, chrom, int(start), int(end))

fig.savefig(OUT, bbox_inches='tight')


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[48]:


import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pygenometracks.tracks as pygtk
import argparse
import glob
import pandas as pd
import cmap_utils

parser = argparse.ArgumentParser()
parser.add_argument('--gene-bed', default=None, dest='gene_bed', help='Path to GTF/gene bed file. Used for gene models.')
parser.add_argument('--out', default='marker-genes-atac.png', help='Name of output plot.')
parser.add_argument('cluster_names', help='Path to file listing cluster_id, cluster_name (tsv)')
parser.add_argument('bigwig', nargs='+', help='Paths to bigwig files. Named {cluster_id}-.*.bw or aggregate.bw')
args = parser.parse_args()

def rgb_to_hex(rgb):
    x = None
    if isinstance(rgb, str):
        x = [int(i) for i in rgb.split(',')]
    else:
        x = rgb
    return '#%02x%02x%02x' % tuple(x)

# CONFIG
GENES = [('GAPDH', 'chr12:6,642,596-6,648,525'), 
         ('MYH1', 'chr17:10,414,886-10,423,629'),
         ('MYH7', 'chr14:23,878,228-23,912,613'), 
         #('CHRNA1', 'chr2:175,608,104-175,633,420'), 
         ('PDGFRA', 'chr4:55,091,467-55,101,711'), 
         #('FBN1', 'chr15:48,641,132-48,997,356'), 
         ('VWF', 'chr12:6,221,410-6,238,772'), 
         #('MYH11', 'chr16:15,948,630-15,953,189'),
         ('ACTA2', 'chr10:90,706,947-90,718,779'),
         #('CD163', 'chr12:7,650,408-7,657,741'),
         ('PTPRC', 'chr1:198,596,593-198,675,597'),
         ('PAX7', 'chr1:18,956,440-18,961,112')]
GENES = [(i[0], i[1].replace(',', '').replace('-', ':')) for i in GENES]
CLUSTER_NAMES = args.cluster_names #'/lab/work/porchard/sn-muscle-project/cluster-names.txt'
GENE_BED = args.gene_bed #'/lab/work/porchard/sn-muscle-project/data/gencode-coding/gencode.v19.annotation.gtf.gz'
BIGWIGS = args.bigwig # glob.glob('/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/process-by-cluster-round-1/atac-bigwigs/*hg19*') + ['/lab/work/porchard/sn-muscle-project/work/browser-sessions/2020-sn-muscle-cell-types-new-features/results/bigwigs/aggregate.bw']
cluster_to_bigwig = {os.path.basename(f).replace('.bw', '').split('-')[0]: f for f in BIGWIGS}


# In[49]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t').astype(str).rename(columns={'old_name': 'cluster_id', 'new_name': 'cluster_name'})
cluster_names['bigwig'] = cluster_names.cluster_id.map(cluster_to_bigwig)
cluster_to_color = cmap_utils.make_colormap_dict(cluster_names.cluster_id.to_list())
cluster_names['color'] = cluster_names.cluster_id.map(cluster_to_color)
if 'aggregate' in cluster_to_bigwig:
    cluster_names = pd.DataFrame({'cluster_id': ['aggregate'] + cluster_names.cluster_id.to_list(),
                                  'cluster_name': ['Aggregate'] + cluster_names.cluster_name.to_list(),
                                  'bigwig': [cluster_to_bigwig['aggregate']] + cluster_names.bigwig.to_list(),
                                  'color': [(1.0,0.5,0)] + cluster_names.color.to_list()
                                 })
print(cluster_names)


# In[29]:


tk = None
if args.gene_bed is not None:
    track_config = dict(file=GENE_BED, file_type='gtf', color='black', merge_transcripts=True, labels=False, display='interleaved',arrow_interval=999999, section_name='', prefered_name='gene_name', font_size=8, style="UCSC", all_labels_inside=True, labels_in_margin=False)
    tk = pygtk.BedTrack(track_config)


# In[50]:


# lay out figure
# columns: first column is the label, then one column per gene
# rows: first row is the gene name, second row in the gene model, then one row per cell type
HEIGHT_RATIOS = [0.5, 1] + [3 for i in range(len(cluster_names))]
fig, axs = plt.subplots(ncols=1+len(GENES), nrows=2+len(cluster_names), gridspec_kw={'hspace':0, 'wspace':0, 'height_ratios':HEIGHT_RATIOS}, figsize=(len(GENES)*1.2, len(cluster_names)*0.6))
axs[0,0].axis('off')
axs[1,0].axis('off')
for ax in axs[:,1]:
    ax.spines['left'].set_visible(False)
for ax in axs[1,:]:
    ax.spines['bottom'].set_visible(False)
for ax in axs[2,:]:
    ax.spines['top'].set_visible(False)


# fill in the gene names
for gene, ax in zip(GENES, axs[0:,1:].flatten()):
    gene_name = gene[0]
    ax.text(x=0.5, y=0, s=gene_name, ha='center', va='bottom', fontsize='large')
    for t in ax.get_yticklabels():
        t.set(visible=False)
    ax.axis('off')

# fill in the gene models
for gene, ax in zip(GENES, axs[1:,1:].flatten()):
    gene_name = gene[0]
    gene_location = gene[1]
    chrom, start, end = gene_location.split(':')
    if tk is not None:
        tk.plot(ax, chrom, int(start), int(end))
    for t in ax.get_yticklabels():
        t.set(visible=False)
    ax.axis('off')


# fill in cell type names
for cell_type_name, ax in zip(cluster_names.cluster_name, axs[2:,0].flatten()):
    ax.text(x=1, y=0.5, s=cell_type_name, ha='right', va='center', fontsize='large')
    ax.axis('off')

# fill in the bigwig tracks
row_index = 1
col_index = 2
for cell_type_index, cell_type in enumerate(cluster_names.index, 2):
    for gene_index, gene in enumerate(GENES, 1):
        properties_dict = {'file': cluster_names.at[cell_type,'bigwig'], 'height': 3, 'color': cluster_names.at[cell_type,'color']}
        bw = pygtk.BigWigTrack(properties_dict)
        gene_location = gene[1]
        chrom, start, end = gene_location.split(':')
        ax = axs[cell_type_index,gene_index]
        bw.plot(ax, chrom, int(start), int(end))
        ax.xaxis.set_ticklabels([])
        ax.set_yticks([], [])
        ax.set_xticks([], [])
        ax.set_ylim(0, 6)
        ax.yaxis.set_ticklabels([])
fig.tight_layout()
fig.savefig(args.out)
fig.clf()


# In[9]:





# In[ ]:





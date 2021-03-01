#!/usr/bin/env python
# coding: utf-8

# In[24]:


import os
import sys
import pandas as pd
import glob
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cmap_utils
import logomaker
import re
import math
import numpy as np


# In[46]:


def prob_matrix_to_ic_matrix(m):
    # Todo: verify correct
    # https://bioconductor.org/packages/release/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf
    # rows are A/C/G/T
    # columns as positions
    U = m.apply(lambda x: -1*(x*np.log2(x)).sum(skipna=True))
    fic = 2 - U
    return m * fic


# In[2]:


#MOTIF_DIR = '/lab/work/porchard/sn-muscle-project/work/gkmexplain-and-fimo-all-40000/work/b2/0e8c6de78fe671aa7fc025f5dc3623'
#GKMEXPLAIN_DIR = '/lab/work/porchard/sn-muscle-project/work/gkmexplain-and-fimo-all-40000/work/b2/0e8c6de78fe671aa7fc025f5dc3623'
#SEQUENCE_NAME = 'chr5:53271420:G:A' # e.g. 2:23934533:G:A
#PLAIN_MOTIF_DIR = '/lab/work/porchard/sn-muscle-project/work/gkmexplain-and-fimo-all-40000/work/b2/0e8c6de78fe671aa7fc025f5dc3623'
MOTIF_DIR, GKMEXPLAIN_DIR, SEQUENCE_NAME, PLAIN_MOTIF_DIR = sys.argv[1:]

# In[3]:


EXPLAINED_FILES = glob.glob(os.path.join(GKMEXPLAIN_DIR, '*reformatted.txt'))
MOTIF_FILES = glob.glob(os.path.join(MOTIF_DIR, '*.fimo.txt'))
SEQUENCE_NAME = SEQUENCE_NAME.replace('chr', '')


# In[4]:


def load_fimo_file(f):
    tmp = pd.read_csv(f, sep='\t')
    tmp = tmp.loc[tmp.motif_id != 'motif_id',['motif_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value', 'matched_sequence']].rename(columns={'p-value': 'pvalue'})
    tmp.start = tmp.start.astype(int)
    tmp.stop = tmp.stop.astype(int)
    tmp.score = tmp.score.astype(float)
    tmp['pvalue'] = tmp['pvalue'].astype(float)
    tmp['ref_or_alt'] = os.path.basename(f).replace('.fimo.txt', '')
    return tmp


# In[6]:


fimo = pd.concat([load_fimo_file(f) for f in MOTIF_FILES]).loc[:,['motif_id', 'sequence_name', 'start', 'stop', 'strand', 'pvalue', 'ref_or_alt']]
fimo = fimo.pivot(index=['motif_id', 'sequence_name', 'start', 'stop', 'strand'], values='pvalue', columns='ref_or_alt')#  tidyr::spread(key=ref_or_alt, value=p.value, fill=1) %>%
fimo = fimo[fimo.alt != fimo.ref].reset_index().fillna(1)
fimo.head()


# In[7]:


# for each of the motifs, plot their pwms against the motif discovered by gkmexplain
# load in the gkmexplain results

def read_explain_file(f):
    FILENAME_RE = '(cluster_.*)\.(ref|alt).explained.reformatted.txt'
    tmp = pd.read_csv(f, header=None, sep='\t', names=['snp', 'pos', 'nuc', 'score'])
    tmp['cluster'] = re.match(FILENAME_RE, os.path.basename(f)).group(1)
    tmp['ref_alt'] = re.match(FILENAME_RE, os.path.basename(f)).group(2)
    return tmp


# In[8]:


explained = pd.concat([read_explain_file(f) for f in EXPLAINED_FILES])
explained.snp = explained.snp.map(lambda x: x.replace('chr', ''))

fimo.sequence_name = fimo.sequence_name.map(lambda x: x.replace('chr', ''))
fimo = fimo[fimo.sequence_name.isin(explained.snp.to_list())]
fimo = fimo[fimo.sequence_name==SEQUENCE_NAME]


# In[9]:


explained.head()


# In[10]:


fimo.head()


# In[65]:


for CLUSTER in explained.cluster.unique():
    POSITIONS_SHOW = 50
    tmp = explained[(explained.cluster==CLUSTER) & (explained.snp==SEQUENCE_NAME)]
    # keep the middle 'positions' positions
    middle = math.ceil(max(tmp.pos) / 2)
    lower = middle - (0.5 * POSITIONS_SHOW)
    upper = middle + (0.5 * POSITIONS_SHOW)
    tmp = tmp[(tmp.pos >= lower) & (tmp.pos <= upper)]
    ref = tmp[tmp.ref_alt=='ref'].sort_values(['pos', 'nuc']).pivot(index='pos', columns='nuc', values='score')
    alt = tmp[tmp.ref_alt=='alt'].sort_values(['pos', 'nuc']).pivot(index='pos', columns='nuc', values='score')
    ref_minus_alt = ref - alt
    
    # add motif
    for i in fimo.index:
        MOTIF = fimo.at[i,'motif_id']
        FIMO_MOTIF_START = fimo.at[i,'start']
        STRAND = fimo.at[i,'strand']
        SNP = fimo.at[i,'sequence_name']

        # position 26 in gkmexplain should be position 51 in fimo...
        motif_file = os.path.join(PLAIN_MOTIF_DIR, MOTIF + '.txt').replace('::', '_')

        if not os.path.exists(motif_file):
            continue
        motif = pd.read_csv(motif_file, header=None, sep='\t')
        motif.columns = [i + FIMO_MOTIF_START - 1 for i in range(len(motif.columns))]
        motif.index = ['A', 'C', 'G', 'T']

        if motif.columns.max() > tmp.pos.max():
            continue
        if motif.columns.min() < tmp.pos.min():
            continue

        # TODO: not 100% sure this is correct
        if STRAND == '-':
            motif.columns = motif.columns.to_series().sort_values(ascending=False)
            motif = motif.loc[['T', 'G', 'C', 'A'],:]
            motif.index = ['A', 'C', 'G', 'T']
            motif = motif.loc[:,motif.columns.to_series().sort_values()]
    
        fig, axs = plt.subplots(nrows=4)

        ref_ax = axs[0]
        alt_ax = axs[1]
        delta_ax = axs[2]
        motif_ax = axs[3]
        y_min = min([
            ref.min().min(),
            alt.min().min(),
            ref_minus_alt.min().min()
        ])
        y_max = max([
            ref.max().max(),
            alt.max().max(),
            ref_minus_alt.max().max()
        ])

        ref_logo = logomaker.Logo(ref,
                                   ax=ref_ax,
                                   baseline_width=0,
                                   show_spines=True,
                                   vsep=0,
                                   width=.95)
        #ref_logo.highlight_position(50, color=(1, 0, 0, 0.5))
        alt_logo = logomaker.Logo(alt,
                                   ax=alt_ax,
                                   baseline_width=0,
                                   show_spines=True,
                                   vsep=0,
                                   width=.95)
        delta_logo = logomaker.Logo(ref_minus_alt,
                                   ax=delta_ax,
                                   baseline_width=0,
                                   show_spines=True,
                                   vsep=0,
                                   width=.95)


        for ax in [ref_ax, alt_ax, delta_ax]:
            ax.set_xticks([])
            ax.set_ylim(bottom=y_min, top=y_max)
        ref_ax.text(x=ref_ax.get_xlim()[1]*0.99, y=ref_ax.get_ylim()[1]*0.95, s='Ref.', va='top', ha='right')
        alt_ax.text(x=alt_ax.get_xlim()[1]*0.99, y=alt_ax.get_ylim()[1]*0.95, s='Alt.', va='top', ha='right')
        delta_ax.text(x=delta_ax.get_xlim()[1]*0.99, y=delta_ax.get_ylim()[1]*0.95, s='Ref. - Alt.', va='top', ha='right')
        fig.subplots_adjust(hspace=0)




        motif_logo = logomaker.Logo(prob_matrix_to_ic_matrix(motif).T,
                                   ax=motif_ax,
                                   baseline_width=0,
                                   show_spines=True,
                                   vsep=0,
                                   width=.95)
        motif_ax.set_xlim(left=ref_ax.get_xlim()[0], right=ref_ax.get_xlim()[1])
        motif_ax.text(x=motif_ax.get_xlim()[1]*0.99, y=motif_ax.get_ylim()[1]*0.95, s=MOTIF, va='top', ha='right')
        motif_ax.set_xticks([])
        fig.savefig(f'{SNP}___{CLUSTER}___{MOTIF}_{STRAND}_{FIMO_MOTIF_START}.pdf')
        fig.clf()


# In[ ]:





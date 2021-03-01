#!/usr/bin/env python
# coding: utf-8

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:
BARCODES, MTX, HORIZONTAL_LINE, VERTICAL_LINE, OUT = sys.argv[1:]

#BARCODES = '/lab/work/porchard/sn-muscle-project/work/rnaseq/results/starsolo/63_20_rna-hg19/Solo.out/GeneFull/raw/barcodes.tsv'
#MTX = '/lab/work/porchard/sn-muscle-project/work/rnaseq/results/starsolo/63_20_rna-hg19/Solo.out/GeneFull/raw/matrix.mtx'

barcodes = pd.read_csv(BARCODES, header=None, names=['barcode'])
barcodes['ind'] = range(1, len(barcodes)+1)

mtx = pd.read_csv(MTX, header=None, skiprows=3, sep=' ')
mtx.columns = ['feature', 'ind', 'count']

umis = mtx.groupby('ind')['count'].sum().reset_index()

both = umis.merge(barcodes, on='ind')

# plot
both = both.sort_values('count', ascending=False).set_index('barcode')
both['rank'] = range(1, len(both)+1)

fig, ax = plt.subplots()
ax.scatter(x='rank', y='count', data=both)
ax.set_xscale('log')
ax.set_yscale('log')
ax.axhline(int(HORIZONTAL_LINE), color='red', linestyle='dashed')
ax.axvline(int(VERTICAL_LINE), color='red', linestyle='dashed')
ax.grid(True)
ax.set_xlabel('Barcode rank')
ax.set_ylabel('# UMIs')
fig.tight_layout()
fig.savefig(OUT)
fig.clf()

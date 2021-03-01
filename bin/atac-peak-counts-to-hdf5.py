#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd
import numpy as np
from scipy.io import mmread

# convert ATAC peak counts to HDF5 format

#ATAC_COUNTS = '/lab/work/porchard/Nova-303/work/counts-in-peaks/results/1846_ATAC-hg19.counts.bed'
#OUT = 'rna.hdf5'
ATAC_COUNTS, INCLUDE_BARCODES, OUT = sys.argv[1:]


# In[ ]:


include_barcodes = pd.read_csv(INCLUDE_BARCODES, header=None).iloc[:,0]


# In[ ]:


atac_counts = pd.read_csv(ATAC_COUNTS, header=None, sep='\t', names=['chrom', 'start', 'end', 'barcode', 'count'])
atac_counts = atac_counts[atac_counts.barcode.isin(include_barcodes)]
atac_counts['feature'] = atac_counts.chrom + ':' + atac_counts.start.astype(str) + ':' + atac_counts.end.astype(str)
#atac_counts['fragments'] = atac_counts['count'].map(lambda x: round(x/2))
atac_counts['fragments'] = atac_counts['count']
atac_counts = atac_counts.pivot(index='barcode', columns='feature', values='fragments').fillna(0)
atac_counts = atac_counts.astype(int)
atac_counts.to_hdf(OUT, key=os.path.basename(OUT).replace('.hdf5', ''))


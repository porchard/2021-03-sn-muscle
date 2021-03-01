#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import sys
import pandas as pd


# In[7]:


#ATAC_THRESHOLDS = '/lab/work/porchard/sn-muscle-project/initial-thresholds-atac.txt'
ATAC_THRESHOLDS = sys.argv[1]

thresholds = pd.read_csv(ATAC_THRESHOLDS, sep='\t')
thresholds = thresholds[~thresholds.library.isin(['133152-hg19', '133154-hg19'])] # remove FANS libraries

LIBRARY_TO_NAME = {
    '125589-hg19': 'Human-rat mix (human nuclei) (ATAC)',
    '125589-rn6': 'Human-rat mix (rat nuclei) (ATAC)',
    '133151-hg19': 'no FANS rep1 (ATAC)',
    '133153-hg19': 'no FANS rep2 (ATAC)',
    '63_20-hg19': '20k nuclei (ATAC)',
    '63_40-hg19': '40k nuclei (ATAC)'
}

thresholds.library = thresholds.library.map(LIBRARY_TO_NAME)
thresholds.drop(columns=['max_hqaa']).rename(columns={'min_hqaa': 'min. reads after filtering', 'min_tss_enrichment': 'min. TSS enrichment', 'max_max_fraction_reads_from_single_autosome': 'max(max fraction reads from single autosome)'}).to_csv(sys.stdout, index=False)


# In[ ]:





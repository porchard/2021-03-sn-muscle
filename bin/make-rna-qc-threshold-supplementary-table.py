#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd


# In[4]:


#RNA_THRESHOLDS = '/lab/work/porchard/sn-muscle-project/initial-thresholds-rna.txt'
RNA_THRESHOLDS = sys.argv[1]

thresholds = pd.read_csv(RNA_THRESHOLDS, sep='\t')

LIBRARY_TO_NAME = {
    '133155-hg19': 'no FANS rep1 (RNA)',
    '133156-hg19': 'FANS rep1 (RNA)',
    '133157-hg19': 'no FANS rep2 (RNA)',
    '133158-hg19': 'FANS rep2 (RNA)',
    '63_20_rna-hg19': '20k nuclei (RNA)',
    '63_40_rna-hg19': '40k nuclei (RNA)'
}

thresholds.library = thresholds.library.map(LIBRARY_TO_NAME)
thresholds.rename(columns={'min_hqaa': 'min. UMIs', 'max_hqaa': 'max. UMIs', 'max_mitochondrial': 'max. fraction mitochondrial'}).to_csv(sys.stdout, index=False)


# In[ ]:





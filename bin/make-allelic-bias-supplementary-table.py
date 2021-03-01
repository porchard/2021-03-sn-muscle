#!/usr/bin/env python
# coding: utf-8

# In[15]:


import os
import sys
import pandas as pd
from statsmodels.stats.multitest import multipletests


# In[10]:


#ALLELIC_BIAS_TABLE = '/lab/work/porchard/sn-muscle-project/work/manuscript-figures-v2/results/allelic-bias/allelic-bias.txt'
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
ALLELIC_BIAS_TABLE, CLUSTER_NAMES = sys.argv[1:]


# In[12]:


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_id_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))


# In[16]:


bias = pd.read_csv(ALLELIC_BIAS_TABLE, sep='\t')
bias = bias[(bias.ref_counts>=1) & (bias.alt_counts>=1)]
bias['cluster_name'] = bias.cluster.map(cluster_id_to_cluster_name)
bias['qvalue'] = bias.groupby('cluster').pvalue.transform(lambda x: multipletests(x, method='fdr_bh')[1])


# In[17]:


bias = bias[['cluster_name', 'variant', 'ref', 'alt', 'ref_counts', 'alt_counts', 'fraction_ref', 'pvalue', 'qvalue']]
bias.to_csv(sys.stdout, sep='\t', index=False)


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[14]:


import os
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import matplotlib.ticker as ticker
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--individuals', help='Path to file of nucleus --> individual assignments.')
parser.add_argument('--rna', nargs='+', help='Path to RNA HDF5 files.')
parser.add_argument('--atac', nargs='+', help='Path to ATAC HDF5 files.')
args = parser.parse_args()

# ATAC
NUCLEUS_TO_INDIVIDUAL_ASSIGNMENTS = args.individuals #'/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files-with-dual-modality/data/nuclei-with-individual-assignments.txt'


# In[21]:


def is_dual_modality_atac(library):
    return ('1846' in library and 'ATAC' in library)


# In[15]:


individuals = pd.read_csv(NUCLEUS_TO_INDIVIDUAL_ASSIGNMENTS, sep='\t', header=None, names=['library', 'barcode', 'individual'])
individuals.loc[individuals.individual.map(lambda x: 'KSM' in x),'individual'] = 'KSM'
nucleus_to_individual = dict(zip(individuals.library + '-' + individuals.barcode, individuals.individual))


#RNA_GLOB = '/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files-with-dual-modality/results/rna-corrected/*hg19*'
#ATAC_GLOB = '/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files-with-dual-modality/results/atac/*hg19*'


# In[24]:


#rna_files = glob.glob(RNA_GLOB)
#atac_files = glob.glob(ATAC_GLOB)
#atac_files = [i for i in atac_files if not is_dual_modality_atac(i)]
rna_files = args.rna
atac_files = args.atac
atac_files = [i for i in atac_files if not is_dual_modality_atac(i)]


individual_matrices = {i: {'ATAC': [], 'RNA': []} for i in individuals.individual.unique()}


for f in rna_files:
    print(f)
    tmp = pd.read_hdf(f)
    # determine individuals present
    individuals_present = tmp.index.to_series().map(nucleus_to_individual).unique()
    for i in individuals_present:
        print(i)
        individual_matrices[i]['RNA'].append(tmp[tmp.index.to_series().map(nucleus_to_individual)==i])


for individual in individual_matrices:
    for modality in ['RNA']:
        if len(individual_matrices[individual][modality]) == 0:
            continue
        tmp = pd.concat(individual_matrices[individual][modality]).fillna(0, downcast='infer')
        name = f'{individual}_{modality}.hdf5'
        tmp.to_hdf(name, key=name.replace('.hdf5', ''))
for individual in list(individual_matrices.keys()):
    del individual_matrices[individual]['RNA']



for f in atac_files:
    print(f)
    tmp = pd.read_hdf(f)
    # determine individuals present
    individuals_present = tmp.index.to_series().map(nucleus_to_individual).unique()
    for i in individuals_present:
        print(i)
        individual_matrices[i]['ATAC'].append(tmp[tmp.index.to_series().map(nucleus_to_individual)==i])


# In[ ]:
for individual in individual_matrices:
    for modality in ['ATAC']:
        if len(individual_matrices[individual][modality]) == 0:
            continue
        tmp = pd.concat(individual_matrices[individual][modality]).fillna(0, downcast='infer')
        name = f'{individual}_{modality}.hdf5'
        tmp.to_hdf(name, key=name.replace('.hdf5', ''))

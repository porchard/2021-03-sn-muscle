#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import glob

#POS_FILE = '/lab/work/porchard/sn-muscle-project/work/get-diamante-cs-deltasvm/pos-file.txt'
#DSVM_FILE = '/lab/work/porchard/sn-muscle-project/work/delta-svm/all-keep-40k/results/cluster_0/snp_scores_cluster_0.txt'
POS_FILE, DSVM_FILE = sys.argv[1:]


def read_dsvm_file(f):
    tmp = pd.read_csv(f, sep='\t', header=None, names=['snp', 'score'], dtype={'snp': str, 'score': float})
    # convert to z-score
    avg = tmp.score.mean()
    sd = tmp.score.std()
    tmp['z-score'] = (tmp.score - avg) / sd
    return tmp

pos = pd.read_csv(POS_FILE, sep='\t', header=None, names=['chrom', 'pos'])
keep = (pos.chrom + ':' + pos.pos.astype(str)).unique()

dsvm = read_dsvm_file(DSVM_FILE)
dsvm['chrom'] = 'chr' + dsvm.snp.map(lambda x: x.split(':')[0])
dsvm['pos'] = dsvm.snp.map(lambda x: x.split(':')[1])
dsvm['chrom_and_pos'] = dsvm.chrom + ':' + dsvm.pos
dsvm.loc[dsvm.chrom_and_pos.isin(keep),['snp', 'score', 'z-score']].to_csv(sys.stdout, sep='\t', index=False)


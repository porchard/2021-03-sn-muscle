#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd

#MODEL_SCORE_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/work/delta-svm-hg19-all-gkm-40000/results/cluster_*/model_scores.tsv')
MODEL_SCORE_FILES = sys.argv[1:]

def load_scores(f):
    cluster = f.split('/')[-2]
    tmp = pd.read_csv(f, sep='\t', header=None, names=['stat', 'value'])
    tmp.stat = tmp.stat.map({'auc': 'AUC', 'accuracy': 'Accuracy', 'spec': 'Specificity', 'sens': 'Sensitivity'})
    tmp = tmp[tmp.stat!='Accuracy']
    tmp['cluster'] = cluster
    return tmp

scores = pd.concat([load_scores(f) for f in MODEL_SCORE_FILES]).pivot(index='cluster', columns='stat', values='value').sort_index().applymap(lambda x: round(x, 3)).reset_index()
scores.to_csv(sys.stdout, sep=',', index=False)

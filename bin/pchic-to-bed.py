#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd

#INTERACTIONS = '/lab/work/porchard/sn-muscle-project/data/pc-hiC/po-interactions.xlsv'
INTERACTIONS = sys.argv[1]

interactions = pd.read_excel(INTERACTIONS, skiprows=1)
interactions['chrom'] = interactions.Interacting_fragment.map(lambda x: x.split('.')[0])
interactions['start'] = interactions.Interacting_fragment.map(lambda x: x.split('.')[1])
interactions['end'] = interactions.Interacting_fragment.map(lambda x: x.split('.')[2])

interactions[['chrom', 'start', 'end', 'Promoter', 'Tissue_type']].to_csv(sys.stdout, sep='\t', header=None, index=False)


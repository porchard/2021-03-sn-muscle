#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import glob
import re

#BINS = '/lab/work/porchard/sn-muscle-project/data/hiC/bins.bed'
#FDR = 0.05
#FITHIC_FILES = glob.glob('/lab/work/porchard/sn-muscle-project/data/hiC/FitHiC_primary_cohort/FitHiC_output.OV_chr*.sparse.matrix.gz')

BINS = sys.argv[1]
FITHIC_FILES = sys.argv[2:]


bins = pd.read_csv(BINS, sep='\t', header=None, names=['chrom', 'start', 'end'])
bins.start = bins.start.astype(int)
bins.end = bins.end.astype(int)


def load_fithic_file(f):
    RE = 'FitHiC_output.(.*)_(.*).sparse.matrix.gz'
    fithic = pd.read_csv(f, sep='\t')
    tissue, chrom = re.match(RE, os.path.basename(f)).groups()
    return fithic.assign(chrom=chrom)


tmp = pd.concat([load_fithic_file(f) for f in FITHIC_FILES])


for CHROM in tmp.chrom.unique():
    chrom_df = tmp[tmp.chrom==CHROM]
    chrom_bins = bins[bins.chrom==CHROM].sort_values('start')
    chrom_bins['bin'] = range(1, len(chrom_bins)+1)

    if chrom_df.ColumnID.max() != chrom_bins.bin.max():
        sys.stderr.write('WARNING: Mismatch in bin numbers')
        chrom_df = chrom_df[chrom_df.RowID.isin(chrom_bins.bin)]
        chrom_df = chrom_df[chrom_df.ColumnID.isin(chrom_bins.bin)]
    chrom_df = chrom_df.merge(chrom_bins.rename(columns={'chrom': 'chrom1', 'start': 'start1', 'end': 'end1', 'bin': 'RowID'}), how='left')
    chrom_df = chrom_df.merge(chrom_bins.rename(columns={'chrom': 'chrom2', 'start': 'start2', 'end': 'end2', 'bin': 'ColumnID'}), how='left')
    chrom_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'QValue']].to_csv(sys.stdout, sep='\t', header=None, index=None)

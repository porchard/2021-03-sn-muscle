#!/usr/bin/env python

import os
import sys
import pybedtools as bt
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BEDS = sys.argv[1:]

N_FILES = len(BEDS)

beds = pd.concat([pd.read_csv(f, usecols=[0, 1, 2], delimiter='\t', header=None, names=['chrom', 'start', 'end']) for f in BEDS])
merged = bt.BedTool().from_dataframe(beds).sort().merge()
master_peak_df = merged.to_dataframe().assign(overlaps=0, peak=lambda df: df.chrom + ':' + df.start.astype(str) + ':' + df.end.astype(str)).set_index('peak')
sys.stderr.write('Found {} master peaks\n'.format(len(master_peak_df)))
sys.stderr.flush()

for f in BEDS:
    sys.stderr.write('Finding overlaps with {}\n'.format(f))
    sys.stderr.flush()
    overlapped = merged.intersect(f, u=True).to_dataframe().applymap(str)
    overlapped = overlapped.chrom + ':' + overlapped.start + ':' + overlapped.end
    master_peak_df.loc[overlapped,'overlaps'] = master_peak_df.loc[overlapped,'overlaps'] + 1

keep = master_peak_df[master_peak_df.overlaps>=float(N_FILES)/2]
keep[['chrom', 'start', 'end']].to_csv(sys.stdout, sep='\t', index=False, header=False)
master_peak_df.overlaps.hist(bins=N_FILES)
plt.xlabel('Number of experiments peak appears in')
plt.ylabel('Number of peaks')
plt.savefig('peak-sharing-distribution.pdf')

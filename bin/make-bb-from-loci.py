#!/usr/bin/env python
import os
import sys
import pandas as pd

credible_set_file = sys.argv[1]
COLORS = {0: '255,0,0', 1: '0,0,255'}

df = pd.read_csv(credible_set_file, delimiter='\t', header=None, names=['chrom', 'start', 'end', 'locus'])
df = df.sort_values(['chrom', 'start'])
locus_to_color = {locus: COLORS[index%2] for index, locus in enumerate(df.locus.unique(), 1)}
df['name'] = df.chrom + ':' + df.start.map(str) + ':' + df.end.map(str)
df['score'] = 1
df['strand'] = '.'
df['count'] = 1
df['sizes'] = df.end - df.start
df['starts'] = 0
df['color'] = df.locus.map(lambda x: locus_to_color[x])
df[['chrom', 'start', 'end', 'name', 'score', 'strand', 'start', 'end', 'color', 'count', 'sizes', 'starts']].to_csv(sys.stdout, sep='\t', header=None, index=False)

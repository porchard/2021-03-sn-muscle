#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pybedtools as bt
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Reduce a cicero co-accessibility file to a file containing only interactions in which at least one TSS is involved', add_help=True)
parser.add_argument('cicero', help='Cicero co-accessibility file')
parser.add_argument('tss', help='File of TSS (chrom, start, end, name, score, strand)')
parser.add_argument('snps', help='File of SNPs (chrom, start, end, name)')
parser.add_argument('chrom_sizes', help='Chrom size file')
parser.add_argument('--tss-distance', dest='tss_distance', type=int, help='Distance to TSS to consider a peak TSS-proximal. Default: 5000')
parser.add_argument('--tss-upstream', dest='tss_upstream', action='store_true', default=False, help='Only annotate a TSS to a peak if the peak is upstream of the TSS (correcting for strand). Default: False')
args = parser.parse_args()

CICERO_IN = args.cicero
TSS_FILE = args.tss
SNP_FILE = args.snps
CHROM_SIZES = args.chrom_sizes
TSS_DISTANCE = args.tss_distance


# first, determine all peaks present
KEEP_PEAKS = set()
ALL_PEAKS = set()

with open(CICERO_IN, 'r') as f:
    for line in f:
        first, second, score = line.rstrip().split()
        ALL_PEAKS.add(first)
        ALL_PEAKS.add(second)

peaks = pd.DataFrame([i.split('_') for i in ALL_PEAKS], columns=['chrom', 'start', 'end'])
peaks.start = peaks.start.astype(int)
peaks.end = peaks.end.astype(int)
peaks_bt = bt.BedTool().from_dataframe(peaks).sort()


tss_bt = bt.BedTool(TSS_FILE).sort()
peak2tss = dict() # chr_start_end --> [tss_1, tss_2, ...]

intersection = peaks_bt.intersect(tss_bt.slop(b=TSS_DISTANCE, g=CHROM_SIZES), wa=True, wb=True) if not args.tss_upstream else peaks_bt.intersect(tss_bt.slop(l=TSS_DISTANCE, r=0, s=True, g=CHROM_SIZES), wa=True, wb=True)
for i in intersection:
    peak_chrom, peak_start, peak_end, tss_chrom, tss_start, tss_end, gene, score, tss_strand = i
    peak = f'{peak_chrom}_{peak_start}_{peak_end}'
    if peak not in peak2tss:
        peak2tss[peak] = set()
    peak2tss[peak].add(gene)



snp = pd.read_csv(SNP_FILE, delimiter='\t', header=None)[[0, 1, 2, 3]]
snp.columns = ['chrom', 'start', 'end', 'name']
snp_bt = bt.BedTool().from_dataframe(snp).sort()
peak2snp = dict() # chr_start_end --> [snp_1, snp_2, ...]

for i in peaks_bt.intersect(snp_bt, wa=True, wb=True):
    peak_chrom, peak_start, peak_end, snp_chrom, snp_start, snp_end, snp = i
    peak = f'{peak_chrom}_{peak_start}_{peak_end}'
    if peak not in peak2snp:
        peak2snp[peak] = set()
    peak2snp[peak].add(snp)


print('\t'.join(['peak_1', 'peak_2', 'score', 'peak_1_tss', 'peak_2_tss', 'peak_1_snps', 'peak_2_snps']))
with open(CICERO_IN, 'r') as f:
    for line in f:
        peak_1, peak_2, score = line.strip().split()
        peak_1_tss = 'None' if peak_1 not in peak2tss else ','.join(peak2tss[peak_1])
        peak_2_tss = 'None' if peak_2 not in peak2tss else ','.join(peak2tss[peak_2])
        peak_1_snps = 'None' if peak_1 not in peak2snp else ','.join(peak2snp[peak_1])
        peak_2_snps = 'None' if peak_2 not in peak2snp else ','.join(peak2snp[peak_2])
        print('\t'.join([peak_1, peak_2, score, peak_1_tss, peak_2_tss, peak_1_snps, peak_2_snps]))


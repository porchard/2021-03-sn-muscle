#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import matplotlib.ticker as ticker
import seaborn as sns
import glob


@ticker.FuncFormatter
def read_count_formatter(x, pos):
    if x >= 1e9:
        return '{}B'.format(x/1e9)
    if x >= 1e6:
        return '{}M'.format(x/1e6)
    if x >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


RNA = '/lab/work/porchard/sn-muscle-project/work/determine-species-dual-modality/results/rna/summarized/rna-processed.txt'
ATAC = '/lab/work/porchard/sn-muscle-project/work/determine-species-dual-modality/results/atac/summarized/atac-processed.txt'
RNA_BARCODES = '~/github/snRNAseq-NextFlow/737K-arc-v1.txt'
ATAC_BARCODES = '~/github/snATACseq-NextFlow/737K-arc-v1.txt'
DROPLET_UTILS = '/lab/work/porchard/sn-muscle-project/work/test-droplet-utils/1846_RNA.test-vs-ambient.txt'
DROPLET_UTILS_FDR_THRESHOLD = 1e-5

droplet_utils = pd.read_csv(DROPLET_UTILS, sep='\t')
droplet_utils = droplet_utils[~droplet_utils['is.cell'].isnull()]
dropletutils_is_cell = {barcode: float(fdr)<= DROPLET_UTILS_FDR_THRESHOLD for barcode, fdr in zip(droplet_utils.barcode, droplet_utils.FDR)}

rna = pd.read_csv(RNA, sep='\t')
rna['fraction_hg19'] = rna['number_mapped_only_to_1846_RNA-hg19'] / (rna['number_mapped_only_to_1846_RNA-hg19']+rna['number_mapped_only_to_1846_RNA-rn6'])

STARSOLO_QC_FILES = glob.glob('/lab/work/porchard/Nova-315/work/starsolo-qc/results/qc/*')
starsolo_umi_counts = pd.concat([pd.read_csv(f, sep='\t').assign(file=f) for f in STARSOLO_QC_FILES])
starsolo_umi_counts = starsolo_umi_counts[['barcode', 'umis', 'fraction_mitochondrial', 'file']]
starsolo_umi_counts['species'] = starsolo_umi_counts.file.map(lambda x: x.replace('.qc.txt', '').split('-')[-1])
starsolo_umi_counts = starsolo_umi_counts[starsolo_umi_counts.umis!='None']
starsolo_umi_counts.umis = starsolo_umi_counts.umis.astype(int)
starsolo_umi_counts.fraction_mitochondrial = starsolo_umi_counts.fraction_mitochondrial.astype(float)
MIN_UMIS = 500
FRAC_MITO = 0.02
RNA_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES = starsolo_umi_counts.loc[(starsolo_umi_counts.umis>=MIN_UMIS) & (starsolo_umi_counts.fraction_mitochondrial<=FRAC_MITO),['barcode', 'species']]
RNA_MEET_QC_THRESHOLD = starsolo_umi_counts[(starsolo_umi_counts.umis>=MIN_UMIS) & (starsolo_umi_counts.fraction_mitochondrial<=FRAC_MITO)].barcode.to_list()


ATAQV_QC_FILES = glob.glob('/lab/work/porchard/Nova-303/work/qc/results/ataqv/*')
ataqv = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['nucleus', 'metric', 'value']) for f in ATAQV_QC_FILES])
ataqv = ataqv.pivot(index='nucleus', columns='metric', values='value')
ataqv['barcode'] = ataqv.index.map(lambda x: x.split('-')[-1])
ataqv['species'] = ataqv.index.map(lambda x: x.split('-')[1])
ataqv = ataqv[ataqv.tss_enrichment!='None']
ataqv.tss_enrichment = ataqv.tss_enrichment.astype(float)
ataqv.hqaa = ataqv.hqaa.astype(int)
ataqv.max_fraction_reads_from_single_autosome = ataqv.max_fraction_reads_from_single_autosome.astype(float)
MIN_HQAA = 20000
MIN_TSS_ENRICH = 2
MAX_MAX_FRAC_READS_FROM_AUTOSOME = 0.15
ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES = ataqv.loc[(ataqv.hqaa>=MIN_HQAA) & (ataqv.tss_enrichment>=MIN_TSS_ENRICH) & (ataqv.max_fraction_reads_from_single_autosome<=MAX_MAX_FRAC_READS_FROM_AUTOSOME),['barcode', 'species']]
ATAC_MEET_QC_THRESHOLD = ataqv[(ataqv.hqaa>=MIN_HQAA) & (ataqv.tss_enrichment>=MIN_TSS_ENRICH) & (ataqv.max_fraction_reads_from_single_autosome<=MAX_MAX_FRAC_READS_FROM_AUTOSOME)].barcode.to_list()
ataqv.head()

PASSED_QC = pd.read_csv('/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files/data/nuclei.txt', header=None, names=['library', 'barcode'], sep='\t')
PASSED_QC.head()
rna_pass_qc = PASSED_QC[PASSED_QC.library.isin(['1846_RNA-hg19', '1846_RNA-rn6'])].barcode.to_list()
atac_pass_qc = PASSED_QC[PASSED_QC.library.isin(['1846_ATAC-hg19', '1846_ATAC-rn6'])].barcode.to_list()
rat_rna_pass_qc = PASSED_QC[PASSED_QC.library=='1846_RNA-rn6'].barcode.to_list()
human_rna_pass_qc = PASSED_QC[PASSED_QC.library=='1846_RNA-hg19'].barcode.to_list()


atac = pd.read_csv(ATAC, sep='\t')
atac['fraction_hg19'] = atac['number_mapped_only_to_1846_ATAC-hg19'] / (atac['number_mapped_only_to_1846_ATAC-hg19']+atac['number_mapped_only_to_1846_ATAC-rn6'])

rna_barcodes = pd.read_csv(RNA_BARCODES, header=None)[0]
atac_barcodes = pd.read_csv(ATAC_BARCODES, header=None)[0]
atac_to_rna = dict(zip(atac_barcodes, rna_barcodes))
rna_to_atac = dict(zip(rna_barcodes, atac_barcodes))


atac_fraction = atac.loc[atac.total_reads>0,['barcode', 'fraction_hg19']].rename(columns={'fraction_hg19': 'atac_fraction_hg19'})
atac_fraction.barcode = atac_fraction.barcode.map(atac_to_rna)
fraction_atac_vs_fraction_rna = rna.loc[rna.total_reads>=0,['barcode', 'fraction_hg19']].rename(columns={'fraction_hg19': 'rna_fraction_hg19'}).merge(atac_fraction, on='barcode')
fraction_atac_vs_fraction_rna.head()


RNA_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES['pass_rna_qc'] = True
ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES = ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES.rename(columns={'barcode': 'atac_barcode'})
ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES['barcode'] = ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES.atac_barcode.map(atac_to_rna)
ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES['pass_atac_qc'] = True
MEET_QC_THRESHOLD_BOTH_MODALITIES = RNA_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES.merge(ATAC_MEET_QC_THRESHOLD_BARCODE_AND_SPECIES, how='outer', on=['species', 'barcode'])
MEET_QC_THRESHOLD_BOTH_MODALITIES = MEET_QC_THRESHOLD_BOTH_MODALITIES.dropna()
MEET_QC_THRESHOLD = MEET_QC_THRESHOLD_BOTH_MODALITIES.barcode.to_list()


fraction_atac_vs_fraction_rna['pass_threshold_in_at_least_one_species'] = fraction_atac_vs_fraction_rna.barcode.isin(MEET_QC_THRESHOLD)
fraction_atac_vs_fraction_rna['assigned_to_species'] = fraction_atac_vs_fraction_rna.barcode.isin(rna_pass_qc)
fraction_atac_vs_fraction_rna['assigned_to'] = 'Mixed'
fraction_atac_vs_fraction_rna.loc[fraction_atac_vs_fraction_rna.barcode.isin(rat_rna_pass_qc),'assigned_to'] = 'Rat'
fraction_atac_vs_fraction_rna.loc[fraction_atac_vs_fraction_rna.barcode.isin(human_rna_pass_qc),'assigned_to'] = 'Human'
fraction_atac_vs_fraction_rna['dropletutils_is_cell'] = fraction_atac_vs_fraction_rna.barcode.map(lambda x: False if x not in dropletutils_is_cell else dropletutils_is_cell[x])
fraction_atac_vs_fraction_rna['assigned_to_dropletutils'] = 'Mixed'
fraction_atac_vs_fraction_rna.loc[(fraction_atac_vs_fraction_rna.dropletutils_is_cell) & (fraction_atac_vs_fraction_rna.atac_fraction_hg19>=0.9) ,'assigned_to_dropletutils'] = 'Human'
fraction_atac_vs_fraction_rna.loc[(fraction_atac_vs_fraction_rna.dropletutils_is_cell) & (fraction_atac_vs_fraction_rna.atac_fraction_hg19<=0.1) ,'assigned_to_dropletutils'] = 'Rat'


fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(12,12))

# plot ATAC unique hg19 v unique rn6
axs[0,0].remove()
ax = axs[1,0]
ax.set_xlim(0, 5e5)
ax.set_ylim(0, 5e5)
ax.axvline(1e5, ymax=1e5/5e5, color='black', linestyle='dashed')
ax.axhline(1e5, xmax=1e5/5e5, color='black', linestyle='dashed')
sns.scatterplot(x='number_mapped_only_to_1846_ATAC-hg19', alpha=0.2, y='number_mapped_only_to_1846_ATAC-rn6', data=atac, ax=ax)
ax.set_xlabel('Fragments unique to human')
ax.set_ylabel('Fragments unique to rat')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_title('ATAC')

ax = axs[2,0]
ax.set_xlim(0, 1e5)
ax.set_ylim(0, 1e5)
sns.scatterplot(x='number_mapped_only_to_1846_ATAC-hg19', alpha=0.2, y='number_mapped_only_to_1846_ATAC-rn6', data=atac, ax=ax)
ax.set_xlabel('Fragments unique to human')
ax.set_ylabel('Fragments unique to rat')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_title('ATAC')

con = ConnectionPatch(xyA=(1e5,0), xyB=(1e5,1e5), axesA=axs[1,0], axesB=axs[2,0], coordsA="data", coordsB="data", color="black", linestyle='dotted')
axs[1,0].add_artist(con)
con = ConnectionPatch(xyA=(0,0), xyB=(0,1e5), axesA=axs[1,0], axesB=axs[2,0], coordsA="data", coordsB="data", color="black", linestyle='dotted')
axs[1,0].add_artist(con)



# plot RNA unique hg19 v unique rn6
axs[0,1].remove()
axs[1,1].remove()
ax = axs[2,1]
ax.set_xlim(0, 2e4)
ax.set_ylim(0, 2e4)
sns.scatterplot(x='number_mapped_only_to_1846_RNA-hg19', alpha=0.2, y='number_mapped_only_to_1846_RNA-rn6', data=rna, ax=ax)
ax.set_xlabel('UMIs unique to human')
ax.set_ylabel('UMIs unique to rat')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_title('RNA')

# plot fraction ATAC vs fraction RNA
# First, try DropletUtils
ax = axs[0,2]
tmp = fraction_atac_vs_fraction_rna[fraction_atac_vs_fraction_rna.dropletutils_is_cell]
assigned_to_counts = tmp.assigned_to_dropletutils.value_counts().reset_index()
assigned_to_counts['label'] = ['{} ({} %)'.format(species, round(100*frac/sum(assigned_to_counts.assigned_to_dropletutils), 1)) for species, frac in zip(assigned_to_counts['index'], assigned_to_counts['assigned_to_dropletutils'])]
update_labels = dict(zip(assigned_to_counts['index'], assigned_to_counts.label))
tmp.assigned_to_dropletutils = tmp.assigned_to_dropletutils.map(update_labels)
g = sns.scatterplot(x='atac_fraction_hg19', y='rna_fraction_hg19', alpha=0.2, hue='assigned_to_dropletutils', data=tmp, ax=ax)
g.legend().set_title('')
ax.set_ylabel('Fraction human unique RNA UMIs')
ax.set_xlabel('Fraction human unique ATAC fragments')
ax.set_title(f'Nuclei called using DropletUtils emptyDrops\n(FDR threshold = {100*DROPLET_UTILS_FDR_THRESHOLD}%)')
ax.axvline(0.9, linestyle='dashed', color='red')
ax.axvline(0.1, linestyle='dashed', color='red')


ax = axs[1,2]
tmp = fraction_atac_vs_fraction_rna[fraction_atac_vs_fraction_rna.pass_threshold_in_at_least_one_species]
assigned_to_counts = tmp.assigned_to.value_counts().reset_index()
assigned_to_counts['label'] = ['{} ({} %)'.format(species, round(100*frac/sum(assigned_to_counts.assigned_to), 1)) for species, frac in zip(assigned_to_counts['index'], assigned_to_counts['assigned_to'])]
update_labels = dict(zip(assigned_to_counts['index'], assigned_to_counts.label))
tmp.assigned_to = tmp.assigned_to.map(update_labels)
g = sns.scatterplot(x='atac_fraction_hg19', y='rna_fraction_hg19', alpha=0.2, hue='assigned_to', data=tmp, ax=ax)
g.legend().set_title('')
ax.set_ylabel('Fraction human unique RNA UMIs')
ax.set_xlabel('Fraction human unique ATAC fragments')
ax.set_title(f'Nuclei passing QC thresholds\nexcept species specificity')

ax = axs[2,2]
sns.scatterplot(x='atac_fraction_hg19', y='rna_fraction_hg19', alpha=0.2, data=fraction_atac_vs_fraction_rna[fraction_atac_vs_fraction_rna.barcode.isin(rna_pass_qc)], ax=ax)
ax.set_ylabel('Fraction human unique RNA UMIs')
ax.set_xlabel('Fraction human unique ATAC fragments')
ax.set_title('Pass QC nuclei')

fig.tight_layout()
fig.savefig('multiome-genome-comparisons.png', bbox_inches='tight')
fig.clf()

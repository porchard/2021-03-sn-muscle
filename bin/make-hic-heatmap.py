#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd
import pybedtools as bt
import matplotlib.pyplot as plt
import glob


# In[27]:


#HIC_ANNOTATED = glob.glob('/lab/work/porchard/sn-muscle-project/work/hic-overlap/results/annotate/*')

# get ALL tss within a certain distance of the snp
#SNP_BED = '/lab/work/porchard/sn-muscle-project/cicero-snps.bed'
#TSS = '/lab/work/porchard/data/gencode-tss/hg19.bed'

SNP_BED = sys.argv[1]
TSS = sys.argv[2]
HIC_ANNOTATED = sys.argv[3:]
DISTANCE_TO_SHOW = 1e6
USE_SIGNED_DISTANCE = True

SNP_TO_PEAK = {
    'rs7132434': 'ITPR2_peak',
    'rs702634': 'ARL15_peak',
    'rs227727': 'NOG_peak'
}

PEAK_TO_SNP = {v: k for k, v in SNP_TO_PEAK.items()}
abbreviation_to_name = {
    'H1': 'Embryonic stem cell',
    'MES': 'Mesendoderm',
    'MSC': 'Mesenchymal stem cell',
    'NPC': 'Neural progenitor cell',
    'TRO': 'Trophoblasts-like cell',
    'IMR90': 'Fibroblast (IMR90)',
    'GM12878': 'Lymphoblast (GM12878)',
    'CO': 'Prefrontal cortex',
    'HC': 'Hippocampus',
    'LG': 'Lung',
    'SB': 'Small bowel',
    'PA': 'Pancreas',
    'LI': 'Liver',
    'BL': 'Bladder',
    'OV': 'Ovary',
    'PO': 'Psoas muscle',
    'AO': 'Aorta',
    'LV': 'Left ventricle',
    'RV': 'Right ventricle',
    'SX': 'Spleen',
    'AD': 'Adrenal gland'
}

def experiment_to_tissue(f):
    return abbreviation_to_name[f.replace('.shmitt.annotated.bed', '').upper()]

# SNP_BED should be BED4, with the name in column 4
# TSS should be BED6
#SNP_BED, TSS = sys.argv[1:]


# In[3]:


HIC_ANNOTATED = [f for f in HIC_ANNOTATED if 'shmitt.annotated.bed' in f]
all_files = [os.path.basename(f) for f in HIC_ANNOTATED]


# In[5]:


def get_connections(f):
    hic = pd.read_csv(f, sep='\t')
    hic = hic[(hic.peak_1_tss!='None') | (hic.peak_2_tss!='None')] # get connections involving a TSS
    hic = hic[(hic.peak_1_snps!='None') | (hic.peak_2_snps!='None')] # get connections involving a peak
    df = []

    for i, r in hic.iterrows():
        if r['peak_1_snps'] != 'None':
            if r['peak_2_tss'] != 'None':
                snps = r['peak_1_snps'].split(',')
                tss = r['peak_2_tss'].split(',')
                for s in snps:
                    for t in tss:
                        df.append([s, t, f])
        if r['peak_2_snps'] != 'None':
            if r['peak_1_tss'] != 'None':
                snps = r['peak_2_snps'].split(',')
                tss = r['peak_1_tss'].split(',')
                for s in snps:
                    for t in tss:
                        df.append([s, t, f])

    df = pd.DataFrame(df, columns=['snp', 'tss', 'file'])
    
    return df


# In[6]:


connections = pd.concat([get_connections(f) for f in HIC_ANNOTATED])
connections['experiment'] = connections.file.map(os.path.basename)
connections = connections[['snp', 'tss', 'experiment']].drop_duplicates()


# In[10]:


nearby_tss = bt.BedTool(SNP_BED).sort().window(bt.BedTool(TSS).sort(), w=DISTANCE_TO_SHOW).to_dataframe()
nearby_tss['distance'] = -1*(nearby_tss.end - nearby_tss.strand)
if not USE_SIGNED_DISTANCE:
    nearby_tss['distance'] = (nearby_tss.end - nearby_tss.thickStart).abs()
nearby_tss = nearby_tss[['name', 'thickEnd', 'distance']].rename(columns={'thickEnd': 'gene', 'name': 'snp'})
if USE_SIGNED_DISTANCE:
    nearby_tss = nearby_tss.groupby(['snp', 'gene']).distance.agg(lambda x: x[x.abs()==x.abs().min()]).reset_index()
else:
    nearby_tss = nearby_tss.groupby(['snp', 'gene']).distance.min().reset_index()


# In[11]:


all_plot = pd.concat([nearby_tss.assign(experiment=e) for e in all_files])


# In[36]:


# for each experiment, note the TSS connections...
for SNP in connections.snp.unique():
    df = connections[connections.snp==SNP].rename(columns={'tss': 'gene'})
    df['connects'] = 1
    df = all_plot.loc[all_plot.snp==PEAK_TO_SNP[SNP],['gene', 'distance', 'experiment']].merge(df, how='left').fillna(0).sort_values('distance')
    df['label'] = [f'{i} ({j/1000} kb)' for i, j in zip(df.gene, df.distance)]
    hm = df[['label', 'experiment', 'connects']].pivot(index='label', columns='experiment', values='connects')
    hm = hm.loc[df.label.unique(),:]
    hm.columns = hm.columns.map(experiment_to_tissue)
    hm = hm.loc[:,list(sorted(hm.columns))]

    fig, ax = plt.subplots(figsize=(len(hm.columns)/4,1+len(hm)/4))
    ax.imshow(hm, aspect='auto', cmap='Reds')
    ax.set_yticks(range(len(hm)))
    ax.set_yticklabels(hm.index)
    ax.set_xticks(range(len(hm.columns)))
    ax.set_xticklabels(hm.columns, rotation=90)
    ax.set_ylim(bottom=ax.get_ylim()[0]+0.5, top=ax.get_ylim()[1]-0.5)
    fig.tight_layout()
    fig.savefig(f'{SNP}-connections.png')
    fig.clf()

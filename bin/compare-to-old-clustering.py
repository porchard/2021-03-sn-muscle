#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import cmap_utils
import alluvial

#OLD_CLUSTERING = '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/liger/round-1/second-louvain-clusters.txt'
#NEW_CLUSTERING = '/lab/work/porchard/daily/2021-01-22/labs/clusters-custom.txt'
#LIBRARY_LABELS = '/lab/work/porchard/sn-muscle-project/library-labels.txt'
#PREFIX = ''

OLD_CLUSTERING, NEW_CLUSTERING, LIBRARY_LABELS, PREFIX, OLD_LABEL, NEW_LABEL = sys.argv[1:]

cluster_names = pd.read_csv('/lab/work/porchard/sn-muscle-project/cluster-names.txt', sep='\t')
cluster_id_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))

library_labels = pd.read_csv(LIBRARY_LABELS, sep='\t')
library_to_modality = dict(zip(library_labels.library, library_labels.modality))


old = pd.read_csv(OLD_CLUSTERING, sep='\t', header=None, names=['library', 'barcode', 'old_cluster'])
new = pd.read_csv(NEW_CLUSTERING, sep='\s+', header=None, names=['library', 'barcode', 'new_cluster'])
new.library = new[['library', 'barcode']].apply(lambda x: x[0].replace(x[1], ''), axis=1)
both = old.merge(new, on=['library', 'barcode'])

assert((~both.old_cluster.isnull()).all())
assert((~both.new_cluster.isnull()).all())
print("Concordance: {}".format(sum(both.old_cluster==both.new_cluster) / len(both)))

MISSING_RNA = [i for i in new.library.unique() if i not in library_to_modality]
for i in MISSING_RNA:
    library_to_modality[i] = 'RNA'

new['species'] = new.library.map(lambda x: x.split('-')[1])
new['modality'] = new.library.map(library_to_modality)


tmp = both.groupby(['old_cluster', 'new_cluster']).size().reset_index().rename(columns={0: 'n'})
tmp['fraction'] = tmp.groupby('new_cluster').n.transform(lambda x: x/sum(x))

for_heatmap = tmp[['old_cluster', 'new_cluster', 'fraction']].pivot(index='old_cluster', columns='new_cluster', values='fraction').fillna(0)
for_heatmap = for_heatmap.sort_index()
for_heatmap = for_heatmap[list(sorted(for_heatmap.columns))]
concordance = for_heatmap.max().mean()

fig, ax = plt.subplots(figsize=(6,6))
ax.imshow(for_heatmap, vmin=0, vmax=1, cmap='Reds')
ax.set_ylabel(OLD_LABEL)
ax.set_xlabel(NEW_LABEL)
for i, r in tmp.iterrows():
    ax.text(y=r['old_cluster'], x=r['new_cluster'], s=int(r['n']), ha='center')
#ax.set_title(f'Concordance = {concordance}')
ax.set_xticks(range(len(for_heatmap.columns)))
ax.set_xticklabels(for_heatmap.columns.map(cluster_id_to_cluster_name), rotation=90)
ax.set_yticks(range(len(for_heatmap)))
ax.set_yticklabels(for_heatmap.index.map(cluster_id_to_cluster_name))
ax.set_ylim(bottom=len(for_heatmap)-0.5, top=-0.5)
fig.tight_layout()
fig.savefig(f'{PREFIX}.compare-to-old.png')
fig.clf()

# make alluvial plots
tmp = both.groupby(['old_cluster', 'new_cluster']).size().reset_index(name='n')

for_alluvial = {'{}. {}'.format(c, cluster_id_to_cluster_name[c]): {} for c in tmp.old_cluster.unique()}
for old_cluster, new_cluster, n in zip(tmp.old_cluster, tmp.new_cluster, tmp.n):
    for_alluvial['{}. {}'.format(old_cluster, cluster_id_to_cluster_name[old_cluster])]['{}. {} '.format(new_cluster, cluster_id_to_cluster_name[new_cluster])] = n
ax = alluvial.plot(for_alluvial, labels=(OLD_LABEL, NEW_LABEL))
fig = ax.get_figure()
fig.set_size_inches(5,7)
fig.savefig(f'{PREFIX}.alluvial-old-to-new.png', bbox_inches='tight')
fig.clf()

for_alluvial = {'{}. {}'.format(c, cluster_id_to_cluster_name[c]): {} for c in tmp.new_cluster.unique()}
for old_cluster, new_cluster, n in zip(tmp.old_cluster, tmp.new_cluster, tmp.n):
    for_alluvial['{}. {}'.format(new_cluster, cluster_id_to_cluster_name[new_cluster])]['{}. {} '.format(old_cluster, cluster_id_to_cluster_name[old_cluster])] = n
ax = alluvial.plot(for_alluvial, labels=(NEW_LABEL, OLD_LABEL))
fig = ax.get_figure()
fig.set_size_inches(5,7)
fig.savefig(f'{PREFIX}.alluvial-new-to-old.png', bbox_inches='tight')
fig.clf()



both['species'] = both.library.map(lambda x: x.split('-')[1])
both['modality'] = both.library.map(library_to_modality)

for s in both.species.unique():
    for m in both.modality.unique():
        tmp = both[(both.species==s) & (both.modality==m)].groupby(['old_cluster', 'new_cluster']).size().reset_index().rename(columns={0: 'n'})
        tmp['fraction'] = tmp.groupby('new_cluster').n.transform(lambda x: x/sum(x))

        for_heatmap = tmp[['old_cluster', 'new_cluster', 'fraction']].pivot(index='old_cluster', columns='new_cluster', values='fraction').fillna(0)
        for_heatmap = for_heatmap.sort_index()
        for_heatmap = for_heatmap[list(sorted(for_heatmap.columns))]
        concordance = for_heatmap.max().mean()

        fig, ax = plt.subplots(figsize=(5,5))
        ax.imshow(for_heatmap, vmin=0, vmax=1, cmap='Reds')
        ax.set_ylabel(OLD_LABEL)
        ax.set_xlabel(NEW_LABEL)
        for i, r in tmp.iterrows():
            ax.text(y=r['old_cluster'], x=r['new_cluster'], s=int(r['n']), ha='center')
        ax.set_title(f'{s}, {m}')
        fig.tight_layout()
        fig.savefig(f'{PREFIX}.compare-to-old-{m}-{s}.png')
        fig.clf()

        # alluvial plots
        for_alluvial = {'{} {}'.format(OLD_LABEL, c): {} for c in tmp.old_cluster.unique()}
        for old_cluster, new_cluster, n in zip(tmp.old_cluster, tmp.new_cluster, tmp.n):
            for_alluvial['{} {}'.format(OLD_LABEL, old_cluster)]['{} {}'.format(NEW_LABEL, new_cluster)] = n
        if len(for_alluvial) == 0:
            continue
        ax = alluvial.plot(for_alluvial)
        fig = ax.get_figure()
        fig.set_size_inches(7,7)
        fig.savefig(f'{PREFIX}.alluvial-old-to-new-{m}-{s}.png', bbox_inches='tight')
        fig.clf()
        
        for_alluvial = {'{} {}'.format(NEW_LABEL, c): {} for c in tmp.new_cluster.unique()}
        for old_cluster, new_cluster, n in zip(tmp.old_cluster, tmp.new_cluster, tmp.n):
            for_alluvial['{} {}'.format(NEW_LABEL, new_cluster)]['{} {}'.format(OLD_LABEL, old_cluster)] = n
        ax = alluvial.plot(for_alluvial)
        fig = ax.get_figure()
        fig.set_size_inches(7,7)
        fig.savefig(f'{PREFIX}.alluvial-new-to-old-{m}-{s}.png', bbox_inches='tight')
        fig.clf()

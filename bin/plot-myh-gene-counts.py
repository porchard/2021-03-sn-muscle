#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import glob
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#CLUSTERING = '/lab/work/porchard/sn-muscle-project/work/process-by-cluster-promising-2021-01-12/clusters.txt'
#CLUSTER_NAMES = '/lab/work/porchard/sn-muscle-project/cluster-names.txt'
#LIGER_INPUT_MATRICES = glob.glob('/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/data/*')

CLUSTERING = sys.argv[1]
CLUSTER_NAMES = sys.argv[2]
LIGER_INPUT_MATRICES = sys.argv[3:]


cluster_names = pd.read_csv(CLUSTER_NAMES, sep='\t')
cluster_to_cluster_name = dict(zip(cluster_names.old_name, cluster_names.new_name))

clusters = pd.read_csv(CLUSTERING, sep='\t', header=None, names=['library', 'barcode', 'cluster'])
clusters['library'] = clusters[['library', 'barcode']].apply(lambda x: x[0].replace(x[1], ''), axis=1)
clusters['nucleus'] = clusters.library + '-' + clusters.barcode
nucleus_to_cluster = dict(zip(clusters.nucleus, clusters.cluster))


def load_interesting_genes(f, interesting_genes=['MYH1', 'MYH2', 'MYH4', 'MYH7']):
    individual, modality = os.path.basename(f).replace('.hdf5', '').split('_')
    tmp = pd.read_hdf(f)
    tmp.columns = [i.upper() for i in tmp.columns]
    return tmp[interesting_genes].assign(individual=individual, modality=modality)
    


count = pd.concat([load_interesting_genes(f) for f in LIGER_INPUT_MATRICES])
count = count[count.index.to_series().isin(nucleus_to_cluster.keys())]

count['MYH1+2+4'] = count[['MYH1', 'MYH2', 'MYH4']].sum(axis=1)
count['cluster'] = count.index.to_series().map(nucleus_to_cluster)
count['individual_and_modality'] = count.individual + ', ' + count.modality

# plot:
# MYH7 vs MYH1+2+4
# MYH7 vs MYH1
# MYH7 vs MYH2
# MYH7 vs MYH4

# one column per cell type
# one row per individual/modality combo

def make_plot(count, x, y):
    cell_types = list(sorted(count.cluster.unique()))
    individual_and_modality = list(sorted(count.individual_and_modality.unique()))
    fig, axs = plt.subplots(nrows=len(individual_and_modality), ncols=len(cell_types), figsize=(3*len(cell_types), 3*len(individual_and_modality)))
    for i1, cell_type in enumerate(cell_types):
        for i2, i_and_m in enumerate(individual_and_modality):
            tmp = count[(count.cluster==cell_type) & (count.individual_and_modality==i_and_m)]
            ax = axs[i2,i1]
            ax.scatter(x=tmp[x], y=tmp[y], alpha=0.1, color='black')
            lim_max = max([ax.get_ylim()[1], ax.get_xlim()[1]])
            ax.set_ylim(bottom=0, top=lim_max*1.1)
            ax.set_xlim(left=0, right=lim_max*1.1)
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            ax.set_title('{}\n{}'.format(cluster_to_cluster_name[cell_type], i_and_m))
    fig.tight_layout()
    return fig, axs


for y in ['MYH1', 'MYH2', 'MYH4', 'MYH1+2+4']:
    # plot just the fiber types
    fig, ax = make_plot(count[count.cluster.isin([0, 1])], 'MYH7', y)
    FIGNAME = 'fiber-types-per-individual-MYH7-vs-{}.png'.format(y).replace('+', '_')
    fig.savefig(FIGNAME)
    fig.clf()
    # plot all clusters
    fig, ax = make_plot(count, 'MYH7', y)
    FIGNAME = 'all-cell-types-per-individual-MYH7-vs-{}.png'.format(y).replace('+', '_')
    fig.savefig(FIGNAME)
    fig.clf()

# plot per species
per_species = count.copy()
per_species['species'] = count.individual.map(lambda x: 'rat' if 'rat' in x else 'human')
per_species['individual_and_modality'] = per_species.species + ', ' + per_species.modality

for y in ['MYH1', 'MYH2', 'MYH4', 'MYH1+2+4']:
    # plot just the fiber types
    fig, ax = make_plot(per_species[per_species.cluster.isin([0, 1])], 'MYH7', y)
    FIGNAME = 'fiber-types-per-species-MYH7-vs-{}.png'.format(y).replace('+', '_')
    fig.savefig(FIGNAME)
    fig.clf()
    # plot all clusters
    fig, ax = make_plot(per_species, 'MYH7', y)
    FIGNAME = 'all-cell-types-per-species-MYH7-vs-{}.png'.format(y).replace('+', '_')
    fig.savefig(FIGNAME)
    fig.clf()

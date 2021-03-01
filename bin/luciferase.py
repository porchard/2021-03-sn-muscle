#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np


SPREADSHEET = '/lab/work/porchard/sn-muscle-project/data/luciferase/ARL15-rs702634_Luciferase_Rani_12-21-20.xlsx'
excel = pd.read_excel(SPREADSHEET)

COLS = {
    'MSCs_REP_1' : 'Unnamed: 1',
    'MSCs_REP_2' : 'Unnamed: 2',
    'Preadipocytes_REP_1' : 'Unnamed: 4',
    'Preadipocytes_REP_2' : 'Unnamed: 5',
    'Adipocytes_REP_1' : 'Unnamed: 7'
}
ROWS = {
    'G (forward)' : list(range(4, 4+4)),
    'A (forward)' : list(range(10, 10+5)),
    'G (reverse)' : list(range(17, 17+5)),
    'A (reverse)' : list(range(24, 24+5)),
    'EV' : list(range(31, 31+4))
}

df = []

for i, col in COLS.items():
    for j, row in ROWS.items():
        vals = excel.loc[row,col].to_list()
        for v in vals:
            df.append([i, j, v])

df = pd.DataFrame(df, columns=['cell_type_and_rep', 'allele_and_direction', 'relative_value'])
df['cell_type'] = df.cell_type_and_rep.map(lambda x: x.split('_')[0])
df['rep'] = df.cell_type_and_rep.map(lambda x: x.split('_')[-1]).astype(int)
df['allele_and_direction'] = pd.Categorical(df['allele_and_direction'], categories=['EV', 'G (forward)', 'A (forward)', 'G (reverse)', 'A (reverse)'], ordered=True)
df['allele'] = df.allele_and_direction.map(lambda x: x.split(' ')[0])


# plot out MSCs, combining rep 1 and rep 2
mscs = df.loc[df.cell_type=='MSCs',['allele_and_direction', 'relative_value', 'rep', 'allele']]

fig, ax = plt.subplots(figsize=(3, 4))
sns.boxplot(x='allele_and_direction', y='relative_value', hue='allele', dodge=False, palette={'G': '#E7A10F', 'A': '#2D7B34', 'EV': 'lightgrey'}, ax=ax, data=mscs)
g = sns.stripplot(x='allele_and_direction', y='relative_value', hue='rep', ax=ax, data=mscs)
handles, labels = ax.get_legend_handles_labels()
g.legend(handles[3:], labels[3:], loc='lower left', bbox_to_anchor=(0, 1), title='Replicate', ncol=2) #.remove()
ax.set_xlabel('Allele (orientation)')
ax.set_ylabel('Rel. luciferace expression')
for t in ax.get_xticklabels():
    t.set(rotation=90)

# add fold change + p-value for forward
g_forward = mscs[mscs.allele_and_direction=='G (forward)'].relative_value
a_forward = mscs[mscs.allele_and_direction=='A (forward)'].relative_value
p_forward = ttest_ind(g_forward, a_forward).pvalue
fold_change = round(np.mean(g_forward) / np.mean(a_forward), 1)
MAX_VALUE = mscs.relative_value.max()*1.1
ax.plot([1, 1, 2, 2], [MAX_VALUE*0.95, MAX_VALUE, MAX_VALUE, MAX_VALUE*0.95], color='black')
ax.text(1.5, MAX_VALUE, s='F.C.={}\np = {:.2e}'.format(fold_change, p_forward), va='bottom', ha='center')
ax.set_ylim(0, MAX_VALUE*1.3)

# now for reverse
g_reverse = mscs[mscs.allele_and_direction=='G (reverse)'].relative_value
a_reverse = mscs[mscs.allele_and_direction=='A (reverse)'].relative_value
p_reverse = ttest_ind(g_reverse, a_reverse).pvalue
fold_change = round(np.mean(g_reverse) / np.mean(a_reverse), 1)
#MAX_VALUE = mscs[mscs.allele_and_direction.map(lambda x: 'reverse' in x)].relative_value.max()*1.1
ax.plot([3, 3, 4, 4], [MAX_VALUE*0.75, MAX_VALUE*0.8, MAX_VALUE*0.8, MAX_VALUE*0.75], color='black')
ax.text(3.5, MAX_VALUE*0.8, s='F.C.={}\np = {:.2e}'.format(fold_change, p_reverse), va='bottom', ha='center')

fig.tight_layout()
fig.savefig('luciferase-msc-reps-combined.pdf')
fig.clf()


# In[83]:


# calculate p-value for G forward vs A forward:
# will do independent t-test
g_forward = mscs[mscs.allele_and_direction=='G (forward)'].relative_value
a_forward = mscs[mscs.allele_and_direction=='A (forward)'].relative_value
ttest_ind(g_forward, a_forward).pvalue


# In[87]:


# plot each rep and each cell type independently
cell_types = ['MSCs', 'Preadipocytes', 'Adipocytes']
reps = [1, 2]

fig, axs = plt.subplots(figsize=(len(reps)*3,len(cell_types)*3), nrows=len(cell_types), ncols=len(reps))
for row, cell_type in enumerate(cell_types):
    for col, rep in enumerate(reps):
        tmp = df[(df.cell_type==cell_type) & (df.rep==rep)]
        ax = axs[row,col]
        if len(tmp) == 0:
            ax.set_visible(False)
            continue
        sns.boxplot(x='allele_and_direction', y='relative_value', color='lightgrey', showfliers=False, ax=ax, data=tmp)
        sns.stripplot(x='allele_and_direction', y='relative_value', color='black', ax=ax, data=tmp)
        ax.set_xlabel('Allele (orientation)')
        ax.set_ylabel('Rel. luciferace expression')
        ax.set_title(f'{cell_type}, rep. {rep}')
        for t in ax.get_xticklabels():
            t.set(rotation=90)
        # add fold change + p-value for forward
        g_forward = tmp[tmp.allele_and_direction=='G (forward)'].relative_value
        a_forward = tmp[tmp.allele_and_direction=='A (forward)'].relative_value
        p_forward = ttest_ind(g_forward, a_forward).pvalue
        fold_change = round(np.mean(g_forward) / np.mean(a_forward), 1)
        MAX_VALUE = tmp.relative_value.max()*1.1
        ax.plot([1, 1, 2, 2], [MAX_VALUE*0.95, MAX_VALUE, MAX_VALUE, MAX_VALUE*0.95], color='black')
        ax.text(1.5, MAX_VALUE, s='F.C.={}\np = {:.2e}'.format(fold_change, p_forward), va='bottom', ha='center')
        ax.set_ylim(0, MAX_VALUE*1.3)

        # now for reverse
        g_reverse = tmp[tmp.allele_and_direction=='G (reverse)'].relative_value
        a_reverse = tmp[tmp.allele_and_direction=='A (reverse)'].relative_value
        p_reverse = ttest_ind(g_reverse, a_reverse).pvalue
        fold_change = round(np.mean(g_reverse) / np.mean(a_reverse), 1)
        ax.plot([3, 3, 4, 4], [MAX_VALUE*0.75, MAX_VALUE*0.8, MAX_VALUE*0.8, MAX_VALUE*0.75], color='black')
        ax.text(3.5, MAX_VALUE*0.8, s='F.C.={}\np = {:.2e}'.format(fold_change, p_reverse), va='bottom', ha='center')

        # get p-value
        g_forward = tmp[tmp.allele_and_direction=='G (forward)'].relative_value
        a_forward = tmp[tmp.allele_and_direction=='A (forward)'].relative_value
        p = ttest_ind(g_forward, a_forward).pvalue
        print('{}, {}, G vs A forward: {}'.format(cell_type, rep, p))
        g_reverse = tmp[tmp.allele_and_direction=='G (reverse)'].relative_value
        a_reverse = tmp[tmp.allele_and_direction=='A (reverse)'].relative_value
        p = ttest_ind(g_reverse, a_reverse).pvalue
        print('{}, {}, G vs A reverse: {}'.format(cell_type, rep, p))
fig.tight_layout()
fig.savefig('luciferase-all.png')
fig.clf()


# In[ ]:





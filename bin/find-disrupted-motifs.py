#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys

REF, ALT = sys.argv[1:]

def load_fimo_file(f):
    tmp = pd.read_csv(f, sep='\t')
    tmp = tmp[tmp.motif_id!='motif_id']
    tmp['motif_instance'] = tmp[['motif_id', 'sequence_name', 'start', 'stop', 'strand']].apply(lambda x: '__'.join(x), axis=1)
    tmp['p-value'] = tmp['p-value'].astype(float)
    tmp['score'] = tmp['score'].astype(float)
    return tmp

ref = load_fimo_file(REF).loc[:,['motif_instance', 'score', 'p-value']].rename(columns={'score': 'ref_score', 'p-value': 'ref_p'})
alt = load_fimo_file(ALT).loc[:,['motif_instance', 'score', 'p-value']].rename(columns={'score': 'alt_score', 'p-value': 'alt_p'})
both = ref.merge(alt, on='motif_instance', how='outer')
both['motif'] = both.motif_instance.map(lambda x: x.split('__')[0])
both['sequence'] = both.motif_instance.map(lambda x: x.split('__')[1])
both[(both.alt_score!=both.ref_score)].to_csv(sys.stdout, sep='\t')


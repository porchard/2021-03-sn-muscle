#!/usr/bin/env python

import os
import sys
import json

ATAC_BAM_FILES = sys.argv[1:]

LIBRARIES = dict()

for b in ATAC_BAM_FILES:
    library = os.path.basename(b).split('.')[0] 
    genome = library.split('-')[1]
    LIBRARIES[library] = {'pruned': b, 'genome': genome, 'modality': 'ATAC'}

print(json.dumps({'libraries': LIBRARIES}, sort_keys=True, indent=4))

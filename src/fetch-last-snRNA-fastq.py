import os
import sys
from john_utilities import symlink
import glob
import re

RUN, TARGET_DIR = sys.argv[1:]

def parse_file_name(f):
    RE = '/lab/data/seqcore/{}/63-NM/Sample_(\d+)-NM-\d+_([a-z])/.*(20|40)K_.*_R(\d+)_.*.fastq.gz'.format(RUN)
    sample, well, nuclei, read = re.match(RE, f).groups()
    return {'sample': sample, 'well': well, 'nuclei': nuclei, 'read': read}

files = [f for f in glob.glob('/lab/data/seqcore/{}/63-NM/*/*.fastq.gz'.format(RUN)) if '_I1_' not in f]

for f in files:
    info = parse_file_name(f)
    info['library'] = '{}_{}'.format(info['sample'], info['nuclei'])
    info['run'] = RUN
    new_file = os.path.join(TARGET_DIR, '{library}_rna_{well}.{read}.fastq.gz'.format(**info))
    symlink(f, new_file, abs = True)

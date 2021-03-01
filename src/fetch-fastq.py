import os
import sys
from john_utilities import symlink
import glob
import re

RUN, TARGET_DIR = sys.argv[1:]

def parse_file_name(f):
    RE = '/lab/data/seqcore/Run_{}/parker/Sample_(\d+)_([a-z])/.*_R(\d+)_.*.fastq.gz'.format(RUN)
    sample, well, read = re.match(RE, f).groups()
    return {'sample': sample, 'well': well, 'read': read}

files = [f for f in glob.glob('/lab/data/seqcore/Run_{}/parker/*/*.fastq.gz'.format(RUN)) if '_I1_' not in f]

for f in files:
    info = parse_file_name(f)
    info['library'] = info['sample']
    info['run'] = RUN
    new_file = os.path.join(TARGET_DIR, '{library}_{well}.{read}.fastq.gz'.format(**info))
    symlink(f, new_file, abs = True)

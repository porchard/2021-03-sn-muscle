#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='List features in a matrix (HDF5 format)', add_help=True)
parser.add_argument('mat', help='Matrix (samples x features), in HDF5 format.')
args = parser.parse_args()

logging.info('Reading {}'.format(args.mat))
mat = pd.read_hdf(args.mat)
logging.info('Finished reading {}'.format(args.mat))

print('\n'.join(mat.columns.to_list()))

logging.info('Done.')


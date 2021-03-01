#!/usr/bin/env python
# coding: utf-8

import os
import sys
import logging
import pybedtools as bt
import argparse

parser = argparse.ArgumentParser(description='Given a bed file, get the center of each element.', add_help=True)
parser.add_argument('bed', help='Bed file.')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

for i in bt.BedTool(args.bed):
    center = round((i.end + i.start) / 2)
    i.start = center
    i.end = center+1
    print(str(i).rstrip())

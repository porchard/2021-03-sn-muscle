#!/usr/bin/env python
# coding: utf-8

import os
import sys
import logging
import csv

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

connections, snp, FDR = sys.argv[1:]
FDR = float(FDR)
region = snp

# given the file of connections and a position of interest (chrom:pos), fetch all regions that interact with that position

# for testing
#connections = '/lab/work/porchard/sn-muscle-project/data/hiC/psoas-hic-connections.txt'
#region = 'chr17:54776955'
#FDR = 0.05

def parse_snp(s):
    """
    s should have format chrom:pos
    """
    chrom, pos = s.split(':')
    return {
        'chrom': chrom,
        'pos': int(pos)
    }


def parse_interacting_region(s):
    """
    s should have format chrom:start:end
    """
    tmp = s.split(':')
    if len(tmp) != 3:
        raise ValueError('{s} is not in format chrom:start:end'.format(**locals()))
    chrom, start, end = s.split(':')
    return {
        'chrom': chrom,
        'start': int(start),
        'end': int(end)
    }


def snp_in_interacting_region(interacting_region, snp):
    """
    interacting_region should have format chrom:start:end
    snp should have format chrom:pos
    """
    interacting_region = parse_interacting_region(interacting_region)
    snp = parse_snp(snp)
    return (snp['chrom'] == interacting_region['chrom']) and (snp['pos'] >= interacting_region['start']) and (snp['pos'] < interacting_region['end'])
    
    

with open(connections, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for line in reader:
        if float(line['QValue']) > FDR:
            continue
        if snp_in_interacting_region(line['bin_1'], region):
            print('{chrom}\t{start}\t{end}'.format(**parse_interacting_region(line['bin_2'])))
        elif snp_in_interacting_region(line['bin_2'], region):
            print('{chrom}\t{start}\t{end}'.format(**parse_interacting_region(line['bin_1'])))


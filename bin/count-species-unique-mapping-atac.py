#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import re
import pysam
import logging

parser = argparse.ArgumentParser()
parser.add_argument('whitelist', help = 'Bam file to adjust tags in')
parser.add_argument('bam_in_1', help = 'Bam file to adjust tags in')
parser.add_argument('bam_in_2', help = 'Bam file to adjust tags in')
parser.add_argument('genome_1', help = 'Bam file to adjust tags in')
parser.add_argument('genome_2', help = 'Bam file to adjust tags in')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

whitelist = set()
with open(args.whitelist, 'r') as f:
    for line in f:
        line = line.rstrip()
        whitelist.add(line)


# ignore PCR duplication, since a read might be marked as a duplicate for one species but not for the other but map in both cases

def get_read_outcomes(bam, whitelist):
    logging.info(f'Getting read outcomes for {bam}')
    read_name_to_outcome = dict()
    with pysam.AlignmentFile(bam, 'rb') as f:
        for read in f.fetch(until_eof=True):
            if not read.has_tag('CB'):
                continue
            read_name = read.query_name
            barcode = read.get_tag('CB')
            if read.is_secondary: # or read.is_supplementary:
                continue
            if not read.is_paired:
                continue
            if not read.is_read1:
                continue
            if barcode not in whitelist:
                continue
            if barcode not in read_name_to_outcome:
                read_name_to_outcome[barcode] = dict()
            if read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair:
                read_name_to_outcome[barcode][read_name] = 0
            else:
                read_name_to_outcome[barcode][read_name] = 1
    return read_name_to_outcome


read_outcomes_1 = get_read_outcomes(args.bam_in_1, whitelist)
read_outcomes_2 = get_read_outcomes(args.bam_in_2, whitelist)

all_barcodes = set(list(read_outcomes_1.keys()) + list(read_outcomes_2.keys()))

print('\t'.join(['barcode', 'total_reads', f'number_mapped_to_{args.genome_1}', f'number_mapped_to_{args.genome_2}', f'number_mapped_only_to_{args.genome_1}', f'number_mapped_only_to_{args.genome_2}']))
for barcode in all_barcodes:
    # for each read, determine the totanumber that mapped to genome_1, number that mapped to genome_2, number that mapped only to genome 1, number that mapped only to genome_2
    genome_1_outcomes = read_outcomes_1[barcode]
    genome_2_outcomes = read_outcomes_2[barcode]
    all_reads = set(list(genome_1_outcomes.keys()) + list(genome_2_outcomes.keys()))
    counters = [len(all_reads), 0, 0, 0, 0] # total, mapped to genome_1, mapped to genome 2, mapped only to genome 1, mapped only to genome 2
    for read in all_reads:
        genome_1_outcome = genome_1_outcomes[read]
        genome_2_outcome = genome_2_outcomes[read]
        if genome_1_outcome == 1:
            counters[1] += 1
            if genome_2_outcome == 0:
                counters[3] += 1
        if genome_2_outcome == 1:
            counters[2] += 1
            if genome_1_outcome == 0:
                counters[4] += 1
    print('\t'.join([barcode] + [str(i) for i in counters]))

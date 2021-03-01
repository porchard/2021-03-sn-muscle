#!/usr/bin/env python
# coding: utf-8

# In[52]:


import os
import sys
import pandas as pd
from ChunkedFileReader import ChunkedFileReader
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('bam_in_1', help = 'Bam file to adjust tags in')
parser.add_argument('bam_in_2', help = 'Bam file to adjust tags in')
parser.add_argument('genome_1', help = 'Bam file to adjust tags in')
parser.add_argument('genome_2', help = 'Bam file to adjust tags in')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

#HUMAN = '/lab/work/porchard/sn-muscle-project/work/determine-species-dual-modality/null/atac/sorted/1846_ATAC-hg19.sorted.txt'


flagstat = dict()

def _interpret_flag(d):
    binary = bin(d).replace('0b', '').rjust(12, '0')
    x = {interpretation: bit == '1' for interpretation, bit in zip(['is_paired', 'is_proper_pair', 'is_unmapped', 'mate_unmapped', 'is_reverse', 'mate_reverse', 'is_first', 'is_second', 'is_secondary', 'fail_qc', 'is_duplicate', 'is_supplementary'], binary[::-1])}
    return x


def interpret_flag(d):
    assert(isinstance(d, int))
    if d not in flagstat:
        flagstat[d] = _interpret_flag(d)
    return flagstat[d]


# In[49]:


# read in first nucleus
human = dict()
rat = dict()
print('\t'.join(['barcode', 'total_reads', f'number_mapped_to_{args.genome_1}', f'number_mapped_to_{args.genome_2}', f'number_mapped_only_to_{args.genome_1}', f'number_mapped_only_to_{args.genome_2}']))

with open(args.bam_in_1, 'r') as fh_human:
    with open(args.bam_in_2, 'r') as fh_rat:
        human_file_reader = ChunkedFileReader(fh_human, 1)
        rat_file_reader = ChunkedFileReader(fh_rat, 1)
        current_barcode = None
        for human in human_file_reader:
            human_read_dict = dict()
            rat_read_dict = dict()
            rat = next(rat_file_reader)
            duplicates = set()
            for read in human:
                read_name, barcode, flag = read.split()
                if current_barcode is None:
                    current_barcode = barcode
                    logging.info(f'Processing barcode {current_barcode}')
                assert(current_barcode == barcode)
                flags = interpret_flag(int(flag))
                if flags['is_secondary'] or not flags['is_paired'] or flags['is_supplementary'] or flags['is_second']:
                    continue
                if flags['is_duplicate']:
                    duplicates.add(read_name)
                human_read_dict[read_name] = flags['is_proper_pair'] and not flags['mate_unmapped'] and not flags['is_unmapped']
            for read in rat:
                read_name, barcode, flag = read.split()
                assert(current_barcode == barcode)
                flags = interpret_flag(int(flag))
                if flags['is_secondary'] or not flags['is_paired'] or flags['is_supplementary'] or flags['is_second']:
                    continue
                if flags['is_duplicate']:
                    duplicates.add(read_name)
                rat_read_dict[read_name] = flags['is_proper_pair'] and not flags['mate_unmapped'] and not flags['is_unmapped']
            all_reads = set(list(human_read_dict.keys()) + list(rat_read_dict.keys()))
            counters = [len(all_reads), 0, 0, 0, 0] # total, mapped to genome_1, mapped to genome 2, mapped only to genome 1, mapped only to genome 2
            for read in all_reads:
                if read in duplicates:
                    continue
                genome_1_outcome = human_read_dict[read]
                genome_2_outcome = rat_read_dict[read]
                if genome_1_outcome == True:
                    counters[1] += 1
                    if genome_2_outcome == False:
                        counters[3] += 1
                if genome_2_outcome == True:
                    counters[2] += 1
                    if genome_1_outcome == False:
                        counters[4] += 1
            print('\t'.join([current_barcode] + [str(i) for i in counters]))
            current_barcode = None


# In[ ]:





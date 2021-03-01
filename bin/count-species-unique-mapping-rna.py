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


def get_barcodes_in_file(f):
    barcodes = set()
    with open(f, 'r') as fh:
        for line in fh:
            barcode = line.rstrip().split()[1]
            barcodes.add(barcode)
    return barcodes


# find barcodes not shared between the files
BARCODES_IN_FILE_1 = get_barcodes_in_file(args.bam_in_1)
BARCODES_IN_FILE_2 = get_barcodes_in_file(args.bam_in_2)
PROCESS_BARCODES = BARCODES_IN_FILE_1.intersection(BARCODES_IN_FILE_2)

# read in first nucleus
print('\t'.join(['barcode', 'total_reads', f'number_mapped_to_{args.genome_1}', f'number_mapped_to_{args.genome_2}', f'number_mapped_only_to_{args.genome_1}', f'number_mapped_only_to_{args.genome_2}']))

with open(args.bam_in_1, 'r') as fh_human:
    with open(args.bam_in_2, 'r') as fh_rat:
        human_file_reader = ChunkedFileReader(fh_human, 1)
        rat_file_reader = ChunkedFileReader(fh_rat, 1)
        for human in human_file_reader:
            current_barcode = human[0].split()[1]
            if current_barcode not in PROCESS_BARCODES:
                logging.info(f'Skipping barcode {current_barcode} (only in human)')
                continue
            human_read_dict = dict()
            rat_read_dict = dict()
            rat = next(rat_file_reader)
            current_rat_barcode = rat[0].split()[1]
            while current_rat_barcode not in PROCESS_BARCODES:
                logging.info(f'Skipping barcode {current_rat_barcode} (only in rat)')
                rat = next(rat_file_reader)
                current_rat_barcode = rat[0].split()[1]
            logging.info(f'Processing barcode {current_barcode}')
            duplicates = set()
            bad_umis = set() # uncorrected umis that dont equal the corrected umi
            for read in human:
                read_name, barcode, uncorrected_umi, corrected_umi, flag = read.split()
                assert(current_barcode == barcode)
                flags = interpret_flag(int(flag))
                if flags['is_secondary'] or flags['is_supplementary']:
                    continue
                if uncorrected_umi != corrected_umi:
                    # umi was fixed in at least one species, so it was a duplicate
                    bad_umis.add(uncorrected_umi)
                if uncorrected_umi not in human_read_dict or not flags['is_unmapped']:
                    # record the UMI if it isn't already recorded, or if it's mapped
                    human_read_dict[uncorrected_umi] = not flags['is_unmapped']
            for read in rat:
                read_name, barcode, uncorrected_umi, corrected_umi, flag = read.split()
                assert(current_barcode == barcode)
                flags = interpret_flag(int(flag))
                if flags['is_secondary'] or flags['is_supplementary']:
                    continue
                if uncorrected_umi != corrected_umi:
                    # umi was fixed in at least one species, so it was a duplicate
                    bad_umis.add(uncorrected_umi)
                if uncorrected_umi not in rat_read_dict or not flags['is_unmapped']:
                    # record the UMI if it isn't already recorded, or if it's mapped
                    rat_read_dict[uncorrected_umi] = not flags['is_unmapped']
            all_umis = set(list(human_read_dict.keys()) + list(rat_read_dict.keys()))
            counters = [len(all_umis) - len(bad_umis), 0, 0, 0, 0] # total, mapped to genome_1, mapped to genome 2, mapped only to genome 1, mapped only to genome 2
            for umi in all_umis:
                if umi in bad_umis:
                    continue
                if umi not in human_read_dict or umi not in rat_read_dict:
                    continue
                genome_1_outcome = human_read_dict[umi]
                genome_2_outcome = rat_read_dict[umi]
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

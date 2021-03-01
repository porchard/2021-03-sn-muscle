import os
import sys
import gzip

gtf = sys.argv[1]

with gzip.open(gtf, 'rt') as f:
    for line in f:
        chrom, source, feature_type, start, end, score, strand, phase, info = line.rstrip().split('\t')
        if feature_type != 'transcript':
            continue
        info = [i.split(' ') for i in info.rstrip(';').split('; ')]
        info = {i[0]: i[1].replace('"', '') for i in info}
        tss_start = int(start) if strand == '+' else int(end) - 1
        tss_end = tss_start + 1
        gene_name = info['gene_name']
        gene_status = info['gene_status']
        gene_type = info['gene_type']
        if gene_status == 'NOVEL':
            continue
        if gene_type not in ['lncRNA', 'rRNA', 'protein_coding', 'retained_intron', 'processed_transcript', 'non_coding', 'ambiguous_orf', 'lincRNA', 'macro_lncRNA', 'bidirectional_promoter_lncRNA']:
            continue
        print('{chrom}\t{tss_start}\t{tss_end}\t{gene_name}\t.\t{strand}'.format(**locals()))



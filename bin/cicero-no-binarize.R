#!/usr/bin/env Rscript

library(cicero)

args <- commandArgs(T)
UMAP <- args[1]
COUNTS <- args[2]
CHROM_SIZE_FILE <- args[3]
WINDOW <- as.numeric(args[4])
OUT <- args[5]

# given the umap and the feature counts file, run cicero
counts <- read.table(COUNTS, head = F, as.is = T, sep = '\t', col.names = c('feature', 'nucleus', 'count'))

umap <- read.table(UMAP, head = F, as.is = T, sep = '\t', col.names = c('nucleus', 'dim1', 'dim2'))
rownames(umap) <- umap$nucleus
umap <- umap[rownames(umap) %in% counts$nucleus,c('dim1', 'dim2')]

input_cds <- make_atac_cds(COUNTS, binarize=F)
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap[,c('dim1', 'dim2')])

# run cicero
conns <- run_cicero(cicero_cds, genomic_coords = CHROM_SIZE_FILE, window=WINDOW)
write.table(conns, file = OUT, append = F, quote = F, sep = '\t', row.names = F, col.names = F)

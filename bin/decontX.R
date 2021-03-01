#!/usr/bin/env Rscript

library(celda)
library(dplyr)
library(ggplot2)
library(glue)


source('/home/porchard/src/load_hdf.R')

#COUNTS <- '/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files/results/rna/133155-hg19.hdf5'
#CLUSTERS <- '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/liger/round-1/second-louvain-clusters.txt'
#UMAP <- '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/liger/round-1/umap.txt'
#PREFIX <- '133155-hg19.'

args <- commandArgs(T)
COUNTS <- args[1]
CLUSTERS <- args[2]
UMAP <- args[3]
PREFIX <- args[4]


counts <- load_hdf(COUNTS)
nuclei <- sort(rownames(counts))
counts <- t(counts[nuclei,])
counts <- as(counts, "dgCMatrix")

clusters <- read.table(CLUSTERS, head=F, as.is=T, sep='\t', col.names=c('nucleus', 'cluster'))
rownames(clusters) <- clusters$nucleus
clusters <- clusters[nuclei,]
decontx_clusters <- clusters$cluster
names(decontx_clusters) <- clusters$nucleus

stopifnot(all(rownames(clusters) == colnames(counts)))

decontaminate <- decontX(counts, z=decontx_clusters)

umap <- read.table(UMAP, head=F, as.is=T, sep='\t', col.names=c('nucleus', 'dim1', 'dim2'))
rownames(umap) <- umap$nucleus
umap <- umap[nuclei,]

contamination <- as.data.frame(decontaminate$estimates$all_cells$contamination)
colnames(contamination) <- c('contamination')
contamination$nucleus <- rownames(contamination)
write.table(contamination[,c('nucleus', 'contamination')], file=glue('{PREFIX}contamination.txt'), append = F, quote=F, row.names=F, col.names=F, sep='\t')


umap <- left_join(umap, contamination)

p <- ggplot(umap) +
geom_point(aes(x=dim1, y=dim2, color=contamination)) +
scale_color_gradient2(low='grey', high='red')
png(glue('{PREFIX}contamination-on-umap.png'), height=5, width=6, units='in', res=300)
print(p)
dev.off()


MEDIAN <- median(contamination$contamination)
p <- ggplot(contamination) +
geom_histogram(aes(x=contamination)) +
geom_vline(xintercept = MEDIAN, color='red') +
xlab('DecontX-estimated contamination') +
ylab('Number of nuclei')
png(glue('{PREFIX}contamination-histogram.png'), height=5, width=6, units='in', res=300)
print(p)
dev.off()


new_counts <- as.data.frame(t(as.matrix(decontaminate$decontXcounts)))
new_counts_rounded <- round(new_counts)
new_counts$nucleus <- rownames(new_counts)
new_counts_rounded$nucleus <- rownames(new_counts)

write.table(new_counts, file=glue('{PREFIX}decontaminated-counts.txt'), append = F, quote=F, row.names=F, col.names=T, sep='\t')
write.table(new_counts_rounded, file=glue('{PREFIX}decontaminated-counts-rounded.txt'), append = F, quote=F, row.names=F, col.names=T, sep='\t')

#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("--rds"), action = "store", type = "character", default = "", help = "[Required] Path to the factorization output")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(liger)
library(glue)


RDS <- opts$rds # '/lab/work/porchard/sn-muscle-project/work/liger-human/results/factorize/30/KSM2_RNA/liger-factorized.Rda'


load(RDS)
sample_specificity <- calcDatasetSpecificity(int.muscle, dataset1 = 'KSM1_RNA', dataset2 = 'KSM2_RNA')
modality_specificity_ksm1 <- calcDatasetSpecificity(int.muscle, dataset1 = 'KSM1_RNA', dataset2 = 'KSM1_ATAC')
modality_specificity_ksm2 <- calcDatasetSpecificity(int.muscle, dataset1 = 'KSM2_RNA', dataset2 = 'KSM2_ATAC')


SAMPLE_SPECIFIC_FACTORS <- seq(1, length(sample_specificity[[3]]))[abs(sample_specificity[[3]]) >= 50]
MODALITY_SPECIFIC_FACTORS_KSM1 <- seq(1, length(modality_specificity_ksm1[[3]]))[abs(modality_specificity_ksm1[[3]]) >= 40]
MODALITY_SPECIFIC_FACTORS_KSM2 <- seq(1, length(modality_specificity_ksm2[[3]]))[abs(modality_specificity_ksm2[[3]]) >= 40]
MODALITY_SPECIFIC_FACTORS <- unique(c(MODALITY_SPECIFIC_FACTORS_KSM1, MODALITY_SPECIFIC_FACTORS_KSM2))
DROP_FACTORS <- c(MODALITY_SPECIFIC_FACTORS)


markers <- getFactorMarkers(int.muscle, factor.share.thresh=10000, num.genes = 10)
shared <- markers$shared
shared$is_ribosomal <- F
RIBOSOMAL_FACTORS <- shared %>% dplyr::group_by(factor_num) %>% dplyr::summarize(number_ribosomal=sum
(is_ribosomal)) %>% dplyr::filter(number_ribosomal>0) %>% dplyr::pull(factor_num) %>% unique()
KEEP_FACTORS <- 1:max(shared$factor_num)
DROP_FACTORS <- c(DROP_FACTORS, RIBOSOMAL_FACTORS)
KEEP_FACTORS <- KEEP_FACTORS[!KEEP_FACTORS %in% DROP_FACTORS]
print('Using factors:')
print(KEEP_FACTORS)
shared$keep_factor <- shared$factor_num %in% KEEP_FACTORS


int.muscle <- suppressWarnings(quantile_norm(int.muscle, dims.use = KEEP_FACTORS, knn_k=5, min_cells=5, ref_dataset = 'KSM1_RNA'))
int.muscle <- runUMAP(int.muscle, distance='cosine', dims.use = KEEP_FACTORS)

umap <- as.data.frame(int.muscle@tsne.coords)
colnames(umap) <- c('dim_1', 'dim_2')
umap$nucleus <- rownames(umap)
umap$library <- gsub('(.*)-(.*)', '\\1', umap$nucleus)
umap$barcode <- gsub('(.*)-(.*)', '\\2', umap$nucleus)
umap <- dplyr::select(umap, library, barcode, dim_1, dim_2)
write.table(umap, file = 'umap.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)

# plot by dataset
int.muscle <- louvainCluster(int.muscle, resolution = 0.1)
p <- plotByDatasetAndCluster(int.muscle, do.shuffle = T, return.plots = T)
png('plot-by-dataset.png', width = 6, height = 6, units='in', res=300)
print(p[[1]])
dev.off()

pdf(glue('plot-gene-loadings.pdf'))
plotGeneLoadings(int.muscle, return.plots = F, factor.share.thresh=10000, log.fc.thresh=0, dataset1 = 'KSM1_RNA', dataset2 = 'rat1_RNA')
dev.off()

pdf(glue('plot-factors.pdf'))
plotFactors(int.muscle, num.genes = 10)
dev.off()

pdf(glue('plot-cluster-proportions.pdf'))
plotClusterProportions(int.muscle)
dev.off()

pdf(glue('plot-cluster-factors.pdf'))
plotClusterFactors(int.muscle, use.aligned=T)
dev.off()


RESOLUTIONS <- seq(0.01, 0.2, by=0.01)
for(r in RESOLUTIONS) {
    int.muscle <- louvainCluster(int.muscle, resolution = r)
    # plot by cluster
    p <- plotByDatasetAndCluster(int.muscle, do.shuffle = T, return.plots = T)
    png(glue('plot-by-cluster-resolution-{r}.png'), width = 6, height = 6, units='in', res=300)
    print(p[[2]])
    dev.off()
    clusters <- as.data.frame(int.muscle@clusters)
    colnames(clusters) <- c('cluster')
    clusters$nucleus <- rownames(clusters)
    clusters$barcode <- gsub('.*-', '', clusters$nucleus)
    clusters$library <- gsub('(.*)-', '\\1', clusters$nucleus)
    num_factors <- max(shared$factor_num)
    write.table(clusters[,c('library', 'barcode', 'cluster')], glue('clusters-{num_factors}-{r}.txt'), append = F, quote = F, row.names = F, col.names = F)
}

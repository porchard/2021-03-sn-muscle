#!/usr/bin/env Rscript
# coding: utf-8

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(glue)
library(rhdf5)

load_hdf <- function(HDF) {
          tmp <- h5ls(HDF)
  GROUP <- paste("/", tmp$name[tmp$group=='/'], sep='')
    df <- h5read(HDF, GROUP, native = F, compoundAsDataFrame=T)
    counts <- t(df$block0_values)
      rownames(counts) <- as.character(df$axis1)
      colnames(counts) <- toupper(as.character(df$axis0))
        return(t(counts))
}

#count_files <- Sys.glob('/lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-5/data/*')
args <- commandArgs(T)
count_files <- args

liger_in <- lapply(count_files, load_hdf)
names(liger_in) <- sapply(count_files, function(f){gsub('.hdf5', '', basename(f))})

# integrate RNA
# then integrate ATAC, using PCA on variable genes
# then project labels from RNA --> ATAC
# then co-embed ATAC and RNA

# load RNA
rna_objects <- list()
for (i in grep('RNA', names(liger_in), value=T)) {
    counts <- liger_in[[i]]
    individual <- strsplit(i, '_')[[1]][1]
    modality <- strsplit(i, '_')[[1]][2]
    print(i)
    metadata <- data.frame(nucleus=colnames(counts))
    rownames(metadata) <- metadata$nucleus
    metadata$individual <- individual
    metadata$modality <- modality
    rna <- CreateSeuratObject(counts = counts, min.cells=5, assay = modality, project = i, metadata=metadata)
    rna[['ACTIVITY']] <- CreateAssayObject(counts = counts)
    rna$tech <- modality
    rna$individual <- individual
    rna$modality <- modality
    DefaultAssay(rna) <- 'ACTIVITY'
    rna_objects[[length(rna_objects)+1]] <- rna
}

# load ATAC
atac_objects <- list()
for (i in grep('ATAC', names(liger_in), value=T)) {
    counts <- liger_in[[i]]
    individual <- strsplit(i, '_')[[1]][1]
    modality <- strsplit(i, '_')[[1]][2]
    print(i)
    metadata <- data.frame(nucleus=colnames(counts))
    rownames(metadata) <- metadata$nucleus
    metadata$individual <- individual
    metadata$modality <- modality
    rna <- CreateSeuratObject(counts = counts, min.cells=5, assay = modality, project = i, metadata=metadata)
    rna[['ACTIVITY']] <- CreateAssayObject(counts = counts)
    rna$tech <- modality
    rna$individual <- individual
    rna$modality <- modality
    DefaultAssay(rna) <- 'ACTIVITY'
    atac_objects[[length(atac_objects)+1]] <- rna
}


# integrate RNA
print('Integrating RNA')
for (i in 1:length(rna_objects)) {
    print(i)
    rna_objects[[i]] <- NormalizeData(rna_objects[[i]], verbose = FALSE)
    rna_objects[[i]] <- FindVariableFeatures(rna_objects[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}


rna.anchors <- FindIntegrationAnchors(object.list = rna_objects, dims = 1:20)
rna.integrated <- IntegrateData(anchorset = rna.anchors, dims = 1:20)
DefaultAssay(rna.integrated) <- "integrated"

rna_variable_features <- VariableFeatures(rna.integrated)

# Run the standard workflow for visualization and clustering
rna.integrated <- ScaleData(rna.integrated, verbose = FALSE)
rna.integrated <- RunPCA(rna.integrated, npcs = 50, verbose = FALSE)

pdf('elbow-plot.pdf', height=5, width=5)
ElbowPlot(rna.integrated, ndims = 100)
dev.off()

DIMS <- 10
PCS <- 1:DIMS
DROP_PCS <- c(9, 11, 13, 15)
PCS <- PCS[!PCS %in% DROP_PCS]

rna.integrated <- RunUMAP(rna.integrated, n.neighbors = 30, reduction = "pca", dims = PCS)
rna.integrated <- FindNeighbors(rna.integrated, k.param=20, dims = PCS)
rna.integrated <- FindClusters(rna.integrated, algorithm = 1, resolution = 0.05)

p1 <- DimPlot(rna.integrated, reduction = "umap", group.by = "individual")
pdf('rna-umap-by-individual.pdf', height=5, width=5)
print(p1)
dev.off()

pdf('rna-umap-by-cluster.pdf', height=5, width=5)
p2 <- DimPlot(rna.integrated, reduction = "umap")
print(p2)
dev.off()


PLOT_GENES <- c('MYH1', "CHRNA1", "MYH7", 'MYH2', 'MYH4', 'PAX7', 'PDGFRB', 'MYH11', 'ACTA2', 'CD163', 'VWF')
png('integrated-marker-genes-rna.png', height=2.5*length(PLOT_GENES), width=2.5*3, units='in', res=300)
FeaturePlot(rna.integrated, c('MYH1', 'PAX7', 'MYH7'), pt.size = 0.2, split.by = 'individual')
dev.off()


# integrate ATAC
print('Integrating ATAC')
for (i in 1:length(atac_objects)) {
    print(i)
    atac_objects[[i]] <- NormalizeData(atac_objects[[i]], verbose = FALSE)
    VariableFeatures(atac_objects[[i]]) <- rna_variable_features
}


atac.anchors <- FindIntegrationAnchors(object.list = atac_objects, dims = 2:20)
atac.integrated <- IntegrateData(anchorset = atac.anchors, dims = 2:20)
DefaultAssay(atac.integrated) <- "integrated"

# visualization
atac.integrated <- ScaleData(atac.integrated, verbose = FALSE)
atac.integrated <- RunPCA(atac.integrated, npcs = 20, verbose = FALSE)
atac.integrated <- RunUMAP(atac.integrated, reduction = "pca", dims = 2:20)


PLOT_GENES <- c('MYH1', "CHRNA1", "MYH7", 'MYH2', 'MYH4', 'PAX7', 'PDGFRB', 'MYH11', 'ACTA2', 'CD163', 'VWF')
png('integrated-marker-genes-atac.png', height=2.5*length(PLOT_GENES), width=2.5*3, units='in', res=300)
FeaturePlot(atac.integrated, c('MYH1', 'PAX7', 'MYH7'), pt.size = 0.2, split.by = 'individual')
dev.off()


# transfer cluster labels
transfer.anchors <- FindTransferAnchors(reference = rna.integrated, query = atac.integrated, features = VariableFeatures(object = rna.integrated), 
    reference.assay = "integrated", query.assay = "integrated", reduction = "cca")


celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna.integrated$seurat_clusters, 
    weight.reduction = atac.integrated[["pca"]], dims=1:20)
atac.integrated <- AddMetaData(atac.integrated, metadata = celltype.predictions)

atac.integrated@meta.data$predicted.id <- factor(atac.integrated$predicted.id, levels = levels(rna.integrated))

p1 <- DimPlot(atac.integrated, reduction = "umap", group.by='predicted.id')
p2 <- DimPlot(atac.integrated, reduction = "umap", group.by='individual')

pdf('atac-umap-by-cluster.pdf', height=5, width=5)
print(p1)
dev.off()

pdf('atac-umap-by-individual.pdf', height=5, width=5)
print(p2)
dev.off()

rna.integrated@meta.data$predicted.id <- rna.integrated$seurat_clusters

# co-embed
genes.use <- transfer.anchors@anchor.features
refdata <- GetAssayData(rna.integrated, assay = "integrated", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac.integrated[["pca"]], dims=1:20)

atac.integrated[["imputed"]] <- imputation
rna.integrated[["imputed"]] <- rna.integrated[['integrated']]
coembed <- merge(x = rna.integrated, y = atac.integrated)

DefaultAssay(coembed) <- 'imputed'

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:20) # 2 : 20 works well

coembed@meta.data$confident <- ifelse(coembed@meta.data$modality=='RNA', T, F)
coembed@meta.data$confident[coembed@meta.data$modality=='ATAC'] <- coembed@meta.data$prediction.score.max[coembed@meta.data$modality=='ATAC'] >= 0.7

pdf('atac-prediction-score-max.pdf', height=5, width=5)
hist(coembed@meta.data$prediction.score.max[coembed@meta.data$modality=='ATAC'])
dev.off()

p1 <- DimPlot(coembed, reduction = "umap", group.by = "modality")
p2 <- DimPlot(coembed, reduction = "umap", group.by = "predicted.id")
p3 <- DimPlot(coembed, reduction = "umap", group.by = "individual")
p4 <- DimPlot(coembed, reduction = "umap", group.by = "confident", split.by='modality')

pdf('joint-umap-by-modality.pdf', height=5, width=5)
print(p1)
dev.off()

pdf('joint-umap-by-cluster.pdf', height=5, width=5)
print(p2)
dev.off()

pdf('joint-umap-by-individual.pdf', height=5, width=5)
print(p3)
dev.off()

pdf('joint-umap-by-confidence.pdf', height=5, width=10)
print(p4)
dev.off()


DROP_NUCLEI <- rownames(coembed@meta.data[!coembed@meta.data$confident,])
KEEP_NUCLEI <- rownames(coembed@meta.data[coembed@meta.data$confident,])


clusters <- as.data.frame(coembed$predicted.id)
colnames(clusters) <- c('cluster')
clusters$nucleus <- rownames(clusters)
clusters$library <- gsub('(.*)-(.*)', '\\1', clusters$nucleus)
clusters$barcode <- gsub('(.*)-(.*)', '\\2', clusters$nucleus)
clusters$keep <- clusters$nucleus %in% KEEP_NUCLEI
clusters$genome <- gsub('.*-', '', clusters$library)
lost_nuclei <- as.data.frame(table(clusters[!clusters$keep,c('genome', 'cluster')]))
write.table(lost_nuclei, 'dropped-nuclei.txt', sep='\t', append=F, quote=F, row.names=F)
write.table(clusters[clusters$keep,c('library', 'barcode', 'cluster')], 'seurat-clusters.txt', sep='\t', append = F, quote = F, row.names = F, col.names = F)


# umap
umap <- as.data.frame(coembed@reductions$umap@cell.embeddings)
colnames(umap) <- c('dim1', 'dim2')
umap$nucleus <- rownames(umap)
umap$library <- gsub('(.*)-(.*)', '\\1', umap$nucleus)
umap$barcode <- gsub('(.*)-(.*)', '\\2', umap$nucleus)
umap$keep <- umap$nucleus %in% KEEP_NUCLEI
write.table(umap[umap$keep,c('library', 'barcode', 'dim1', 'dim2')], 'seurat-umap.txt', sep='\t', append = F, quote = F, row.names = F, col.names = F)

save(coembed, file="coembed.Rda")
save(atac.integrated, file="atac.integrated.Rda")
save(rna.integrated, file="rna.integrated.Rda")

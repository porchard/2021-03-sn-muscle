#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(rhdf5)
library(dplyr)
library(glue)

#NUCLEI_INDIVIDUAL_ASSIGNMENTS <- '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/nucleus-qc/nuclei-with-individual-assignments.txt'
#RNA_HDFS <- Sys.glob('/lab/work/porchard/sn-muscle-project/work/quality-nuclei-count-files/results/rna/*hg19*')

args <- commandArgs(T)
NUCLEI_INDIVIDUAL_ASSIGNMENTS <- args[1]
RNA_HDFS_GLOB <- args[2]
RNA_HDFS <- Sys.glob(RNA_HDFS_GLOB)



# In[2]:


individual_assignments <- read.table(NUCLEI_INDIVIDUAL_ASSIGNMENTS, sep='\t', head=F, col.names=c('library', 'barcode', 'individual'))
individual_assignments$nucleus <- with(individual_assignments, paste(library, barcode, sep='-'))
individual_assignments <- individual_assignments[,c('individual', 'nucleus')]
head(individual_assignments)


# In[4]:


load_hdf <- function(HDF) {
    tmp <- h5ls(HDF)
    GROUP <- paste("/", tmp$name[tmp$group=='/'], sep='')
    df <- h5read(HDF, GROUP, native = F, compoundAsDataFrame=T)
    counts <- t(df$block0_values)
    rownames(counts) <- as.character(df$axis1)
    colnames(counts) <- as.character(df$axis0)
    return(counts)
}


# In[5]:


rna_objects <- list()
for(f in RNA_HDFS) {
    print(f)
    rna_counts <- load_hdf(f)
    lib <- gsub('.hdf5', '', basename(f))
    lib <- gsub('-', '', lib)
    lib <- gsub('_', '', lib)
    metadata <- data.frame(library=rep(lib, nrow(rna_counts)))
    rownames(metadata) <- rownames(rna_counts)
    metadata$nucleus <- rownames(metadata)
    metadata <- left_join(metadata, individual_assignments)
    rna <- CreateSeuratObject(counts = t(rna_counts), min.cells=5, assay = "RNA", project = lib, metadata=metadata)
    rna$tech <- "rna"
    rna$library <- lib
    rna$individual <- metadata$individual
    rna_objects[[length(rna_objects)+1]] <- rna
}


# In[6]:

if (length(rna_objects) > 1) {
  rna_additional <- c()
  for(i in 2:length(rna_objects)) {
      rna_additional <- c(rna_additional, rna_objects[[i]])
  }
  rna <- merge(rna_objects[[1]], y = rna_additional, project = "RNA")
  rna <- SplitObject(rna, split.by = "individual")
} else {
  rna <- rna_objects
  names(rna) <- 'rat1'
}


for (i in 1:length(rna)) {
    print(i)
    rna[[i]] <- NormalizeData(rna[[i]], verbose = FALSE)
    rna[[i]] <- FindVariableFeatures(rna[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}


# In[10]:


# try each individual separately
for (i in 1:length(rna)) {
    print(i)
    rna[[i]] <- ScaleData(rna[[i]], verbose = FALSE)
    rna[[i]] <- RunPCA(rna[[i]], npcs = 20, verbose = FALSE)
    rna[[i]] <- RunUMAP(rna[[i]], reduction = "pca", dims = 1:20)
}


# In[11]:


#PLOT_GENES <- c('MYH1', "MYH2", "MYH7",'CHRNA1', 'VWF', 'PAX7')
PLOT_GENES <- c('Myh1', 'Tnnt3', "Tnnt1", "Myh1",'Myh4', 'Chrna1', 'Lrrk2')
for(i in 1:length(rna)) {
    individual <- names(rna)[i]
    png(glue('{individual}-marker-genes.png'), height=8, width=5, units='in', res=300)
    FeaturePlot(rna[[i]], PLOT_GENES, pt.size = 0.2, order=F)
    dev.off()
}

# TODO: change if we end up with data from more than one rat
rna.integrated <- rna[[1]]

p2 <- DimPlot(rna.integrated, reduction = "umap", group.by = "individual")
png('integrated-individuals.png', height=12, width=6, units='in', res=300)
p2
dev.off()


# In[16]:


#PLOT_GENES <- c('MYH1', 'MYH4', "CHRNA1", "MYH7",'LRRK2', 'TNNT1', 'TNNT3', 'PAX7', 'PDGFRB', 'MYH11', 'ACTA2', 'CD163', 'VWF')
PLOT_GENES <- c('Myh1', 'Tnnt3', "Tnnt1", "Myh1",'Myh4', 'Chrna1', 'Lrrk2')
png('integrated-marker-genes.png', height=2.5*length(PLOT_GENES), width=6, units='in', res=300)
FeaturePlot(rna.integrated, PLOT_GENES, pt.size = 0.2, split.by = 'individual')
dev.off()


# In[17]:


rna.integrated <- FindNeighbors(rna.integrated, dims = 1:20, k.param = 10)
rna.integrated <- FindClusters(rna.integrated, resolution = 0.1, n.start = 100)
png('integrated-clusters.png', height=6, width=6, units='in', res=300)
DimPlot(rna.integrated, reduction = "umap")
dev.off()


# In[18]:


clusters = as.data.frame(rna.integrated@active.ident)
colnames(clusters) <- c('cluster')
clusters$nucleus <- rownames(clusters)
clusters$barcode <- gsub('.*-(.*)', '\\1', clusters$nucleus)
write.table(clusters[,c('nucleus', 'cluster')], file = 'clusters.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)


umap <- as.data.frame(rna.integrated@reductions$umap@cell.embeddings)
colnames(umap) <- c('dim1', 'dim2')
umap$nucleus <- rownames(umap)
write.table(umap[,c('nucleus', 'dim1', 'dim2')], file='umap.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)

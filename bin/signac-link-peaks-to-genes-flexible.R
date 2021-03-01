#!/usr/bin/env Rscript
# coding: utf-8

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(rhdf5)

set.seed(1234)

ATAC_BARCODE_LIST = '/home/porchard/github/snATACseq-NextFlow/737K-arc-v1.txt'
RNA_BARCODE_LIST = '/home/porchard/github/snRNAseq-NextFlow/737K-arc-v1.txt'
LIBRARY <- '1846_RNA-hg19' # NOTE: if change this, must change the 'annotation' information further down as it may no longer be hg19

args <- commandArgs(T)

CLUSTERS <- args[1]
TEST_CLUSTER <- args[2] # a cluster ID, or 'all'
TEST_GENES <- args[3] # file listing genes to test
ATAC_HDF5 <- args[4]
RNA_HDF5 <- args[5]
METHOD <- args[6] # spearman, pearson
DISTANCE <- as.numeric(args[7])
MIN_CELLS <- as.numeric(args[8])
OUT <- args[9]

TEST_GENES <- read.table(TEST_GENES, header=F)$V1

# for testing
#CLUSTERS <- '/lab/work/porchard/Nova-315/work/process-by-cluster-grant-2021-01/clusters.txt'
#TEST_CLUSTER <- 3 # a cluster ID, or 'all'
#TEST_GENES <- c('MYH1', 'ARL15', 'FST', 'FBN1')
#ATAC_HDF5 <- '/lab/work/porchard/Nova-315/work/hd5-count-matrices/atac_hg19.hdf5'
#RNA_HDF5 <- '/lab/work/porchard/Nova-315/work/hd5-count-matrices/rna_hg19.hdf5'


# load in the RNA

load_hdf <- function(HDF) {
  tmp <- h5ls(HDF)
  GROUP <- paste("/", tmp$name[tmp$group=='/'], sep='')
  df <- h5read(HDF, GROUP, native = F, compoundAsDataFrame=T)
  counts <- t(df$block0_values)
  rownames(counts) <- toupper(as.character(df$axis1))
  colnames(counts) <- as.character(df$axis0)
  return(t(counts))
}


clusters <- read.table(CLUSTERS, header=F, col.names=c('library', 'barcode', 'cluster'))
clusters$cluster <- as.character(clusters$cluster)
clusters <- clusters[clusters$library==LIBRARY,]
KEEP_BARCODES <- c()
if (TEST_CLUSTER == 'all') {
  KEEP_BARCODES <- clusters$barcode
} else {
  KEEP_BARCODES <- clusters$barcode[clusters$cluster==TEST_CLUSTER]  
}


atac_barcodes <- read.table(ATAC_BARCODE_LIST, header=F)$V1
rna_barcodes <- read.table(RNA_BARCODE_LIST, header=F)$V1
barcodes <- data.frame(atac=atac_barcodes, rna=rna_barcodes)
rownames(barcodes) <- barcodes$atac

ATAC <- load_hdf(ATAC_HDF5)
# switch ATAC barcodes --> RNA barcodes
colnames(ATAC) <- barcodes[colnames(ATAC),'rna']

RNA <- load_hdf(RNA_HDF5)
rownames(RNA) <- gsub('.*\t', '', rownames(RNA))
colnames(RNA) <- gsub('.*-', '', colnames(RNA))


# In[4]:


# TODO: use GTF rather than ensembl??
#GTF <- '/lab/work/porchard/sn-muscle-project/data/star/hg19/hg19.gtf'
#annotation <- makeGRangesFromGTF(
#  GTF,
#  level = c("genes", "transcripts"),
#  ignoreTxVersion = TRUE,
#  .checkAgainstTxDb = FALSE
#)


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg19"

# create a Seurat object containing the RNA data
pbmc <- CreateSeuratObject(
  counts = RNA,
  assay = "RNA"
)


# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = ATAC,
  sep = c(":", ":"),
  #fragments = fragpath,
  annotation = annotation
)


DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
DefaultAssay(pbmc) <- "ATAC"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg19)

cells_of_interest <- subset(pbmc, cells=KEEP_BARCODES)

# link peaks to genes
cells_of_interest <- LinkPeaks(
  object = cells_of_interest,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = TEST_GENES,
  pvalue_cutoff=1,
  score_cutoff=0,
  min.cells=MIN_CELLS,
  method=METHOD,
  distance=DISTANCE
)

links <- as.data.frame(cells_of_interest@assays$ATAC@links)
links$peak <- gsub('-', ':', links$peak)
write.table(links, file=OUT, append = F, quote = F, sep = '\t', row.names = F, col.names = T)


#!/usr/bin/env Rscript
# dyn.load('/home/porchard/anaconda3/lib/libicui18n.so.58')
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(tidyr)
library(Matrix)

args <- commandArgs(T)
COUNTS_GLOB <- args[1]
MOTIF_GLOB <- args[2]
CLUSTERS <- args[3]
MOTIF_TO_TF <- args[4]

# for testing
#COUNTS_GLOB <- '/lab/work/porchard/sn-muscle-project/work/chromvar/results/fragment-counts/*'
#MOTIF_GLOB <- '/lab/work/porchard/sn-muscle-project/work/chromvar/results/motifs-in-peaks/*'
#MOTIF_TO_TF <-'/lab/work/porchard/sn-muscle-project/data/cisbp-motifs/motif-id-to-tf-name.txt'
#CLUSTERS <- '/lab/work/porchard/sn-muscle-project/work/downstream-new-features/results/liger/round-1/second-louvain-clusters.txt'


motif2tf <- read.table(MOTIF_TO_TF, head=F, as.is=T, sep='\t', col.names=c('motif', 'tf'))

nucleus2cluster <- read.table(CLUSTERS, head=F, as.is=T, sep='\t', col.names=c('library', 'barcode', 'cluster'))
nucleus2cluster$nucleus <- with(nucleus2cluster, paste(library, barcode, sep='-'))
nucleus2cluster <- nucleus2cluster[,c('nucleus', 'cluster')]

COUNT_FILES <- Sys.glob(COUNTS_GLOB)
MOTIF_FILES <- Sys.glob(MOTIF_GLOB)

# Make sparse matrix of fragment counts
parse_count_file_name <- function(f) {
  return(gsub('.counts.bed', '', basename(f)))
}

load_count_file <- function(f) {
  tmp <- read.table(f, head=F, as.is=T, sep='\t', col.names=c('chrom', 'start', 'end', 'barcode', 'count'), colClasses = c('character', 'character', 'character', 'character', 'numeric'))
  tmp$peak <- with(tmp, paste(chrom, start, end, sep=':'))
  tmp$nucleus <- paste(parse_count_file_name(f), tmp$barcode, sep='-')
  return(tmp[,c('peak', 'nucleus', 'count')])
}

counts <- bind_rows(lapply(COUNT_FILES, load_count_file))
counts <- counts[counts$nucleus %in% nucleus2cluster$nucleus,]
peak_indices <- counts %>% dplyr::select(peak) %>% unique() %>% mutate(peak_index=1:n())
nucleus_indices <- counts %>% dplyr::select(nucleus) %>% unique() %>% mutate(nucleus_index=1:n())
counts <- left_join(counts, peak_indices) %>% left_join(nucleus_indices)
count_matrix <- sparseMatrix(i=counts$peak_index, j=counts$nucleus_index, x=counts$count)
colnames(count_matrix) <- nucleus_indices$nucleus

nucleus_indices <- left_join(nucleus_indices, nucleus2cluster)

# get peaks, from the counts matrix
peaks <- dplyr::select(peak_indices, peak) %>% tidyr::separate(col=peak, into=c('chrom', 'start', 'end'), convert=T)
PEAK_TMP_FILE <- 'peaks.tmp.bed'
write.table(peaks, PEAK_TMP_FILE, quote=F, append = F, sep = '\t', row.names = F, col.names = F)
peaks <- getPeaks(PEAK_TMP_FILE, sort_peaks=F)
file.remove(PEAK_TMP_FILE)

fragment_counts <- SummarizedExperiment(assays = list(counts = count_matrix), rowRanges = peaks)
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg19)
# fragment_counts <- filterPeaks(fragment_counts, non_overlapping = TRUE)

# Get motif locations
parse_motif_file_name <- function(f) {
  return(gsub('.in-peaks.bed', '', basename(f)))
}

load_motif_file <- function(f) {
  tmp <- read.table(f, head=F, as.is=T, sep='\t', col.names=c('chrom', 'start', 'end'), colClasses = c('character', 'numeric', 'numeric'))
  tmp$motif <- parse_motif_file_name(f)
  return(tmp)
}

motifs <- bind_rows(lapply(MOTIF_FILES, load_motif_file))
MOTIF_TMP_FILE <- 'motifs.tmp.bed'
write.table(motifs, MOTIF_TMP_FILE, quote=F, append = F, sep = '\t', row.names = F, col.names = F)
motifs <- getAnnotations(MOTIF_TMP_FILE, rowRanges = rowRanges(fragment_counts), column = 4)
file.remove(MOTIF_TMP_FILE)

# compute deviations
dev <- computeDeviations(object = fragment_counts, annotations = motifs)
saveRDS(dev, file = 'chromvar.RDS')
deviation_zscores <- as.data.frame(deviationScores(dev))
deviation_zscores$motif <- rownames(deviation_zscores)
deviation_zscores <- deviation_zscores %>% tidyr::gather(key = 'nucleus', value='deviation_z', -motif)
write.table(deviation_zscores, file='deviation-zscores.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

variability <- computeVariability(dev) %>% left_join(motif2tf %>% dplyr::rename(name=motif))
variability <- variability[order(variability$variability),]
write.table(variability, 'variability.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
png(file = 'variability.png', width = 10, height = 5, units = 'in', res=300)
plotVariability(variability, use_plotly = FALSE)
dev.off()

diff_acc <- differentialDeviations(dev, nucleus_indices$cluster)
diff_acc$motif <- rownames(diff_acc)
diff_acc <- left_join(diff_acc, motif2tf)
diff_acc <- diff_acc[order(diff_acc$p_value_adjusted),]
write.table(diff_acc, 'diff-deviation.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

diff_var <- differentialVariability(dev, nucleus_indices$cluster)
diff_var$motif <- rownames(diff_var)
diff_var <- left_join(diff_var, motif2tf)
diff_var <- diff_var[order(diff_var$p_value_adjusted),]
write.table(diff_var, 'diff-variation.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

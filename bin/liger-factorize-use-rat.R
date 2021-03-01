#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--mats"), action = "store", type = "character", default = "", help = "[Required] (comma-separated) paths to the *.txt files for input (columns are nuclei, rows are genes"),
  make_option(c("--factorization_k"), action = "store", type = "numeric", default = 20, help = "[Required] k parameter for factorization"),
  make_option(c("--factorization_lambda"), action = "store", type = "numeric", default = 5, help = "[Required] lambda parameter for factorization"),
  make_option(c("--select_genes"), action = "store", type = "character", default = '', help = "[Required] Regex for matrices used to select genes."),
  make_option(c("--out"), action = "store", type = "character", default = '', help = "[Required] Output .Rda file.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# for testing
# opts <- list(
#   'mats' = paste(list.files('/lab/work/porchard/sn-muscle-project/work/downstream-final/results/liger/round-1/input-merged/', full.names=T), collapse=','),
#   'factorization_k' = 25,
#   'factorization_lambda' = 5,
#   'prefix' = 'test'
# )

PREFIX <- opts$prefix

count_files <- sort(unlist(strsplit(opts$mats, ',')))

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(liger)
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

remove_genes_with_no_signal_in_n_cells <- function(m, n) {
	  nonzero <- apply(m, 1, function(x){sum(x>0)})
  return(m[nonzero>=n,])
}

liger_in <- lapply(count_files, load_hdf)
names(liger_in) <- sapply(count_files, function(f){gsub('.hdf5', '', basename(f))})
liger_in <- lapply(liger_in, function(m){remove_genes_with_no_signal_in_n_cells(m, 10)})

print('Selecting genes with:')
SELECT_GENES <- opts$select_genes
print(SELECT_GENES)

int.muscle <- createLiger(liger_in)
int.muscle <- normalize(int.muscle)
USE <- unique(sort(c(grep(SELECT_GENES, names(liger_in)), grep('rat1_RNA', names(liger_in))))) # use RNA to select genes and ignore ATAC
print(names(liger_in)[USE])
int.muscle <- selectGenes(int.muscle, datasets.use = USE) # use RNA to select genes and ignore ATAC
int.muscle <- scaleNotCenter(int.muscle)

# factorization
print('Factorizing')
int.muscle <- optimizeALS(int.muscle, k = opts$factorization_k, lambda = opts$factorization_lambda, nrep=4)
print('Finished factorizing')
save(int.muscle, file=opts$out)
print('Finished saving object')

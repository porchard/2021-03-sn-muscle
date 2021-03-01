library(dplyr)

options(scipen=999)

args <- commandArgs(T)
OUT <- args[1]
CHROM_SIZES <- args[2]
HIC_FILES <- args[3:length(args)]


# for testing
#CHROM_SIZES <- '/lab/work/porchard/data/chrom_sizes/hg19.chrom_sizes'
#HIC_FILES <- list.files('/lab/work/porchard/sn-muscle-project/data/hiC/FitHiC_primary_cohort', pattern = 'FitHiC_output.PO_', full.names = T)

hic_file_to_chromosome <- function(f) {
  return(gsub('FitHiC_output.*_(chr.*).sparse.matrix.gz', '\\1', basename(f)))
}

load_hic_file <- function(f) {
  tmp <- read.table(gzfile(f), head = T, as.is = T, sep = '\t')
  tmp$chrom <- hic_file_to_chromosome(f)
  return(tmp)
}


hiC <- bind_rows(lapply(HIC_FILES, load_hic_file))
chrom_sizes <- read.table(CHROM_SIZES, head = F, as.is = T, sep = '\t', col.names = c('chrom', 'size'))
chrom_sizes <- chrom_sizes[chrom_sizes$chrom %in% hiC$chrom,]

# make the bins...

bins <- bind_rows(apply(chrom_sizes,
              1,
              function(x){
                chrom <- x[1]
                size <- as.numeric(x[2])
                bin_starts <- seq(0, size, by=40e3)
                bin_ends <- c(seq(40e3, size, by=40e3), size)
                tmp <- data.frame(
                  chrom=chrom,
                  start=bin_starts,
                  end=bin_ends,
                  bin=1:length(bin_starts)
                )
                return(tmp)
              }))
bins$bin_label <- with(bins, paste(chrom, start, end, sep = ':'))

hiC <- left_join(
  hiC,
  bins %>%
    dplyr::select(chrom, bin, bin_label) %>%
    dplyr::rename(
      RowID=bin,
      bin_1=bin_label
    )
)

hiC <- left_join(
  hiC,
  bins %>%
    dplyr::select(chrom, bin, bin_label) %>%
    dplyr::rename(
      ColumnID=bin,
      bin_2=bin_label
    )
)

hiC <- dplyr::select(hiC, bin_1, bin_2, ObservedCount, ExpectedCount, OERatio, PValue, QValue)

write.table(hiC, file = OUT, append = F, quote = F, sep = '\t', row.names = F, col.names = T)

#!/usr/bin/env Rscript

# set up dropletUtils (in seurat-v4 environment)
library(DropletUtils)
library(glue)

# SPARSE_DIR <- '/lab/work/porchard/sn-muscle-project/work/test-droplet-utils'
# SAMPLE <- '133155-hg19'

args <- commandArgs(T)
SPARSE_DIR <- args[1]
SAMPLE <- args[2]
FDR_TRESHOLD <- 0.0001

samp <- read10xCounts(samples=c(SPARSE_DIR), col.names=T, sample.names=SAMPLE, version='3', type='sparse')


# make knee plot. using code from vignette.
br.out <- barcodeRanks(counts(samp))
FILENAME <- glue('{SAMPLE}.knee.png')
png(FILENAME, width = 5, height = 5, units='in', res = 300)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))
dev.off()


# test for empty droplets
set.seed(2020)
e.out <- emptyDrops(counts(samp), niters=1e5)
e.out$is.cell <- e.out$FDR <= FDR_TRESHOLD
stopifnot(sum(e.out$Limited & !is.na(e.out$is.cell) & !e.out$is.cell, na.rm=T)==0)
e.out$barcode <- rownames(e.out)
FILENAME <- glue('{SAMPLE}.test-vs-ambient.txt')
write.table(e.out, file = FILENAME, append = F, quote = F, sep = '\t', row.names = F, col.names = T)

is.cell <- e.out$FDR <= FDR_TRESHOLD
keep <- e.out[!is.na(e.out$FDR) & e.out$FDR<=FDR_TRESHOLD,]
keep$barcode <- rownames(keep)
FILENAME <- glue('{SAMPLE}.keep.txt')
write.table(keep[,c('barcode'),drop=F], file = FILENAME, append = F, quote = F, sep = '\t', row.names = F, col.names = F)


# vignette: If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
RERUN <- sum(e.out$Limited & !is.na(is.cell) & !is.cell)
stopifnot(RERUN==0)


XMAX <- 1000
YLIM <- 1.2 * max(-1*e.out$LogProb[e.out$Total<=XMAX], na.rm=T)
FILENAME <- glue('{SAMPLE}.logprob.png')
png(FILENAME, width = 5, height = 5, units='in', res = 300)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability", xlim=c(0, XMAX), ylim=c(0, YLIM))
dev.off()

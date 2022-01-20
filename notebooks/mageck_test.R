library("MAGeCKFlute")
# also needs sva
countMat <- read.table('mageck_dirty_counts.txt', header=TRUE, sep='\t') # can't have NA
batMat <- read.table('mageck_dirty_14_2_batch.txt', header=TRUE)
x <- BatchRemove(countMat, batMat, outdir=".", pca=TRUE)
x$data[x$data<0] <- 0
write.table(x$data, "14_2.dirty.batch.corrected.txt", row.names=FALSE, quote=FALSE, sep="\t")
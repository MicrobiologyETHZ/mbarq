library(MAGeCKFlute)
args <- commandArgs(trailingOnly=TRUE)
# count file
count_file <-args[1]
# batch file
batch_file <- args[2]
library <- args[3]
outDir <- args[4]

batchMat <- data.frame(read.table(batch_file, header=TRUE))
countMat <- data.frame(read.table(count_file, header=TRUE))
x <- MAGeCKFlute::BatchRemove(mat=countMat, batchMat=batchMat, pca=TRUE, outdir=outDir)
x$data[x$data < 0] <- 0
write.table(x$data, file.path(outDir, paste0(library, "_count.batchcorrected.txt")), quote=FALSE,
row.names=FALSE, sep='\t')
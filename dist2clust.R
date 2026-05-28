#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3) {
	cat("Usage: Rscript dist2clust.R distance.matrix hclust.method cluster.number\n")
	quit()
}
# read distance matrix
dt <- read.table(args[1], header=TRUE, row.names=1)
dt <- as.dist(dt, diag=TRUE, upper=TRUE)
# hclust
num = as.numeric(args[3])
cl <- hclust(dt, method=args[2])
gp <- cutree(cl, k=num)
for(i in 1:num) {
	cat(c(labels(gp[gp == i])), sep="\t", fill=TRUE)
}

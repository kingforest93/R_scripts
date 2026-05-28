#!/usr/bin/env Rscript
#load package
library(ggplot2)
library(RColorBrewer)

#read data
args <- commandArgs(TRUE)
if (length(args) != 1) {
	cat("Rscript stat2hist.R ccs.read.fq.gz.stat", fill=TRUE)
	quit()
}
dat <- read.table(args[1], sep="\t", header=FALSE, skip=1, row.names=1, col.names=c("id", "len", "rq", "np"))

#draw histogram of read length
hist.p <- ggplot(dat, aes(x=len)) +
			geom_histogram(bins=50, fill="blue", color="grey", size=0.1)+
			labs(x="Read length (bp)", y="Number of reads") +
			theme_bw()
ggsave(paste(args[1], ".len.pdf", sep=""), hist.p, device="pdf", width=12, height=8, units="cm")

#draw hitogram of read quality
hist.p <- ggplot(dat, aes(x=rq)) +
            geom_histogram(bins=50, fill="blue", color="grey", size=0.1)+
            labs(x="Read accuracy", y="Number of reads") +
            theme_bw()
ggsave(paste(args[1], ".rq.pdf", sep=""), hist.p, device="pdf", width=12, height=8, units="cm")

#draw histogram of pass number
hist.p <- ggplot(dat, aes(x=np)) +
            geom_histogram(bins=50, fill="blue", color="grey", size=0.1)+
            labs(x="Number of passes", y="Number of reads") +
            theme_bw()
ggsave(paste(args[1], ".np.pdf", sep=""), hist.p, device="pdf", width=12, height=8, units="cm")


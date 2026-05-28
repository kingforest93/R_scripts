#!/usr/bin/env Rscript
#load package
library(ggplot2)
library(RColorBrewer)

#read data
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1){
	cat("Rscript depth_line.R assembly.duplicate_busco.tsv", fill=TRUE)
	quit()
}
dat <- read.table(args[1], sep="\t", header=TRUE)
dat$gene_contig <- as.factor(dat$gene_contig)

#draw line plot of read depth and Duplicate BUSCOs
line.p <- ggplot(dat, aes(x=position, y=depth)) +
			geom_line(color="blue", size=1.5)+
			labs(x="Contig position", y="Read depth") +
			facet_wrap(~ gene_contig, ncol=2, scales="free_x") +
			theme_bw() + theme(legend.position="none")
args[1] <- sub("tsv", "", args[1])
ggsave(paste(args[1], "read_depth.pdf", sep=""), line.p, device="pdf", width=60, height=1500, units="cm", limitsize=FALSE)


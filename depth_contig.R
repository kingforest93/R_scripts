#!/usr/bin/env Rscript
#load package
library(ggplot2)
library(RColorBrewer)

#read data
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
	cat("Rscript depth_contig.R contig.base.cov xaxis_tick_interval(bp)\n")
	quit()
}
dat <- read.table(args[1], sep="\t", header=FALSE)
colnames(dat) <- c("contig", "position", "depth")
dat$contig <- as.factor(dat$contig)
win <- as.numeric(args[2])
max_pos <- (as.numeric(dat$position[nrow(dat)]) %/% win + 1) * win

#draw line plot of read depth and Duplicate BUSCOs
line.p <- ggplot(dat, aes(x=position, y=depth)) +
			geom_line(color="blue", size=1.5)+
			labs(x="Contig position", y="Read depth") +
			scale_x_continuous(limits=c(0, max_pos), breaks=seq(0, max_pos, win), expand=c(0, 0)) +
			facet_wrap(~ contig, ncol=1, scales="free") +
			theme_bw() + theme(text=element_text(size=20), legend.position="none")
args[1] <- sub("tsv", "", args[1])
ggsave(paste(args[1], ".pdf", sep=""), line.p, device="pdf", width=1500, height=30, units="cm", limitsize=FALSE)

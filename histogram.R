#!/usr/bin/env Rscript

#load package
library(ggplot2)
library(RColorBrewer)

#read data
args <- commandArgs(TRUE)
if (length(args) != 1) {
    cat("Rscript histogram.R PB.base.cov.50kb.dp\n")
	quit()
}
dat <- read.table(args[1], sep="\t", header=TRUE)
colnames(dat)=c("contig", "start", "stop", "depth")
dat$contig <- cut(dat$depth, breaks=c(-Inf, 3, seq(21, 91, 14), Inf), labels=c("repeat", "monoplotig", "diplotig", "triplotig", "tetraplotig", "pentaplotig", "hexaplotig", "repeat"))

#draw histogram
hist.p <- ggplot(dat, aes(x=depth)) +
			geom_histogram(aes(fill=contig), binwidth=1)+
			#geom_density(fill="blue", color="blue", alpha=0.8)+
			#xlim(0, 100)+
			#labs(x="Average depth window", y="Number of windows") +
			scale_x_continuous("Average depth of window", breaks=seq(0, 100, 14), limits=c(0, 100)) +
			#scale_y_continuous("Number of windows", breaks=seq(0, 100000, 20000), limits=c(0, 100000)) +
			labs(y="Number of windows") +
			theme_bw()
ggsave(paste0(args[1], ".hist.pdf"), hist.p, device="pdf", width=20, height=10, units="cm")


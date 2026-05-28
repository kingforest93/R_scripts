#!/public/software/R-3.3.1/bin/Rscript
#load package
library(ggplot2)
library(RColorBrewer)

#read data
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3){
	cat("Rscript line_plot.R gff.stat column axis_tick", fill=TRUE)
	quit()
}
dat <- read.table(args[1], sep="\t", header=FALSE, col.names=c("id", "tlen", "start", "end", "number", "length"))
col <- as.integer(args[2])
tick <- as.integer(args[3])

#draw line plot of feature length in each window along the contig
line.p <- ggplot(dat, aes(x=end, y=dat[,col])) +
			geom_line(color="blue", size=1.5)+
			labs(x="Contig position", y=paste("Feature", as.character(colnames(dat)[col]), sep=" ")) +
			scale_x_continuous(limits=c(0, max(dat$end) + 1), breaks=seq(tick, max(dat$end), tick), expand=c(0, 0)) +
			theme_bw() + theme(legend.position="none")
ggsave(paste(args[1], ".pdf", sep=""), line.p, device="pdf", width=200, height=40, units="cm", limitsize=FALSE)


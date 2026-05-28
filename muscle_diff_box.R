#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
    cat("Usage: Rscript muscle_diff_box.R m1.muscle.diff m2.muscle.diff m3.muscle.diff ...\n")
    quit()
}

# read multiple muscle.diff
rk <- c()
iden <- c()
for (i in args) {
	if (file.info(i)$size == 0) {
		next
	}
	tem <- read.table(i, sep = "\t", header = FALSE)
	tem <- c(tem$V10)
	iden <- c(iden, tem[order(tem, decreasing = TRUE)])
	rk <- c(rk, c(1:length(tem)))
}
dt <- data.frame(rk = rk, iden = iden)

# draw violin and boxplot
gp <- ggplot(data = dt, aes(x = factor(rk), y = iden)) +
	geom_violin(trim = TRUE, fill = "grey") +
	geom_boxplot(width = 0.1, fill = "white") +
	labs(title = "Alignment identity of sequences", x = "Relative rank (decreasing)", y = "Identity (%)") +
	theme_bw()
ggsave(filename = "muscle.diff.violin_boxplot.pdf", plot = gp, device = "pdf", width = 14, height = 10)


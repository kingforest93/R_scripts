#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
    cat("Usage: Rscript muscle.diff_Ks_box.R OG1.muscle.diff.axt.KaKs OG2.muscle.diff.axt.KaKs OG3.muscle.diff.axt.KaKs ...\n")
    quit()
}

# read multiple muscle.diff
max_ks = args[1]
group <- c()
ks <- c()
num <- 0
for (i in args[2:length(args)]) {
	tem <- read.table(i, sep = "\t", header = FALSE, na.strings = "NA")
	tem <- tem[, c(1,4)]
	colnames(tem) <- c('pair', 'ks')
	tem[is.na(tem)] = 0
	if (max(tem$ks) > max_ks) {
		next
	}
	#tem <- subset(tem, ks <= max_ks)
	t <- unlist(strsplit(tem$pair, "&"))
	for (i in 1:nrow(tem)) {
		tem$pair[i] = paste(substr(t[2 * i - 1], 12, 13), substr(t[2 * i], 12, 13), sep = "&")
	}
	group <- c(group, tem$pair)
	ks <- c(ks, tem$ks)
	num <- num + 1
}
dt <- data.frame(group = group, ks = ks)

# order chromosome pairs
tem <- aggregate(dt$ks, by = list(type = dt$group), median)
tem <- tem[order(tem$x, decreasing = FALSE), ]
dt$group <- factor(dt$group, levels=tem$type)

# draw violin and boxplot
gp <- ggplot(data = dt, aes(x = group, y = ks)) +
	geom_violin(trim = TRUE, fill = "grey") +
	geom_boxplot(width = 0.1, fill = "white") +
	labs(title = paste0("Max Ks: ", max_ks, " Num OGs: ", num), x = "Chromosome pair", y = "Synonymous substitution rate (Ks)") +
	theme_bw()
ggsave(filename = "combined.muscle.diff.Ks.violin_boxplot.pdf", plot = gp, device = "pdf", width = 12, height = 8)


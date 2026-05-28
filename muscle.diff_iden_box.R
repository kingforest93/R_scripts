#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
    cat("Usage: Rscript muscle.diff_iden_box.R OG1.muscle.diff OG2.muscle.diff OG3.muscle.diff ...\n")
    quit()
}

# read multiple muscle.diff
min_iden = 90
group <- c()
identity <- c()
num <- 0
for (i in args[1:length(args)]) {
	tem <- read.table(i, sep = "\t", header = FALSE)
	tem <- tem[, c(1,10)]
	colnames(tem) <- c('pair', 'iden')
	if (min(tem$iden) <= min_iden) {
		next
	}
	#tem <- subset(tem, iden >= min_iden)
	t <- unlist(strsplit(tem$pair, "_vs_"))
	for (i in 1:nrow(tem)) {
		tem$pair[i] = paste(substr(t[2 * i - 1], 12, 13), substr(t[2 * i], 12, 13), sep = "&")
	}
	group <- c(group, tem$pair)
	identity <- c(identity, tem$iden)
	num <- num + 1
}
dt <- data.frame(group = gsub("&", "vs", group), identity = identity)

# order chromosome pairs
#cat(num, "OGs\n")
#aggregate(dt$identity, by = list(type = dt$group), median)
#quit()

#tem <- aggregate(dt$identity, by = list(type = dt$group), median)
#tem <- tem[order(tem$x, decreasing = TRUE), ]
#dt$group <- factor(dt$group, levels=tem$type)

# draw violin and boxplot
gp <- ggplot(data = dt, aes(x = group, y = identity)) +
	geom_violin(trim = TRUE, fill = "grey") +
	geom_boxplot(width = 0.1, fill = "white") +
	labs(title = paste0("Min identity: ", min_iden, " Num OGs:", num), x = "Chromosome pair", y = "Identity (%)") +
	ylim(90, 100) +
	scale_x_discrete(limits = c("H1vsH2", "H3vsH4", "H5vsH6", "H1vsH3", "H2vsH3", "H1vsH4", "H2vsH4", "H1vsH5", "H2vsH5", "H3vsH5", "H4vsH5", "H1vsH6", "H2vsH6", "H3vsH6", "H4vsH6")) +
	theme_bw()
ggsave(filename = "combined.muscle.diff.identity.boxplot.pdf", plot = gp, device = "pdf", width = 12, height = 8)


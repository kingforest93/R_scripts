#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
    cat("Usage: Rscript coords_iden_box.R align1.1coords align2.1coords align3.1coords ...\n")
    quit()
}
min_aln <- 100
min_iden <- 80

# read mummer.delta.1v1.coords.tab and filter
pair <- c()
identity <- c()
for (i in args) {
	tem1 <- read.table(i, sep = "\t", header = FALSE)
	tem1 <- tem1[, 5:7]
	colnames(tem1) <- c('taln', 'qaln', 'iden')
	tem1 <- subset(tem1, taln >= min_aln & qaln >= min_aln & iden >= min_iden)
	tem2 <- sub(".delta.1coords", "", i)
	tem2 <- sub("Htub.Chr[0-9][0-9].", "", tem2)
	tem2 <- sub("_to_", "vs", tem2)
	pair <- c(pair, rep(tem2, nrow(tem1)))
	identity <- c(identity, tem1$iden)
}
dt <- data.frame(pair = pair, identity = identity)

# draw violin and boxplot
gp <- ggplot(data = dt, aes(x = pair, y = identity)) +
	geom_violin(trim = TRUE, fill = "grey") +
	geom_boxplot(width = 0.1, fill = "white") +
	labs(title = "Alignment identity of sequences", x = "Sequence pairs", y = "Identity (%)") +
	#scale_x_discrete(limits = c("H2vsH1", "H4vsH3", "H6vsH5", "H4vsH2", "H4vsH1", "H3vsH2", "H3vsH1", "H6vsH4", "H6vsH3", "H6vsH2", "H6vsH1", "H5vsH4", "H5vsH3", "H5vsH2", "H5vsH1")) +
	scale_x_discrete(limits = c("H1vsHa", "H2vsHa", "H3vsHa", "H4vsHa", "H5vsHa", "H6vsHa")) +
	theme_bw()
f <- args[1]
f <- sub(".H[1-6]_to_H[1-6a].delta.1coords", "", f)
ggsave(filename = paste0(f, ".6H_Ha.1coords.boxplot.pdf"), plot = gp, device = "pdf", width = 14, height = 8)


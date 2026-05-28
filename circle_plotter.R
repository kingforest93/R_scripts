#!/usr/bin/env Rscript

# load package
suppressPackageStartupMessages(library(circlize))
library(scales)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[2]) ) {
	cat("Usage: Rscript circle_plotter.R G.cor_genome_chr.fa.len G.cor_genome_chr.pep.collinearity.filtered.list\n")
	quit()
}

# read chromosome.len and collinearity.list
dflink <- read.table(args[2], sep = "\t", header = FALSE)
dflink <- dflink[as.integer(row.names(dflink)) %% 2 == 1, c(2,4,5,7,9,10)]
dflink_from <- data.frame(chr = dflink$V2, start = dflink$V4, end = dflink$V5)
dflink_to <- data.frame(chr = dflink$V7, start = dflink$V9, end = dflink$V10)
dflen <- read.table(args[1], sep = "\t", header = FALSE)
dflen <- data.frame(name = dflen$V1, start = rep(1, nrow(dflen)), end = dflen$V2)

# initialize circos
pdf(file = paste0(args[2], ".circle.pdf"))
circos.par("track.height" = 0.1, start.degree = 90)
circos.genomicInitialize(dflen, sector.names = dflen$name, plotType = c("labels"))

# draw chromosome
circos.track(ylim = c(0, 0.05), panel.fun = function(x, y) {
			 	xlim = CELL_META$xlim
				ylim = CELL_META$ylim
				circos.rect(xleft = xlim[1], ybottom = 0, xright = xlim[2], ytop = 0.05, col = "grey")
			},
			track.height = 0.05, bg.border = NA)

# draw colinear link
color <- data.frame(group = rep("", nrow(dflink)), code = rep("", nrow(dflink)))
for (i in 1:nrow(dflink)) {
  color$group[i] = paste(dflink_from$chr[i], "-", dflink_to$chr[i])
}
temp <- data.frame(group = levels(factor(color$group)), code = hue_pal(l = 75, c = 100)(length(levels(factor(color$group)))))
for (i in 1:nrow(dflink)) {
  color$code[i] <- temp[temp$group == color$group[i], 'code']
}
circos.genomicLink(region1 = dflink_from, region2 = dflink_to, col = color$code)
circos.clear()
dev.off()


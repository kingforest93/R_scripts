#!/usr/bin/env Rscript

# load package
suppressPackageStartupMessages(library(circlize))
library(RColorBrewer)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
	cat("Usage: Rscript circos_plot.R H.tub_genome_chr.fa.len H.tub_genome_chr.genes.gff.mRNA.1M.stat H.tub_genome_chr.TEs.gff.1M.stat H.tub_genome_chr.TRs.gff.1M.stat H.tub_genome_chr.fa.GC.1M.stat\n")
	quit()
}

# read chromosome.len and feature density
chrlen <- read.table(args[1], sep = "\t", header = FALSE)
chrlen <- data.frame(chr = chrlen$V1, start = rep(1, nrow(chrlen)), end = chrlen$V2)
genebed <- read.table(args[2], sep = "\t", header = FALSE)
genebed <- data.frame(chr = genebed[,1], start = genebed[,3], end = genebed[,4], value = genebed[,5])
tebed <- read.table(args[3], sep = "\t", header = FALSE)
tebed <- data.frame(chr = tebed[,1], start = tebed[,3], end = tebed[,4], value = tebed[,6])
trbed <- read.table(args[4], sep = "\t", header = FALSE)
trbed <- data.frame(chr = trbed[,1], start = trbed[,3], end = trbed[,4], value = trbed[,5])
gcbed <- read.table(args[5], sep = "\t", header = FALSE)
colnames(gcbed) <- c("chr", "start", "end", "value")

# initialize circos
pdf(file = paste0(args[1], ".mRNA.TE.TR.GC.circos.pdf"))
circos.par(track.height = 0.1, start.degree = 87, gap.degree = c(rep(0.4, nrow(chrlen) - 1), 3), cell.padding = c(0, 0, 0, 0))
colors = brewer.pal(8, "Set2")

# draw chromosome track
circos.genomicInitialize(chrlen, sector.names = chrlen$chr, plotType = c("labels", "axis"), major.by = 50000000, axis.labels.cex = 0.25*par("cex"), labels.cex = 0.5*par("cex"))
circos.track(ylim = c(0, 1), bg.col = colors[8], bg.border = "black", track.height = 0.025)

# draw gene density track
circos.genomicTrack(genebed,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area = TRUE, type = "line", lwd = 0.1, col = colors[1], border = colors[1])
}, ylim = c(0, 200), track.height = 0.05, bg.col = NA, bg.border = "black")

# draw TE density track
circos.genomicTrack(tebed,
	panel.fun = function(region, value, ...) {
		circos.genomicLines(region, value, area = TRUE, type = "line", lwd = 0.1, col = colors[2], border = colors[2])
}, ylim = c(0, 1200000), track.height = 0.05, bg.col = NA, bg.border = "black")

# draw TR density track
circos.genomicTrack(trbed,
	panel.fun = function(region, value, ...) {
		circos.genomicLines(region, value, area = TRUE, type = "line", lwd = 0.1, col = colors[3], border = colors[3])
}, ylim = c(0, 600), track.height = 0.05, bg.col = NA, bg.border = "black")

# draw GC% track
circos.genomicTrack(gcbed,
	panel.fun = function(region, value, ...) {
		circos.genomicLines(region, value, area = TRUE, type = "line", lwd = 0.1, col = colors[4], border = colors[4])
}, ylim = c(30, 50), track.height = 0.05, bg.col = NA, bg.border = "black")

# reset, save and close
circos.clear()
dev.off()


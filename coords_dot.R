#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
    cat("Usage: Rscript coords_dot.R mummer.delta.1v1.coords.tab minimum_alignment_length(1) minimum_identity(1) plot_size(10)\n")
    quit()
}
min_aln = args[2]
min_iden = args[3]
plot_size = args[4]
if (is.na(min_aln)) {min_aln = 100}
if (is.na(min_iden)) {min_iden = 50}
if (is.na(plot_size)) {plot_size = 10}

# read mummer.delta.1v1.coords.tab and filter
dt <- read.table(args[1], sep = "\t", header = FALSE)
dt <- dt[, -c(10, 11)]
colnames(dt) <- c('tstart', 'tend', 'qstart', 'qend', 'taln', 'qaln', 'iden', 'tlen', 'qlen', 'tname', 'qname')
dt <- subset(dt, taln >= min_aln & qaln >= min_aln & iden >= min_iden)

# convert to 2D coordinate system
newq <- dt[, c('qname', 'qlen')]
newq <- newq[!duplicated(newq$qname), ]
if (nrow(newq) > 1) {
	newq <- newq[order(newq$qlen, decreasing = TRUE), ]
	newq <- data.frame(qname = newq$qname, qlen = cumsum(newq$qlen), qadd = c(0, cumsum(newq$qlen)[-nrow(newq)]))
	for (i in 1:nrow(dt)) {
		dt$qstart[i] = dt$qstart[i] + newq[newq$qname == dt$qname[i], 'qadd']
		dt$qend[i] = dt$qend[i] + newq[newq$qname == dt$qname[i], 'qadd']
	}
}

newt <- dt[, c('tname', 'tlen')]
newt <- newt[!duplicated(newt$tname), ]
if (nrow(newt) > 1) {
	newt <- newt[order(newt$tlen, decreasing = TRUE), ]
	newt <- data.frame(tname = newt$tname, tlen = cumsum(newt$tlen), tadd = c(0, cumsum(newt$tlen)[-nrow(newt)]))
	for (i in 1:nrow(dt)) {
		dt$tstart[i] = dt$tstart[i] + newt[newt$tname == dt$tname[i], 'tadd']
		dt$tend[i] = dt$tend[i] + newt[newt$tname == dt$tname[i], 'tadd']
	}
}

# draw dot plot
xpos <- rep(0, nrow(newt))
xpos[1] <- newt$tlen[1] / 2
if (length(xpos) > 1) {
	for (i in 2:length(xpos)) {
  		xpos[i] = (newt$tlen[i] + newt$tlen[i - 1]) / 2
	}
}

ypos <- rep(0, nrow(newq))
ypos[1] <- newq$qlen[1] / 2
if (length(ypos) > 1) {
	for (i in 2:length(ypos)) {
  		ypos[i] = (newq$qlen[i] + newq$qlen[i - 1]) / 2
	}
}

gp <- ggplot(data = dt) +
	geom_segment(aes(x = tstart, xend = tend, y = qstart, yend = qend, color = iden), size = 1) +
	geom_hline(data = newq, aes(yintercept = qlen), color = "grey", linetype = 5) +
	geom_vline(data = newt, aes(xintercept = tlen), color = "grey", linetype = 5) +
	scale_x_continuous(breaks = xpos, labels = newt$tname, expand = c(0, 0)) +
	scale_y_continuous(breaks = ypos, labels = newq$qname, expand = c(0, 0)) +
	scale_color_gradient(low = "green", high = "red") +
	xlab(NULL) + ylab(NULL) +
	theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.ticks = element_line(color = "white"))
ggsave(filename = paste0(args[1], ".dotplot.pdf"), plot = gp, device = "pdf", width = plot_size * 1.2, height = plot_size)


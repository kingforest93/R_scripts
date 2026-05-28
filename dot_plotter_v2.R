#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[2]) ) {
    cat("Usage: Rscript dot_plotter.R Htub.Chr03_04.len Htub.Chr03_04.collinearity.filtered2.list FEH.ids\n")
    quit()
}

# read chromosome.len, collinearity.list and gene.ids
dflink <- read.table(args[2], sep = "\t", header = FALSE)
if (substr(dflink$V2[1], 1, 4) == substr(dflink$V7[1], 1, 4)) {
	dflink <- dflink[, c(1,2,4,5,6,7,9,10)]
} else {
	dflink <- dflink[as.integer(row.names(dflink)) %% 2 == 0, c(1,2,4,5,6,7,9,10)]	
}
dflen <- read.table(args[1], sep = "\t", header = FALSE)
dflen <- data.frame(name = dflen$V1, len = dflen$V2)
ids <- read.table(args[3], sep = "\t", header = FALSE)
ids <- c(ids$V1)

# convert to 2D coordinate system
dflink <- data.frame(qid = dflink$V1, qname = dflink$V2, qlen = rep(0, nrow(dflink)), qstart = dflink$V4, qend = dflink$V5, 
                     tid = dflink$V6, tname = dflink$V7, tlen = rep(0, nrow(dflink)), tstart = dflink$V9, tend = dflink$V10)

newq <- dflink[, c('qname', 'qlen')]
newq <- newq[!duplicated(newq$qname), ]
rownames(newq) <- 1:nrow(newq)
for (i in 1:nrow(newq)) {
	newq$qlen[i] = dflen$len[which(dflen$name == newq$qname[i])]
}
if (nrow(newq) > 1) {
	newq <- newq[order(newq$qname, decreasing = FALSE), ]
	newq <- data.frame(qname = newq$qname, qlen = cumsum(newq$qlen), qadd = c(0, cumsum(newq$qlen)[-nrow(newq)]))
}

newt <- dflink[, c('tname', 'tlen')]
newt <- newt[!duplicated(newt$tname), ]
rownames(newt) <- 1:nrow(newt)
for (i in 1:nrow(newt)) {
	newt$tlen[i] = dflen$len[which(dflen$name == newt$tname[i])]
}
if (nrow(newt) > 1) {
	newt <- newt[order(newt$tname, decreasing = FALSE), ]
	newt <- data.frame(tname = newt$tname, tlen = cumsum(newt$tlen), tadd = c(0, cumsum(newt$tlen)[-nrow(newt)]))
}

if (nrow(newq) > 1) {
	for (i in 1:nrow(dflink)) {
  		dflink$qstart[i] = dflink$qstart[i] + newq[newq$qname == dflink$qname[i], 'qadd']
  		dflink$qend[i] = dflink$qend[i] + newq[newq$qname == dflink$qname[i], 'qadd']
	}
}

if (nrow(newt) > 1) {
	for (i in 1:nrow(dflink)) {
		dflink$tstart[i] = dflink$tstart[i] + newt[newt$tname == dflink$tname[i], 'tadd']
		dflink$tend[i] = dflink$tend[i] + newt[newt$tname == dflink$tname[i], 'tadd']
	}
}

# mark selected genes
gene <- data.frame(qid = c(""), qname = c(""), qlen = c(0), qstart = c(0), qend = c(0),
                   tid = c(""), tname = c(""), tlen = c(0), tstart = c(0), tend = c(0))
for (i in 1:nrow(dflink)) {
	if (dflink$qid[i] %in% ids | dflink$tid[i] %in% ids) {
		gene <- rbind(gene, dflink[i,])
	}
}
gene <- gene[-1,]

# draw dot plot
xpos <- rep(0, nrow(newq))
xpos[1] <- newq$qlen[1] / 2
if (nrow(newq) > 1) {
	for (i in 2:length(xpos)) {
  		xpos[i] = (newq$qlen[i] + newq$qlen[i - 1]) / 2
	}
}

ypos <- rep(0, nrow(newt))
ypos[1] <- newt$tlen[1] / 2
if (nrow(newt) > 1) {
	for (i in 2:length(ypos)) {
  		ypos[i] = (newt$tlen[i] + newt$tlen[i - 1]) / 2
	}
}

gp <- ggplot(data = dflink) +
	geom_point(aes(x = qstart, y = tstart, color = factor(paste0(qname, tname))), size = 1) +
  	#geom_point(aes(x = qend, y = tend, color = factor(paste0(qname, tname))), size = 1) +
	#geom_segment(aes(x = qstart, xend = qend, y = tstart, yend = tend, color = factor(paste0(qname, tname))), size = 1) +
	#geom_segment(data = gene, aes(x = qstart, xend = qend, y = tstart, yend = tend), color = "black", size = 1) +
	geom_point(data = gene, aes(x = qstart, y = tstart), color = "black", size = 2, shape = 1) +
	geom_hline(data = newt, aes(yintercept = tlen), colour = "grey", linetype = 5) +
	geom_vline(data = newq, aes(xintercept = qlen), colour = "grey", linetype = 5) +
	#geom_hline(data = gene, aes(yintercept = tstart), colour = "black", linetype = 5) +
	#geom_hline(data = gene, aes(yintercept = tend), colour = "grey", linetype = 1) +
	#geom_vline(data = gene, aes(xintercept = qstart), colour = "black", linetype = 5) +
	#geom_vline(data = gene, aes(xintercept = qend), colour = "grey", linetype = 1) +
	scale_x_continuous(breaks = xpos, labels = newq$qname, expand = c(0, 0)) +
	scale_y_continuous(breaks = ypos, labels = newt$tname, expand = c(0, 0)) +
	xlab(NULL) + ylab(NULL) +
	theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
      	axis.ticks = element_line(colour = "white"), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(filename = paste0(args[2], ".dotplot.pdf"), plot = gp, device = "pdf", width = 17, height = 16)

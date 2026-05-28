#!/usr/bin/env Rscript

# load package
library(ggplot2)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[1]) ) {
    cat("Usage: Rscript paf_aln_dot.R minimap2.aln.paf minimum_sequence_length(10000) minimum_alignment_length(100) minimum_identity(50) plot_size(10)\n")
    quit()
}
min_seq = args[2]
min_aln = args[3]
min_iden = args[4]
plot_size = args[5]
if (is.na(min_seq)) {min_seq = 10000}
if (is.na(min_aln)) {min_aln = 100}
if (is.na(min_iden)) {min_iden = 50}
if (is.na(plot_size)) {plot_size = 20}

cat(paste("Minimum sequence length:", min_seq, "\n"))
cat(paste("Minimum alignment length: ", min_aln, "\n"))
cat(paste("Minimum alignment identity: ", min_iden, "\n"))

# read minimap2.paf.tab, filter and flip alignments
dt <- read.table(args[1], sep = "\t", header = FALSE, fill = TRUE)
dt <- dt[, c(1:11)]
colnames(dt) <- c('qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'match', 'align')
#colnames(dt) <- c('qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'match', 'align', 'iden')
dt$iden <- round(dt$match / dt$align * 100, 2)
#dt$iden <- round((1 - dt$iden) * 100, 2)
cat(paste("Number of alignments before filtering: ", nrow(dt), "\n"))

dt <- subset(dt, qlen >= as.numeric(min_seq) & tlen >= as.numeric(min_seq))
dt <- subset(dt, align >= as.numeric(min_aln) & iden >= as.numeric(min_iden))
tem <- dt$qstart[which(dt$strand == '-')]
dt$qstart[which(dt$strand == '-')] <- dt$qend[which(dt$strand == '-')]
dt$qend[which(dt$strand == '-')] <- tem
cat(paste("Number of alignments after filtering: ", nrow(dt), "\n"))

# convert to 2D coordinate system
newq <- dt[, c('qname', 'qlen')]
newq <- newq[!duplicated(newq$qname), ]
if (nrow(newq) > 1) {
	newq <- newq[order(newq$qname), ]
	newq <- data.frame(qname = newq$qname, qlen = cumsum(as.numeric(newq$qlen)), qadd = c(0, cumsum(as.numeric(newq$qlen))[-nrow(newq)]))
	for (i in 1:nrow(newq)) {
		dt$qstart[which(dt$qname == newq$qname[i])] = dt$qstart[which(dt$qname == newq$qname[i])] + newq$qadd[i]
		dt$qend[which(dt$qname == newq$qname[i])] = dt$qend[which(dt$qname == newq$qname[i])] + newq$qadd[i]
	}
}

newt <- dt[, c('tname', 'tlen')]
newt <- newt[!duplicated(newt$tname), ]
if (nrow(newt) > 1) {
	newt <- newt[order(newt$tname), ]
	newt <- data.frame(tname = newt$tname, tlen = cumsum(as.numeric(newt$tlen)), tadd = c(0, cumsum(as.numeric(newt$tlen))[-nrow(newt)]))
	for (i in 1:nrow(newt)) {
		dt$tstart[which(dt$tname == newt$tname[i])] = dt$tstart[which(dt$tname == newt$tname[i])] + newt$tadd[i]
		dt$tend[which(dt$tname == newt$tname[i])] = dt$tend[which(dt$tname == newt$tname[i])] + newt$tadd[i]
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
	geom_segment(aes(x = tstart, xend = tend, y = qstart, yend = qend, color = iden), linewidth = 1) +
	geom_hline(data = newq, aes(yintercept = qlen), color = "black", linetype = "dashed", linewidth = 0.3) +
	geom_vline(data = newt, aes(xintercept = tlen), color = "black", linetype = "dashed", linewidth = 0.3) +
	scale_x_continuous(breaks = xpos, labels = newt$tname, expand = c(0, 0)) +
	scale_y_continuous(breaks = ypos, labels = newq$qname, expand = c(0, 0)) +
	scale_color_gradient(low = "yellow", high = "red") +
	xlab(NULL) + ylab(NULL) +
	theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.ticks = element_line(color = "white"), axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(args[1], ".dotplot.pdf"), plot = gp, device = "pdf", width = plot_size * 1.2, height = plot_size)
cat("Ploting done!\n")

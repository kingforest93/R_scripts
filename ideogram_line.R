#!/usr/bin/env Rscript

# load package
library(RIdeogram)
library(RColorBrewer)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( length(args) < 2 ) {
    cat("Usage: Rscript ideogram_annotation.R G.cor_genome_chr.fa.len G.cor_genome_chr.genes.gff.mRNA.1M.stat G.cor_genome_chr.TEs.gff.intact_copia.gff.1M.stat G.cor_genome_chr.TEs.gff.intact_gypsy.gff.1M.stat\n")
    quit()
}

# read chromosome.len and annotation.gff
karydt <- read.table(args[1], sep = "\t", header = FALSE)
overdt <- read.table(args[2], sep = "\t", header = TRUE)
label1 <- read.table(args[3], sep = "\t", header = TRUE)
#label2 <- read.table(args[4], sep = "\t", header = TRUE)

#convert to karyotype and density object
karydt <- data.frame(Chr = karydt$V1, Start = rep(1, nrow(karydt)), End = karydt$V2)
overdt <- data.frame(Chr = overdt[,1], Start = overdt[,3], End = overdt[,4], Value = overdt[,5])
#labeldt <- data.frame(Chr = label1[,1], Start = label1[,3], End = label1[,4], Value_1 = label1[,6], Color_1 = rep("e78ac3", nrow(label1)))
labeldt <- data.frame(Chr = label1[,1], Start = label1[,3], End = label1[,4], Value_1 = label1[,6], Color_1 = rep("66c2a5", nrow(label1)))
#labeldt <- data.frame(Chr = label1[,1], Start = label1[,3], End = label1[,4], Value_1 = label1[,6], Color_1 = rep("ffd92f", nrow(label1)), Value_2 = label2[,6], Color_2 = rep("a6d854", nrow(label2)))
#labeldt <- data.frame(Chr = label1[,1], Start = label1[,3], End = label1[,4], Value_1 = label1[,6], Color_1 = rep("8da0cb", nrow(label1)), Value_2 = label2[,6], Color_2 = rep("fc8d62", nrow(label2)))

# draw ideogram and density plot
ideogram(karyotype = karydt, overlaid = overdt, label = labeldt, label_type = "line", colorset1 = c("#ffffff", "#000080"))
convertSVG("chromosome.svg", device = "pdf")
file.rename("chromosome.pdf", paste0(args[1], ".density.pdf"))
file.remove("chromosome.svg")


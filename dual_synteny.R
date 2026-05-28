#!/usr/bin/env Rscript

# load package
library(RIdeogram)
library(scales)

# parse input
args <- commandArgs(trailingOnly = TRUE)
if( is.na(args[2]) ) {
    cat("Usage: Rscript dual_synteny.R G.cor_vs_C.set_genome_chr.fa.len G.cor_vs_C.set_genome_chr.pep.collinearity.filtered.list\n")
    quit()
}

# read chromosome.len and collinearity.list
dflink <- read.table(args[2], sep = "\t", header = FALSE)
dflink <- dflink[as.integer(row.names(dflink)) %% 2 == 1, c(2,4,5,7,9,10)]
dflen <- read.table(args[1], sep = "\t", header = FALSE)

#convert to synteny and karyotype object
karydt <- data.frame(Chr = substring(dflen$V1, 6), Start = rep(1, nrow(dflen)), End = dflen$V2,
  					 fill = rep(substr(tolower(col2hcl("grey")), 2, 7), nrow(dflen)),
					 species = substring(dflen$V1, 1, 4), size = rep(12, nrow(dflen)),
					 color = rep(substr(tolower(col2hcl("black")), 2, 7), nrow(dflen)))

syndt <- data.frame(Species_1 = substring(dflink$V2, 6), Start_1 = dflink$V4, End_1 = dflink$V5,
                    Species_2 = substring(dflink$V7, 6), Start_2 = dflink$V9, End_2 = dflink$V10,
                    fill = rep("", nrow(dflink)))
gp <- factor(paste(syndt$Species_1, syndt$Species_2, sep = "-"))
tem <- data.frame(group = levels(gp), code = hue_pal(c = 100, l = 75)(length(levels(gp))))
sp1 <- karydt[karydt$species == substring(dflink$V2[1], 1, 4), 'Chr']
sp2 <- karydt[karydt$species == substring(dflink$V7[1], 1, 4), 'Chr']

for (i in 1:nrow(dflink)) {
	syndt$fill[i] = substr(tolower(tem[tem$group == paste(syndt$Species_1[i], syndt$Species_2[i], sep = "-"), 'code']), 2, 7)
	syndt$Species_1[i] <- which(sp1 == syndt$Species_1[i])
	syndt$Species_2[i] <- which(sp2 == syndt$Species_2[i])
}
syndt$Species_1 <- as.integer(syndt$Species_1)
syndt$Species_2 <- as.integer(syndt$Species_2)

# draw dual synteny plot
ideogram(karyotype = karydt, synteny = syndt)
convertSVG("chromosome.svg", device = "pdf")
file.rename("chromosome.pdf", paste0(args[2], ".synteny.pdf"))
file.remove("chromosome.svg")


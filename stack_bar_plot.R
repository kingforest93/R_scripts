#/usr/bin/env Rscript

#load package
library(ggplot2)
library(tidyr)
library(RColorBrewer)

#read data
args <- commandArgs(TRUE)
if (length(args) != 1) {
    cat("Rscript stack_bar_plot.R file.tsv\n")
    quit()
}

dat <- read.csv(args[1], header = TRUE, encoding = "UTF-8")
dat[,1] <- paste(dat[,1], dat[,2], sep = "_")
colnames(dat)[1] <- "Genome"
dat <- dat[order(dat$All),]
dat <- dat[,-c(2,9)]
genome <- dat$Genome
type <- colnames(dat)[-1]
dat <- gather(dat, key = "Type", value = "Perc", -Genome)

#draw stacked bar plot
dat$Genome <- factor(dat$Genome, levels = genome)
dat$Type <- factor(dat$Type, levels = rev(type))
bar.p <- ggplot(dat, aes(x = Genome, y = Perc)) +
  geom_bar(aes(fill = Type), position = "stack", stat = "identity",
           color = "black", width = 0.6, size = 0.2) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme(line=element_line(color="black", size=0.2),
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.2),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank())
ggsave(paste0(args[1], ".stacked.barplot.pdf"), bar.p, width = 14, height = 8, units = "cm")

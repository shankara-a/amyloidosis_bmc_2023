#!/usr/bin/Rscript

# -- import packages: --------------------------------------------------------------------
suppressMessages(library('ggpubr'))
suppressMessages(library('ggdendro'))
suppressMessages(library('dplyr'))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

source("funcs/rfuncs.R")

# Figure 2
pdf("figures/full_na_dataset_enrichments.pdf", w=10,h=8)
plotE(
    "data/clustering/full_na_dataset/contrasts_qvars.tsv",
    "data/clustering/full_na_dataset/contrasts_fe.tsv"
)
dev.off()

# Figure S3-4
p1 <- plotE1("data/clustering/full_na_dataset/subgroup_comparisons/contrasts_qvars_cardiac_stage_2.tsv")
p2 <- plotE1("data/clustering/full_na_dataset/subgroup_comparisons/contrasts_qvars_renal_stage_1.tsv")

pdf("figures/full_stage2_renal1.pdf", w=10,h=8)
options(repr.plot.width=10, repr.plot.height=8)
grid.arrange(p1,p2,nrow=1)
dev.off()
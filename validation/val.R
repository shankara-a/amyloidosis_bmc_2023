#!/usr/bin/Rscript

library(ConsensusClusterPlus)
library(mice)
library(cluster)
library(purrr)
library(factoextra)
library(clValid)

# Random Seed
set.seed(123)

# Load Data
X.mice <- read.table("../data/imputed/mice_qvars_05.tsv", sep="\t", header=T)
X.mice.scaled <- scale(complete(X.mice))

data.df <- read.table("../data/processed/AL_with_ccp_03.tsv", sep="\t", header=T, row.names=1)

# Run Internal Validation
i.val <- clValid(
        t(X.mice.scaled), 
        2:6, 
        clMethods=c("hierarchical","kmeans","pam","fanny","clara"),
        validation = c("internal","stability"),
        metric="correlation",
        maxitems=ncol(X.mice.scaled),
        method="average",
        verbose=TRUE
)

# Save
saveRDS(i.val, "results/01_2074_ival.rds")


#!/usr/bin/Rscript


suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(mice))
suppressMessages(library(cluster))
suppressMessages(library(purrr))
suppressMessages(library(factoextra))
suppressMessages(library(clValid))
suppressMessages(library(dplyr))
suppressMessages(library(gridExtra))
suppressMessages(library(usedist))
source("funcs/rfuncs.R")

# Params
n_iterations <- 50
n_samples <- 1750

#-------------------------------
# Load Data
#-------------------------------
data.df <- read.table("data/processed/AL_with_ccp_03.tsv", sep='\t', row.names=1, header=T)
data.df <- data.df[data.df$fna3_cluster %in% c(1,2,3),]

FILENAME <- "data/processed/AL_for_ccp_02.tsv"
OUT_DIR <- "data/clustering"

data <- getDistMatrix(FILENAME)

X <- data$Xfilt
dt <- data$dt

metrics_list = list()
i <- 1

clusters.df <- runClusterings(X, dt)
clusters.df$fna3_cluster <- data.df[rownames(X),"fna3_cluster"]

metrics_list[[i]] <- computeClusterMetrics(dt, clusters.df)
metrics_list[[i]]$run <- "full"
i <- i+1

for (n in 1:n_iterations) {
        print(paste("Running iteration:", n))
        idx <- sample(1:nrow(X), n_samples, replace = FALSE)

        clusters.df <- runClusterings(X[idx,], dist_subset(dt, idx))
        clusters.df$fna3_cluster <- data.df[rownames(X[idx,]),"fna3_cluster"]

        metrics_list[[i]] <- computeClusterMetrics(dist_subset(dt, idx), clusters.df)
        metrics_list[[i]]$run <- paste("sample", n, sep="_")
        i <- i+1
}

# Combine and save
metrics.df <- do.call(rbind, metrics_list)
file.name <- paste(OUT_DIR, "cluster_algorithm_comparison.tsv", sep="/")
write.table(metrics.df, file.name, sep="\t")      

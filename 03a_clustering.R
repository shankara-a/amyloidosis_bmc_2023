#!/usr/bin/Rscript

# -- import packages: --------------------------------------------------------------------
library(ConsensusClusterPlus)
library(mice)

# Clustering inputs
maxK <- 6
reps <- 100
pItem = 0.8
pFeature = 1
seed <- 123

# ------------------------------
# Imputed files used for clustering are found in data/imputed/**
# Loop through all datasets with different imputed values
# ------------------------------
for (filename in Sys.glob("data/imputed/**/*.tsv")){
    # Save filenames
    dataset.name <- strsplit(filename, "/")[[1]][3]
    imputed.name <- strsplit(strsplit(filename, "/")[[1]][4],".tsv")[[1]][1]
    out.dir <- paste("data/clustering", dataset.name, imputed.name, sep="/")

    print(paste("* Running for dataset: ", dataset.name, " for imputed data: ", imputed.name))
    print(paste("   Writing out to: ",out.dir))

    # Load data
    X <- read.table(filename, sep="\t", header=T)

    if ("Code.ID" %in% colnames(X)){
        rownames(X) <- X$Code.ID
        X$Code.ID <- NULL
    }
    
    # Scale data for clustering
    X.scaled <- scale(complete(X))

    # Write out input matrix
    write.table(X.scaled, paste(out.dir, "input_matrix.tsv", sep="/"), sep="\t")

    # Run clustering
    ccp = ConsensusClusterPlus(
        as.matrix(t(X.scaled)),
        maxK=maxK,
        reps=reps,
        pItem=pItem,
        pFeature=pFeature,
        clusterAlg="hc",
        distance="spearman",
        title=out.dir,
        seed=seed,
        plot="png"
    )

    # Save files
    icl = calcICL(ccp, title=out.dir, plot="png")
    write.table(icl[["itemConsensus"]], file.path(out.dir, "item_consensus.tsv"), sep='\t')
    write.table(icl[["clusterConsensus"]], file.path(out.dir, "cc_result.tsv"), sep='\t')
    saveRDS(ccp, file = file.path(out.dir, "ccp.rds"))
}

# ------------------------------
# Run consensus clustering w/ missing values
# ------------------------------
# Clustering inputs
maxK <- 6
reps <- 25
pItem = 0.8
pFeature = 1
seed <- 123

FILENAME <- "data/processed/AL_for_ccp_02.tsv"

# Load file
X <- read.table(FILENAME, sep="\t", header=T)
rownames(X) <- X$Code.ID
X$Code.ID <- NULL

# Scale data
X.scaled <- scale(complete(X))

# Test w/ missing data
out.dir <- "data/clustering/full_na_dataset"

# Compute distances with pairwise complete observations
# Note, we lose ~7 samples which we are unable to compare this way
dt <- as.dist(1-cor(as.matrix(t(X.scaled)), method="spearman", use="pairwise.complete.obs"))
x = as.matrix(dt)
x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
dt = as.dist(x)

# Save input
write.table(X[colnames(x),], paste(out.dir, "input_matrix.tsv", sep="/"), sep="\t")

# Run clustering
ccp = ConsensusClusterPlus(
    dt,
    maxK=maxK,
    reps=reps,
    pItem=pItem,
    pFeature=pFeature,
    clusterAlg="hc",
    distance="spearman",
    title=out.dir,
    seed=seed,
    plot="png",
    corUse = "complete.obs",
    verbose=F
)

# Save results
icl = calcICL(ccp, title=out.dir, plot="png")
write.table(icl[["itemConsensus"]], file.path(out.dir, "item_consensus.tsv"), sep='\t')
write.table(icl[["clusterConsensus"]], file.path(out.dir, "cc_result.tsv"), sep='\t')
saveRDS(ccp, file = file.path(out.dir, "ccp.rds"))
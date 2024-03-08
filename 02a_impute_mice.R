#!/usr/bin/Rscript

# -- import packages: --------------------------------------------------------------------
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(forestplot))
suppressMessages(library(ggsci))
suppressMessages(library(mice))

# -------------------------
# Full dataset
# -------------------------
X <- read.table('data/processed/AL_for_ccp_02.tsv', sep='\t', header=T, row.names='Code.ID')
Ximp <- mice(data = X, m = 5, maxit = 100, seed = 500)

# Save object
saveRDS(Ximp, file = "data/imputed/full_dataset/mice_qvars_01.rds")

for (i in 1:5){
    # Create new dataframe
    X.imputed <- complete(Ximp,i)
    rownames(X.imputed) <- rownames(X)
    
    write.table(X.imputed, paste("data/imputed/full_dataset/mice_qvars_0",i,".tsv", sep=""), sep='\t')
}

# -------------------------
# 2004 dataset
# -------------------------
X <- read.table('data/processed/AL_2004_for_ccp_02.tsv', sep='\t', header=T, row.names='Code.ID')
Ximp <- mice(data = X, m = 5, maxit = 100, seed = 500)

# Save object
saveRDS(Ximp, file = "data/imputed/2004_dataset/mice_qvars_01.rds")

for (i in 1:5){
    # Create new dataframe
    X.imputed <- complete(Ximp,i)
    rownames(X.imputed) <- rownames(X)
    
    write.table(X.imputed, paste("data/imputed/2004_dataset/mice_qvars_0",i,".tsv", sep=""), sep='\t')
}
# -------------------------
# 2008 dataset
# -------------------------
X <- read.table('data/processed/AL_2008_for_ccp_02.tsv', sep='\t', header=T, row.names='Code.ID')
Ximp <- mice(data = X, m = 5, maxit = 100, seed = 500)

# Save object
saveRDS(Ximp, file = "data/imputed/2008_dataset/mice_qvars_01.rds")

for (i in 1:5){
    # Create new dataframe
    X.imputed <- complete(Ximp,i)
    rownames(X.imputed) <- rownames(X)
    
    write.table(X.imputed, paste("data/imputed/2008_dataset/mice_qvars_0",i,".tsv", sep=""), sep='\t')
}
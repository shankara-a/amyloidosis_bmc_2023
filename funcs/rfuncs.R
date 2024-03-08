suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(forestplot))
suppressMessages(library(ggsci))
suppressMessages(library(dplyr))
suppressMessages(library(adjustedCurves))
suppressMessages(library(tidyverse))
suppressMessages(library(cmprsk))
suppressMessages(library(caret))
suppressMessages(library(ggplot2))
suppressMessages(library(ggalluvial))
suppressMessages(library(gridExtra))

getDistMatrix <- function(file.name, clustering_metrics=FALSE){
    # Load file
    X <- read.table(file.name, sep="\t", header=T)
    rownames(X) <- X$Code.ID
    X$Code.ID <- NULL

    # Scale data
    X.scaled <- scale(complete(X))

    # Compute distances with pairwise complete observations
    dt <- as.dist(1-cor(as.matrix(t(X.scaled)), method="spearman", use="pairwise.complete.obs"))
    x = as.matrix(dt)
    x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
    dt = as.dist(x)

    if (clustering_metrics==TRUE){
            # Plots
            p1 <- fviz_nbclust(X.scaled[rownames(x),], FUN = cluster::pam, diss=dt, method = "silhouette") + labs(title='PAM')
            p2 <- fviz_nbclust(X.scaled[rownames(x),], FUN = hcut, diss=dt, method = "silhouette") + labs(title='Hierarchical')
            p3 <- fviz_nbclust(X.scaled[rownames(x),], FUN = cluster::fanny, diss=dt, method = "silhouette") + labs(title='Fanny')
            p4 <- fviz_nbclust(X.scaled[rownames(x),], FUN = cluster::clara, diss=dt, method = "silhouette", correct.d=T) + labs(title='CLARA')

            p11 <- fviz_nbclust(X.scaled[rownames(x),], FUN = cluster::pam, diss=dt, method = "wss") + labs(title='PAM')
            p21 <- fviz_nbclust(X.scaled[rownames(x),], FUN = hcut, diss=dt, method = "wss") + labs(title='Hierarchical')
            p31 <- fviz_nbclust(X.scaled[rownames(x),], FUN = cluster::fanny, diss=dt, method = "wss") + labs(title='Fanny')
            p41 <- fviz_nbclust(X.scaled[rownames(x),], FUN = cluster::clara, diss=dt, method = "wss", correct.d=T) + labs(title='CLARA')

            options(repr.plot.width=16, repr.plot.height=6)
            p <- grid.arrange(p1,p2,p3,p4,p11,p21,p31,p41, nrow=2)

            return(list("X" = X.scaled, "dt" = dt, "Xfilt" = X.scaled[rownames(x),], "plot"=p))
    } else{
        return(list("X" = X.scaled, "dt" = dt, "Xfilt" = X.scaled[rownames(x),]))
    }
}

getCoef <- function(mod){
    df <- data.frame(summary(mod)$coef)
    df$coef_exp <- signif(exp(df$coef),3)
    df$ci_lower <- signif(exp(df$coef - 1.96 * df$se.coef.),3)
    df$ci_upper <- signif(exp(df$coef + 1.96 * df$se.coef.),3)
    df$HR <- paste(df$coef_exp, " (", df$ci_lower, "-", df$ci_upper, ")", sep="")
    return(data.frame(df))
}

coxSummary <- function(cox_models, var, filt=TRUE){
    i <- 1
    result_list = list()

    for (i in seq_along(cox_models)) {
        formula <- Reduce(paste,deparse(formula(cox_models[[i]])))
        hazard_ratios <- getCoef(cox_models[[i]])$HR
        names(hazard_ratios) <- rownames(getCoef(cox_models[[i]]))

        if (filt==T){
            hazard_ratios <- hazard_ratios[names(hazard_ratios)[grepl(paste0("^", var), names(hazard_ratios))]]
            names(hazard_ratios) <- sub(paste0("^", var), "", names(hazard_ratios))
        }

        df <- data.frame(rbind(hazard_ratios))
        df$Model <- as.character(formula)
        result_list[[i]] <- df
    }
    # Combine all dataframes into one
    results_df <- do.call(rbind, result_list)
    rownames(results_df) <- results_df$Model
    results_df$Model <- NULL

    # Return the dataframe
    return(results_df)

}

getRiskDiff <- function(clust.adjsurv, groupings, times.to.use=c(1,5,10)){
    # Risk Difference Table
    clust.adjci <- clust.adjsurv
    clust.adjci$adjsurv$surv <- 1-clust.adjci$adjsurv$surv

    # Define the desired order of the groups
    order <- groupings

    # Print the reordered dataframe
    result_list = list()

    # For each grouping
    for (i in seq_along(groupings)){
        group_1 = str_split(groupings[[i]], " vs. ")[[1]][1]
        group_2 = str_split(groupings[[i]], " vs. ")[[1]][2]
        result_list[[i]] <- adjusted_curve_diff(clust.adjci, times=times.to.use, group_1=group_1, group_2=group_2, conf_int=TRUE)
        result_list[[i]]$group <- groupings[[i]]
    }

    cd.df <- do.call(rbind, result_list)

    # Percentage
    cd.df$diff <- signif(cd.df$diff*100,3)
    cd.df$ci_lower <- signif(cd.df$ci_lower*100,3)
    cd.df$ci_upper <- signif(cd.df$ci_upper*100,3)
    cd.df$RD <- paste(cd.df$diff, " (", cd.df$ci_lower, "-", cd.df$ci_upper, ")", sep="")
    
    r <- spread(cd.df[,c("time","RD","group")], key = time, value = RD) %>% 
        mutate(group = factor(group, levels = order)) %>%
        arrange(group)

    return(r)
}

getIncidenceGroup <- function(clust.adjsurv, times.to.use=c(1,5,10)){
  clust.adjci <- clust.adjsurv
  clust.adjci$adjsurv$surv <- 1-clust.adjci$adjsurv$surv

  # CI Difference Table
  closest_rows <- list()

  for (i in seq_along(times.to.use)){
    closest_rows[i] <- clust.adjsurv$boot_adjsurv[which.min(abs(clust.adjsurv$boot_adjsurv$time - times.to.use[[i]])), "time"]
  }
  ci.df <- clust.adjsurv$boot_adjsurv[clust.adjsurv$boot_adjsurv$time %in% closest_rows,]
  ci.df$boot_surv <- signif((1-ci.df$boot_surv)*100,3)
  ci.df$ci_upper <- signif((1-ci.df$ci_upper)*100,3)
  ci.df$ci_lower <- signif((1-ci.df$ci_lower)*100,3)
  ci.df$CI <- paste(ci.df$boot_surv, " (", ci.df$ci_upper, "-", ci.df$ci_lower, ")", sep="")

  result.df <- spread(ci.df[,c("time","CI","group")], key = time, value = CI) %>% 
    mutate(group = factor(group)) %>%
    arrange(group)
  
  rownames(result.df) <- result.df$group
  result.df$group <- NULL
  colnames(result.df) <- times.to.use
  return(result.df)
}

# ----------------------------------------
# Competing Risks
# ----------------------------------------
getDummy <- function(df, vars){
    dmy <- dummyVars(" ~ .", data = df[,vars,drop=F], fullRank=TRUE)
    cov <- data.frame(predict(dmy, newdata = df[,vars,drop=F]))
    return(cov)
}

runFineGray <- function(df, covariates){
    cov1 <- getDummy(df, covariates)

    eskd.crr <- crr(ftime=df$CR_time, fstatus=df$CR_event, cov1=cov1, failcode=1, cencode=0)
    death.crr <- crr(ftime=df$CR_time, fstatus=df$CR_event, cov1=cov1, failcode=2, cencode=0)

    dcr.df <- getCoef(death.crr)
    dcr.df$Risk <- "Death"
    ecr.df <- getCoef(eskd.crr)
    ecr.df$Risk <- "ESKD"

    dcr.df <- tibble::rownames_to_column(dcr.df, "covariate")
    ecr.df <- tibble::rownames_to_column(ecr.df, "covariate")
    
    coef.df <- rbind(dcr.df,ecr.df)
    coef.df$formula <- as.character(paste("~ ",paste(covariates, collapse=" + ")))

    result <- list("death.crr" = death.crr, "eskd.crr" = eskd.crr, "coef"=coef.df)

    return(result)
}

# ----------------------------------------
# Plotting
# ----------------------------------------
plotSurv <- function(fit, data.df, legend.title, legend.labs, ...){
    p <- ggsurvplot(
        fit,                    
        data = data.df,
        risk.table = TRUE,
        pval = TRUE,
        conf.int = F,
        xlab = "Time (Yr)",
        risk.table.y.text.col = T,
        risk.table.y.text = F,
        legend.title = legend.title,
        legend.labs = legend.labs,
        ...
    )
    
    return(p)
}

plotAlluvial <- function(df, a, b, title="", colors=NULL){
    if (!a %in% colnames(df) || !b %in% colnames(df)) {
        stop("Cluster column names not found in dataframe.")
    }

    # Filter out missing values
    data.df <- df[,c(a,b)]
    data.df <- na.omit(data.df)
        
    # Create alluvial plot
    p <- ggplot(data.df, aes(axis1 = !!sym(a), axis2 = !!sym(b))) +
        geom_alluvium(aes(fill = !!sym(a))) +
        geom_stratum() +
        geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        theme_minimal() +
        labs(title = title, y = "Patients") +
        theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) + 
            #panel.grid.major.y = element_line(color = "transparent")) + 
        scale_fill_manual(values = colors)

    return(p)
}

plotEnrichment <- function(contrasts.df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black', fix_id=T){
    contrasts.df$sig <- contrasts.df$fdr_bh < pval.thresh
    contrasts.df$logq <- -log10(contrasts.df$fdr_bh)
    ### Order axis by dendrogram
    # Load data
    X <- contrasts.df[,c('X','cluster','statistic')]
    X <- reshape(X[,c('X','cluster','statistic')], timevar='cluster', idvar='X', direction='wide',)
    rownames(X) <- X$X
    X$X <- NULL

    X[is.na(X)] <- 0

    # Build the dendrogram
    dend <- as.dendrogram(hclust(d = dist(x = X)))
    dendro.plot <- ggdendrogram(dend,rotate = TRUE)

    # Use dendrogram order to order colomn
    order <- order.dendrogram(dend) # dendrogram order
    contrasts.df$X <- factor(x = contrasts.df$X, levels = unique(contrasts.df$X)[order], ordered = TRUE)

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        contrasts.df,
        x="cluster",
        y="X",
        fill = "statistic",
        size="logq",
        color=ifelse(contrasts.df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(contrasts.df$statistic)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=0))

    return(p)
}

plotFisherExactEnrichment <- function(contrasts.df, pval.thresh=0.1, filter=NULL, palette='Blues', h=13, w=15, s_color='black', fix_id=T){
    contrasts.df$sig <- contrasts.df$pval_adj < pval.thresh
    contrasts.df$logq <- -log10(contrasts.df$pval_adj)
    
    # Remove non-significants
    contrasts.df <- contrasts.df[contrasts.df$feat %in% contrasts.df[contrasts.df$sig,]$feat,]

    ### Order axis by dendrogram
    # Load data
    X <- contrasts.df[,c('feat','cluster','odds_r')]
    X <- reshape(X[,c('feat','cluster','odds_r')], timevar='cluster', idvar='feat', direction='wide',)
    rownames(X) <- X$feat
    X$feat <- NULL

    X[is.na(X)] <- 0

    # Build the dendrogram
    dend <- as.dendrogram(hclust(d = dist(x = X)))
    dendro.plot <- ggdendrogram(dend,rotate = TRUE)

    # Use dendrogram order to order colomn
    order <- order.dendrogram(dend) # dendrogram order
    contrasts.df$feat <- factor(x = contrasts.df$feat, levels = unique(contrasts.df$feat)[order], ordered = TRUE)

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        contrasts.df,
        x="cluster",
        y="feat",
        fill = "odds_r",
        size="logq",
        color=ifelse(contrasts.df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(contrasts.df$odds_r)) * c(0, 1))+
        labs(x="", y="", fill="log Odds Ratio", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=0))

    return(p)
}

plotE <- function(f1, f2){
    # Continuous
    contrasts.df <- read.table(f1, sep="\t", header=T)
    contrasts.df$cluster <- factor(contrasts.df$fna3_cluster_n, levels=c("Low","Intermediate","High"), ordered=T)

    # Fisher Exact
    contrasts.fe.df <- read.table(f2, sep="\t", header=T)
    contrasts.fe.df$odds_r <- -log10(contrasts.fe.df$odds_r)
    contrasts.fe.df <- contrasts.fe.df %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x))
    contrasts.fe.df$cluster <- factor(contrasts.fe.df$fna3_cluster_n, levels=c("Low","Intermediate","High"), ordered=T)

    # Create plots
    p1 <- plotFisherExactEnrichment(contrasts.fe.df)
    p2 <- plotEnrichment(contrasts.df)

    options(repr.plot.width=10, repr.plot.height=8)
    return(grid.arrange(p1, p2, nrow=1))
}

plotE1 <- function(f1){
    # Continuous
    contrasts.df <- read.table(f1, sep="\t", header=T)
    contrasts.df$cluster <- factor(contrasts.df$fna3_cluster_n, levels=c("Low","Intermediate","High"), ordered=T)
    p <- plotEnrichment(contrasts.df) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))

    options(repr.plot.width=5, repr.plot.height=8)
    return(p)
}

# ----------------------------------------
# Clustering
# ----------------------------------------
computeClusterMetrics <- function(dt, clusters.df){
    sil_list = list()
    connect_list = list()
    dunn_list = list()
    i=1
    j=1
    k=1

    for (col in names(clusters.df)){
        
        try({
            # Silhouette score
            sil_list[[i]] <- mean(silhouette(clusters.df[,col], dt)[,3])
            names(sil_list)[j] <- col
            i <- i+1

            # Connectivity
            connect_list[[j]] <- connectivity(clusters.df[,col], distance=dt)
            names(connect_list)[j] <- col
            j <- j+1

            # Dunn
            dunn_list[[k]] <- dunn(clusters.df[,col], distance=dt)
            names(dunn_list)[k] <- col
            k <- k+1
        })
    }

    metrics.df <- data.frame(cbind(cbind(t(data.frame(sil_list)),t(data.frame(connect_list))), t(data.frame(dunn_list))))
    names(metrics.df) <- c("Silhouette","Connectivity","Dunn")

    return(metrics.df)
}

runClusterings <- function(X, dt, ks=c(2,3,4,5,6)){
    # Clusterings with missing valuess
    result_list = list()
    i <- 1

    for (k in ks){
        # Heirarchical Clustering
        df <- data.frame(cutree(hclust(dt, method="average"),k))
        names(df)[1] <- paste("hc",k, sep="_")
        result_list[[i]]  <- df
        i <- i+1

        # Partition Around Medoids
        df <- data.frame(pam(dt,k)$cluster)
        names(df)[1] <- paste("pam",k, sep="_")

        result_list[[i]]  <- df
        i <- i+1

        # Fanny
        df <- data.frame(fanny(dt, k, memb.exp=1, maxit=1000)$cluster)
        names(df)[1] <- paste("fanny",k, sep="_")

        result_list[[i]]  <- df
        i <- i+1

        # Clara
        df <- data.frame(clara(X, k, correct.d=T)$cluster)
        names(df)[1] <- paste("clara",k, sep="_")

        result_list[[i]]  <- df
        i <- i+1
    }

    return(do.call(cbind, result_list))
}

computeClusterMetrics <- function(dt, clusters.df){
    sil_list = list()
    connect_list = list()
    i=1
    j=1

    for (col in names(clusters.df)){
        
        try({
            # Silhouette score
            sil_list[[i]] <- mean(silhouette(clusters.df[,col], dt)[,3])
            names(sil_list)[j] <- col
            i <- i+1

            # Connectivity
            connect_list[[j]] <- connectivity(clusters.df[,col], distance=dt)
            names(connect_list)[j] <- col
            j <- j+1
        })
    }

    metrics.df <- data.frame(cbind(cbind(t(data.frame(sil_list)),t(data.frame(connect_list)))))
    names(metrics.df) <- c("Silhouette","Connectivity")

    return(metrics.df)
}
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(forestplot))
suppressMessages(library(ggsci))
suppressMessages(library(dplyr))
suppressMessages(library(adjustedCurves))
suppressMessages(library(tidyverse))
suppressMessages(library(cmprsk))
suppressMessages(library(caret))

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
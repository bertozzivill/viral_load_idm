####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Find summary measures (coefficients, SE, AIC, Rubin error, RMSE[oos]) for all models tested

## Input:  main_dir: Directory in which imputed model data live
##                    
## Output: 
##
## Run from within the "vectorize" folder!
####################################################################################################


rm(list=ls())

library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)
library(stringr)

source("LinearSurvivalModel.R")
source("rubin.method.R")
source("output.AIC.R")
source("summarize_models.r")

#automatically finds correct directory name so we can stop commenting and uncommenting lines
dir.exists <- function(d) {
  de <- file.info(d)$isdir
  ifelse(is.na(de), FALSE, de)
}
root_dir <- ifelse(dir.exists("/home/cselinger/"), "/home/cselinger/HIV-Cascade/merge/data/", "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/")
main_dir <- ifelse(length(commandArgs())>2, commandArgs()[3], root_dir)


##---------------------------
## Load data
##---------------------------

#load imputed data and indices for data transforms and survial models
print("loading imputed data")
load(paste0(main_dir, "imputed_survival_data.rdata"))
load(paste0(main_dir, "indexes.rdata"))
load(paste0(main_dir, "survival_model_output.rdata"))

# fill a named list with model outputs
imputation_count <- 10
model_spec_names <- apply(index.survival.models,1,function(x)paste(x, collapse="-"))

survival.model.summaries<- NULL
list_idx <- 1

for (transform_imp in names(survival.model.output)){
  print(transform_imp)
  
  #find data transform and imputation county
  pattern <- "(.*)-imp_count_([0-9]*)"
  pattern_match <- str_match(transform_imp, pattern)
  data_transform <- pattern_match[1,2]
  imp_number <- as.numeric(pattern_match[1,3])

  this_transform<-survival.model.output[[transform_imp]]
  
  for (model_spec in colnames(this_transform)){
    this_model <- this_transform[, model_spec]
    model_summary <- summary(this_model$lm)$coefficients
    
    model_summary <- data.table(model_summary, 
                                 covariate=rownames(model_summary),
                                 imputation=imp_number,
                                 model_spec=model_spec,
                                 data_transform=data_transform)
    survival.model.summaries[[list_idx]] <- model_summary
    list_idx<- list_idx+1
    
  }
  
}

survival.model.summaries <- rbindlist(survival.model.summaries)

##-------------------------------
## Compute summary coefficients
##-------------------------------

#get summary coefficients for model outputs (across imputations)
by_names <- c("data_transform", "model_spec","covariate")

setnames(survival.model.summaries, c("Estimate", "Std. Error"), c("beta", "se"))

# means 
mean.survival.model.summaries <- survival.model.summaries[, list(beta=mean(beta)), by=by_names] 

#calculate standard errors using the method outlined on page 6 of the AMELIA documentation: https://cran.r-project.org/web/packages/Amelia/vignettes/amelia.pdf

#calculate variance of the mean estimates
survival.model.summaries[, across_beta_var:= var(beta), by=by_names]

#calculate mean of the variances
survival.model.summaries[, mean_variance:= mean(se^2), by=by_names]

#calculate joint variance of imputations
variance.survival.model.summaries<-unique(survival.model.summaries[, list(var=mean_variance + across_beta_var*(1+1/imputation_count)), by=by_names])
variance.survival.model.summaries[, se:=sqrt(var)]

#merge means and variances to get a full summary of results
summary.survival.model.summaries <- merge(mean.survival.model.summaries, variance.survival.model.summaries, by=by_names, all=T)

##-------------------------------
## Load and merge error metrics
##-------------------------------

# get rubin error
load(paste0(main_dir, "rubin_method_error.rdata"))
mean_rubin <- data.frame(mean.rubin.method.error)
setnames(mean_rubin, names(mean.rubin.method.error))
mean_rubin$data_transform <- rownames(mean_rubin)
mean_rubin <- data.table(melt(mean_rubin, id.vars="data_transform", variable.name="model_spec", value.name="mean_rubin_error"))
summary.survival.model.summaries <- merge(summary.survival.model.summaries, mean_rubin, by=c("data_transform", "model_spec"), all=T)

# get AIC
load(paste0(main_dir, "AIC.rdata"))
mean_aic <- data.frame(meanAIC)
setnames(mean_aic, names(meanAIC))
mean_aic$data_transform <- rownames(mean_aic)
mean_aic <- data.table(melt(mean_aic, id.vars="data_transform", variable.name="model_spec", value.name="mean_aic"))
summary.survival.model.summaries <- merge(summary.survival.model.summaries, mean_aic, by=c("data_transform", "model_spec"), all=T)


# get RMSE from OOS file
load(paste0(main_dir, "compiled_rmse.rdata"))
oos_rmse <- data.frame(compiled_rmse)
setnames(oos_rmse, colnames(compiled_rmse))
oos_rmse$data_transform <- rownames(oos_rmse)
oos_rmse <- data.table(melt(oos_rmse, id.vars="data_transform", variable.name="model_spec", value.name="oos_rmse"))
summary.survival.model.summaries <- merge(summary.survival.model.summaries, oos_rmse, by=c("data_transform", "model_spec"), all=T)


#get ranking as per RMSE
ranking <- unique(summary.survival.model.summaries[order(oos_rmse), list(data_transform, model_spec, oos_rmse)])
ranking[, ranking:= as.numeric(rownames(ranking))]
summary.survival.model.summaries <- merge(summary.survival.model.summaries, ranking, by=c("data_transform", "model_spec", "oos_rmse"), all=T)
summary.survival.model.summaries <- summary.survival.model.summaries[order(ranking, data_transform, model_spec)]

##-------------------------------
## Visualize
##-------------------------------

#format age types
summary.survival.model.summaries[, age_bins:= ifelse(model_spec %like% "c\\(", "binned", "continuous")]
summary.survival.model.summaries[age_bins=="binned", age_bins:= ifelse(model_spec %like% "20", "binned_10yr", "binned_15yr")]

#format covariate names
summary.survival.model.summaries[, covariate_combined:=gsub("spvl_(fraser|hybrid|model)", "spvl", covariate)]



## density plot of coefficient value across covariates, by model specification
pdf(paste0(main_dir, "summary_distributions.pdf"), width=14, height=8)

for (age_type in unique(summary.survival.model.summaries$age_bins)){
  
  subset <- summary.survival.model.summaries[age_bins==age_type]
  
  density <- ggplot(subset, aes(x=beta, group=model_spec)) +
    geom_density(aes(group=model_spec, color=model_spec, fill=model_spec), alpha=0.5, size=1) +
    facet_wrap(~covariate_combined, scales="free")
  
  print(density)

}


#density plot of aic/rmse/rubin across models
errors <- unique(summary.survival.model.summaries[, list(data_transform, model_spec, mean_rubin_error, mean_aic, oos_rmse)])
errors <- melt(errors, id.vars = c("data_transform", "model_spec"), variable.name = "error_type", value.name="error")

ggplot(errors, aes(x=error)) +
  geom_density(alpha=0.5, size=1) +
  facet_wrap(~error_type, scales="free")


graphics.off()





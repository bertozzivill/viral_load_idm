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

#run all models
print ("running survival models")

# fill a named list with model outputs
imputation_count <- unique(index.data.transform$imputation_count)
model_spec_names <- apply(index.survival.models,1,function(x)paste(x, collapse="-"))

survival.model.output<- NULL
list_idx <- 1

for (transform_name in colnames(data.for.survival)){
  print(transform_name)
  
  for (model_idx in 1:length(model_spec_names)){
    model_info <- index.survival.models[model_idx,]
    model_name <- model_spec_names[model_idx]
    print(model_name)

    for (imputation_name in rownames(data.for.survival)){
      print(imputation_name)
      this_data=data.for.survival[imputation_name, transform_name][[1]]
      
      output= LinearSurvivalModel(data=this_data,
                          spvl_method=model_info$spvl_method,
                          interaction=model_info$interaction,
                          bins=model_info$bins,
                          return.modelobject = 1)
      output_summary <- summary(output)$coefficients
      output_summary <- data.table(output_summary, 
                                   covariate=rownames(output_summary),
                                   imputation=as.numeric(gsub(".*=([0-9]*)", "\\1", imputation_name)),
                                   model_spec=model_name,
                                   data_transform=transform_name)
      survival.model.output[[list_idx]] <- output_summary
      list_idx<- list_idx+1
    }
  }
}

survival.model.output <- rbindlist(survival.model.output)

##-------------------------------
## Compute summary coefficients
##-------------------------------

#get summary coefficients for model outputs (across imputations)
by_names <- c("data_transform", "model_spec","covariate")

setnames(survival.model.output, c("Estimate", "Std. Error"), c("beta", "se"))

# means 
mean.survival.model.output <- survival.model.output[, list(beta=mean(beta)), by=by_names] 

#calculate standard errors using the method outlined on page 6 of the AMELIA documentation: https://cran.r-project.org/web/packages/Amelia/vignettes/amelia.pdf

#calculate variance of the mean estimates
survival.model.output[, across_beta_var:= var(beta), by=by_names]

#calculate mean of the variances
survival.model.output[, mean_variance:= mean(se^2), by=by_names]

#calculate joint variance of imputations
variance.survival.model.output<-unique(survival.model.output[, list(var=mean_variance + across_beta_var*(1+1/imputation_count)), by=by_names])
variance.survival.model.output[, se:=sqrt(var)]

#merge means and variances to get a full summary of results
summary.survival.model.output <- merge(mean.survival.model.output, variance.survival.model.output, by=by_names, all=T)

##-------------------------------
## Load and merge error metrics
##-------------------------------

# get rubin error
load(paste0(main_dir, "rubin_method_error.rdata"))
mean_rubin <- data.frame(mean.rubin.method.error)
setnames(mean_rubin, names(mean.rubin.method.error))
mean_rubin$data_transform <- rownames(mean_rubin)
mean_rubin <- data.table(melt(mean_rubin, id.vars="data_transform", variable.name="model_spec", value.name="mean_rubin_error"))
summary.survival.model.output <- merge(summary.survival.model.output, mean_rubin, by=c("data_transform", "model_spec"), all=T)

# get AIC
load(paste0(main_dir, "AIC.rdata"))
mean_aic <- data.frame(meanAIC)
setnames(mean_aic, names(meanAIC))
mean_aic$data_transform <- rownames(mean_aic)
mean_aic <- data.table(melt(mean_aic, id.vars="data_transform", variable.name="model_spec", value.name="mean_aic"))
summary.survival.model.output <- merge(summary.survival.model.output, mean_aic, by=c("data_transform", "model_spec"), all=T)


# get RMSE from OOS file
load(paste0(main_dir, "compiled_rmse.rdata"))
oos_rmse <- data.frame(compiled_rmse)
setnames(oos_rmse, colnames(compiled_rmse))
oos_rmse$data_transform <- rownames(oos_rmse)
oos_rmse <- data.table(melt(oos_rmse, id.vars="data_transform", variable.name="model_spec", value.name="oos_rmse"))
summary.survival.model.output <- merge(summary.survival.model.output, oos_rmse, by=c("data_transform", "model_spec"), all=T)

##-------------------------------
## Visualize
##-------------------------------

#format age types
summary.survival.model.output[, age_bins:= ifelse(model_spec %like% "c\\(", "binned", "continuous")]
summary.survival.model.output[age_bins=="binned", age_bins:= ifelse(model_spec %like% "20", "binned_10yr", "binned_15yr")]

#format covariate names
summary.survival.model.output[, covariate_combined:=gsub("spvl_(fraser|hybrid|model)", "spvl", covariate)]



## density plot of coefficient value across covariates, by model specification
pdf(paste0(main_dir, "summary_distributions.pdf"), width=14, height=8)

for (age_type in unique(summary.survival.model.output$age_bins)){
  
  subset <- summary.survival.model.output[age_bins==age_type]
  
  density <- ggplot(subset, aes(x=beta, group=model_spec)) +
    geom_density(aes(group=model_spec, color=model_spec, fill=model_spec), alpha=0.5, size=1) +
    facet_wrap(~covariate_combined, scales="free")
  
  print(density)

}


#density plot of aic/rmse/rubin across models
errors <- unique(summary.survival.model.output[, list(data_transform, model_spec, mean_rubin_error, mean_aic, oos_rmse)])
errors <- melt(errors, id.vars = c("data_transform", "model_spec"), variable.name = "error_type", value.name="error")

ggplot(errors, aes(x=error)) +
  geom_density(alpha=0.5, size=1) +
  facet_wrap(~error_type, scales="free")


graphics.off()





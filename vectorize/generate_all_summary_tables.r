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


#load imputed data and indices for data transforms and survial models
print("loading imputed data")
load(paste0(main_dir, "imputed_survival_data.rdata"))
load(paste0(main_dir, "indexes.rdata"))

#run all models
print ("running survival models")

# fill a nesting named list with model outputs: top layer is data transforms, middle layer is model specifications, 
# deepest layer is imputation count
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

#get summary coefficients for model outputs (across imputations)






#TODO: get rubin error (what data transformation, across imputations, has the lowest rubin?)
# (keeps crashing rstudio at the moment)
#rubin.method.error<-rubin.method(survival.model.output,imputation_count=imputation_count,output.type = 'error')
#mean.rubin.method.error<-lapply(rubin.method.error,function(x)rowMeans(x))

## get AIC for each model specification
# AIC<-output.AIC(survival.model.output,imputation_count=imputation_count)
# meanAIC<-lapply(AIC,rowMeans)


# TODO: get RMSE from OOS file (when my internet lives again)



#ok, putting this on hold until I can figure out why I can't summarize anything


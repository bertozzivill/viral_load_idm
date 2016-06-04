####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Survival model coupling nonlinear viral load progression with informative censoring
##              based events.
## Input:  main_dir: Directory in which the file "alldata.rdata" (either original or the training data)
##                    lives. If not running cross-validation, this will be one of:
##                    "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
##                    "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"
##                    
## Output: Saves a dataset called "survival.model.output" to main_dir, which contains a matrix of 
##          "lm" elements from which we can later predict and compare to the testing dataset if 
##            we're validating. Also generates a number of plots for whatever model is deemed best
##            by the AIC method of model selection. 
##
## Run from within the "vectorize" folder!
####################################################################################################

rm(list=ls())


#test.run=c(1,4)#run only selected rows of data.transform.index
test.run=0#run full data.transform.index

library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)

#automatically finds correct directory name so we can stop commenting and uncommenting lines
dir.exists <- function(d) {
  de <- file.info(d)$isdir
  ifelse(is.na(de), FALSE, de)
}
root_dir <- ifelse(dir.exists("/home/cselinger/"), "/home/cselinger/HIV-Cascade/merge/data/", "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/")
main_dir <- ifelse(length(commandArgs())>2, commandArgs()[3], root_dir)
validation <- ifelse(length(commandArgs())>2, T, F)

#load data
load(paste0(main_dir, "prepped_data.rdata"))
setnames(surv, "event_timeNew", "event_time_debiased")

############################################################################################################
## Run data transformations 
##############################################################################################################
imputation_count <- 10

index.data.transform<-expand.grid(upper_bound=c(2.9, 3.0, 3.1),
                                 debias=c(F,T),
                                 pre_1996_only=c(F,T),
                                 observed_only=c(F))
observed.only.index <- expand.grid(upper_bound=c(2.9), 
                                   debias=c(T,F),
                                   pre_1996_only=c(F,T),
                                   observed_only=c(T))
index.data.transform <- rbind(index.data.transform, observed.only.index)

#source("data.transform.R")

############################################################################################################
## Run Models 
##############################################################################################################

index.survival.models<-expand.grid(
  spvl_method=paste0('spvl_',c('model','fraser')),
  interaction_type=c("none", "two_way", "three_way"),
  include.age=T)
index.survival.models$spvl_method<-as.character(index.survival.models$spvl_method)
age.only <- list("none", "none", T)
spvl.only <- expand.grid(
  spvl_method=paste0('spvl_',c('model','fraser')),
  interaction_type="none",
  include.age=F)
null.model <- list("none", "none", F)
index.survival.models <- rbind(index.survival.models, age.only, spvl.only)

save(index.data.transform, index.survival.models, file=paste0(main_dir, "indexes.rdata"))

source("model.specification.R")

############################################################################################################
## Select best model (in-sample)
##############################################################################################################

if (!validation){
 source("model.selection.in_sample.R")
}



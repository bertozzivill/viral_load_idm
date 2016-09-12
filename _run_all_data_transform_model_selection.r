####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Survival model coupling nonlinear viral load progression with informative censoring
##              based events.
## Input:  this_main_dir: Directory in which the file "alldata.rdata" (either original or the training data)
##                    lives. If not running cross-validation, this will be one of:
##                    "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
##                    "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"
##                    
## Output: Saves a dataset called "survival.model.output" to this_main_dir, which contains a matrix of 
##          "lm" elements from which we can later predict and compare to the testing dataset if 
##            we're validating. Also generates a number of plots for whatever model is deemed best
##            by the AIC method of model selection. 
##
## Run from within the "vectorize" folder!
####################################################################################################

library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)

## is this a validation run? if yes, a variable called "validation" should exist, and it should be TRUE.
## otherwise, set "validation" to F
if (!(exists("validation"))) { validation <- F}

#set main directory
this_main_dir <- ifelse(validation, paste0(main_dir, iteration, "/", split, "/"), # if validation==T, variables called "iteration" and "split" should also exist
                        "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/")

#set age quintile cutoffs
age_quints <- c(15.4, 28.2, 41, 53.8, 66.6, Inf)
#age_quarts <- c(0, 30.4, 40.3, 50.2, Inf)

#load data
load(paste0(this_main_dir, "prepped_data.rdata"))
setnames(surv, "event_timeNew", "event_time_debiased")

############################################################################################################
## Run data transformations 
##############################################################################################################

## generate index of imputations we want to run

imputation_count <- 10

index.data.transform<-expand.grid(upper_bound=c(2.9, 3.0, 3.1),
                                 debias=c(F,T),
                                 pre_1996_only=c(F,T),
                                 observed_only=c(F),
                                 age_type=c("cont", "bin_10", "quint"))
observed.only.index <- expand.grid(upper_bound=c(2.9), 
                                   debias=c(F,T),
                                   pre_1996_only=c(F,T),
                                   observed_only=c(T),
                                   age_type=c("cont", "bin_10", "quint"))
index.data.transform <- rbind(index.data.transform, observed.only.index)

## run imputations based on inputs from index.data.transform
run_imputation <- T

if (run_imputation){ # we only ever want to impute on validation datasets now
  print("Running imputation")
  source("TransformData.R")
  
  run=apply(index.data.transform,1,function(x) paste(x,collapse='-'))
  
  data.for.survival<-mapply(TransformData,
                            upper_bound=index.data.transform$upper_bound,
                            debias=index.data.transform$debias,
                            pre_1996_only=index.data.transform$pre_1996_only,
                            observed_only=index.data.transform$observed_only,
                            age_type=index.data.transform$age_type,
                            MoreArgs=list(surv=surv)
  )
  
  rownames(data.for.survival)<-paste0("imputation_number=",c(1:imputation_count))
  colnames(data.for.survival)<-apply(index.data.transform,1,function(x)paste(x,collapse="-"))
  
  save(data.for.survival, index.data.transform, file=paste0(this_main_dir, "imputed_survival_data.rdata"))  
}else{
  load(file=paste0(this_main_dir,"imputed_survival_data.rdata"))
}

############################################################################################################
## Run Models 
##############################################################################################################

## generate index of models we want to run

index.survival.models<-expand.grid(
  spvl_method=paste0('spvl_',c('model','fraser')),
  interaction_type=c("none", "two_way", "three_way"),
  include.age=T)

index.survival.models$spvl_method<-as.character(index.survival.models$spvl_method)

age.only <- expand.grid(
  spvl_method="none",
  interaction_type="none",
  include.age=T,
  age.type= c("cont", "bin_10", "quint"))

spvl.only <- expand.grid(
  spvl_method=paste0('spvl_',c('model','fraser')),
  interaction_type="none",
  include.age=F,
  age.type= "none")

null.model <- list("none", "none", F, "none")

index.survival.models <- rbind(index.survival.models, age.only, spvl.only, null.model)

save(index.data.transform, index.survival.models, file=paste0(this_main_dir, "indexes.rdata"))

## run models based on inputs from index.survival.models

source("LinearSurvivalModel.R")

survival.model.output<-list()
print("running survival models")
for (k in 1:length(data.for.survival)){
  orig_data=data.table(data.for.survival[[k]])
  survival.model.output[[k]]<-mapply(LinearSurvivalModel,
                                     spvl_method=index.survival.models$spvl_method,
                                     interaction_type=index.survival.models$interaction_type,
                                     include.age=index.survival.models$include.age,
                                     age.type=index.survival.models$age.type,
                                     MoreArgs=list(orig_data=orig_data))
  colnames(survival.model.output[[k]]) <- apply(index.survival.models,1, function(x) paste(x, collapse="-"))
}

#generate names for this list
data_transform_names <- apply(index.data.transform,1, function(x) paste(x, collapse="-"))
data_transform_names <- data.table(expand.grid(1:10, data_transform_names))
data_transform_names <- data_transform_names[, list(Var2, imp=paste0("imp_count_",Var1))]
data_transform_names <- apply(data_transform_names, 1, function(x) paste(x, collapse="-"))
names(survival.model.output) <- data_transform_names

if (validation){
  print("calculating rmse")
 source("calculate_rmse.r")
}else{
  #save regression outputs for cross-validation, as well as the index values telling you what each list element means
  print("saving imputed survival data")
  save(survival.model.output, index.survival.models, index.data.transform, file=paste0(this_main_dir, "survival_model_output.rdata"))
}


############################################################################################################
## Select best model (in-sample)
##############################################################################################################

# if (!validation){
#  source("model.selection.in_sample.R")
# }



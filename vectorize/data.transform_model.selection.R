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


  ############################################################################################################
  ## loop over data transformations/imputations
  ##############################################################################################################
  
  
  index.data.transform=expand.grid(upper_bound=c(3.0, 3.2),
                                   debias=c(0,1),
                                   impute_type=c("with_vars", "no_vars"),
                                   impute.with.aids=c(F,T),
                                   imputation_count=10)
 
  run=apply(index.data.transform,1,function(x) paste(x,collapse='-'))

  if(test.run!=0){index.data.transform<-index.data.transform[c(test.run),]}
  
  source("data.transform.R")

  data.for.survival<-mapply(data.transform,
                            upper_bound=index.data.transform$upper_bound,
                            debias=index.data.transform$debias,
                            impute_type=index.data.transform$impute_type,
                            impute.with.aids=index.data.transform$impute.with.aids,
                            imputation_count=index.data.transform$imputation_count,
                            MoreArgs=list(surv=surv)
                            )

  rownames(data.for.survival)<-paste0("imputation_number=",c(1:index.data.transform$imputation_count[1]))
  colnames(data.for.survival)<-apply(index.data.transform,1,function(x)paste(x,collapse="-"))
  
  #names(out)<-paste0('imputed_data: ',"upper bound=", upper_bound," | debias=", debias," | impute_type=",impute_type,"| impute.with.aids=", impute.with.aids," || imputation_number=",1:imputation_count)
  save(data.for.survival, index.data.transform, file=paste0(main_dir, "imputed_survival_data.rdata"))  

  ############################################################################################################
  ## loop over survival models, per data transform and multiple imputations
  ##############################################################################################################
  
  index.survival.models<-expand.grid(
                              spvl_method=paste0('spvl_',c('model','fraser','hybrid')),
                              interaction_type=c("none", "two_way", "three_way"),
                              bins=list(0,c(seq(15,60,15),100),c(15,20,30,40,100)))
  
  index.survival.models$spvl_method<-as.character(index.survival.models$spvl_method)
  
  save(index.data.transform, index.survival.models, file=paste0(main_dir, "indexes.rdata"))
  
  load(file=paste0(main_dir,"imputed_survival_data.rdata"))
  source("LinearSurvivalModel.R")
  

  survival.model.output<-list()
  for (k in 1:length(data.for.survival)){
            data=data.table(data.for.survival[[k]])
            survival.model.output[[k]]<-mapply(LinearSurvivalModel,
            spvl_method=index.survival.models$spvl_method,
            interaction_type=index.survival.models$interaction_type,
            bins=index.survival.models$bins,
            MoreArgs=list(data=data))
  }
  
  #save regression outputs for cross-validation, as well as the index values telling you what each list element means
  save(survival.model.output, index.survival.models, index.data.transform, file=paste0(main_dir, "survival_model_output.rdata"))
  

##Rubins's method for multiple imputations
imputation_count=nrow(data.for.survival)

source("rubin.method.R")
rubin.method.error<-rubin.method(survival.model.output,imputation_count=imputation_count,output.type = 'error')
mean.rubin.method.error<-lapply(rubin.method.error,function(x)rowMeans(x))
save(rubin.method.error, mean.rubin.method.error, file=paste0(main_dir, "rubin_method_error.rdata"))

############################################################################################################
## modelselection per data transformation
##############################################################################################################

if (!validation){
  source("output.AIC.R")
  AIC<-output.AIC(survival.model.output,imputation_count=imputation_count)
  
  
  #get min AIC per data.transform (average over imputations)
  meanAIC<-lapply(AIC,rowMeans)
  A<-sapply(1:ncol(data.for.survival),function(k) {order(sapply(meanAIC,function(x){x[k]}))[1]})
  save(AIC, meanAIC, file=paste0(main_dir, "AIC.rdata"))
  
  #now: for data transform with min Rubin, take the model with min AIC
  
  delta.hats<-unlist(lapply(mean.rubin.method.error,function(x) order(x)[1]))
  mu.hat<-order(mapply(function(x,k){x[k]},meanAIC,unlist(lapply(mean.rubin.method.error,function(x) order(x)[1]))))[1]
  delta.hat<-delta.hats[mu.hat]
  
  smallest5.delta.hats<-unlist(lapply(mean.rubin.method.error[mu.hat],function(x) order(x)[1:5]))
  smallest5.mean.rubin.method.error<-unlist(lapply(mean.rubin.method.error[mu.hat],function(x) sort(x)[1:5]))
  
  
  ##results: best model with smallest imputation error average
  modelnames<-apply(index.survival.models,1,function(x) paste0(x,collapse="-"))
  
  print(paste0("BEST PERFORMER: DATA TRANSFORMATION ",
               colnames(data.for.survival)[delta.hat],
               " and MODEL ",
               modelnames[mu.hat],
               " with the 5 smallest RUBIN errors for the model",
               modelnames[mu.hat],
               colnames(data.for.survival)[smallest5.delta.hats],
               ": ",
               smallest5.mean.rubin.method.error
  ))
  
  #rerun bestmodel for ALL imputations, get combined coefficients/standard errors
  bestmodel<-lapply(data.for.survival[,delta.hat], function(this_data){
    output <- LinearSurvivalModel(
      data=this_data,
      return.modelobject=1,
      spvl_method=index.survival.models$spvl_method[mu.hat],
      interaction_type=index.survival.models$interaction_type[mu.hat],
      bins=unlist(index.survival.models$bins[mu.hat])
    )
    summary <- summary(output)$coefficients
    return(data.table(summary, covariate=rownames(summary)))
  }
  )
  
  bestmodel <- lapply(1:10, function(imp){
    bestmodel[[imp]][, imputation:=imp]
  })
  
  bestmodel <- do.call("rbind", bestmodel)
  setnames(bestmodel, c("Estimate", "Std. Error"), c("beta", "se"))
  
  #run function that calculates summary means and standard errors for a dataset of combined model results
  source("summarize_models.r")
  bestmodel_summary <- summarize_models(bestmodel)
  write.csv(bestmodel_summary, file=paste0(main_dir,"coefficients.csv", sep=''), row.names=F)
  
  bestmodel<-list('lm'=bestmodel_summary,'data'=data.for.survival,'name'=modelnames[mu.hat],'data.name'=colnames(data.for.survival)[delta.hat])
  save(bestmodel,file=paste0(main_dir,'bestmodel_in_sample.Rdata'))
  
  #do the same thing, but for the out-of-sample model: data transform 3.2-1-no_vars- TRUE-10 and model specification spvl_model-1-cont_age
  oos_transform <- "3.2-1-no_vars- TRUE-10"
  data_transform_index <- which(colnames(data.for.survival)==oos_transform)
  data_specification_index <- which(modelnames=="spvl_model-1-0")
  
  oos_bestmodel <-lapply(data.for.survival[,data_transform_index], function(this_data){
    output <- LinearSurvivalModel(
      data=this_data,
      return.modelobject=1,
      spvl_method=index.survival.models$spvl_method[data_specification_index],
      interaction_type=index.survival.models$interaction_type[data_specification_index],
      bins=unlist(index.survival.models$bins[data_specification_index])
    )
    summary <- summary(output)$coefficients
    return(data.table(summary, covariate=rownames(summary)))
  }
  )
  
  oos_bestmodel <- lapply(1:10, function(imp){
    oos_bestmodel[[imp]][, imputation:=imp]
  })
  
  oos_bestmodel <- do.call("rbind", oos_bestmodel)
  setnames(oos_bestmodel, c("Estimate", "Std. Error"), c("beta", "se"))
  
  #run function that calculates summary means and standard errors for a dataset of combined model results
  source("summarize_models.r")
  oos_bestmodel_summary <- summarize_models(oos_bestmodel)
  oos_bestmodel<-list('lm'=oos_bestmodel_summary,'data'=data.for.survival,'name'=modelnames[data_specification_index],'data.name'=colnames(data.for.survival)[data_transform_index])
  save(oos_bestmodel,file=paste0(main_dir,'bestmodel_out_of_sample.rdata'))
}



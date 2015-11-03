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
library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)

#main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
#main_dir <- "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"
  
main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/cross_validation/1/1/"

#load data
load(paste0(main_dir, "surv.rdata"))

  ############################################################################################################
  ## loop over data transformations/imputations
  ##############################################################################################################
  
  
  index.data.transform=expand.grid(upper_bound=3.7,
                                   debias=c(0,1),
                                   impute_type=c("with_vars", "no_vars"),
                                   impute.with.aids=c(F,T),
                                   imputation_count=3)

  run=apply(index.data.transform,1,function(x) paste(x,collapse='-'))

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
  save(data.for.survival, file=paste0(main_dir, "imputed_survival_data.rdata"))  

  ############################################################################################################
  ## loop over survival models, per data transform and multiple imputations
  ##############################################################################################################
  
  index.survival.models<-expand.grid(
                              spvl_method=paste0('spvl_',c('model','fraser','hybrid')),
                              interaction=c(0,1),
                              bins=list(c(0),c(seq(15,60,15),100)))
  
  index.survival.models$spvl_method<-as.character(index.survival.models$spvl_method)
  
  source("LinearSurvivalModel.R")
  
  survival.model.output<-list()
  for (k in 1:length(data.for.survival)){
            data=data.for.survival[[k]]
            survival.model.output[[k]]<-mapply(LinearSurvivalModel,
            spvl_method=index.survival.models$spvl_method,
            interaction=index.survival.models$interaction,
            bins=index.survival.models$bins)
  }
  
  #save regression outputs for cross-validation, as well as the index values telling you what each list element means
  save(survival.model.output, index.survival.models, index.data.transform, file=paste0(main_dir, "survival_model_output.rdata"))
  
  ##Rubins's method for multiple imputations
  source("rubin.method.R")
  rubin.method.error<-rubin.method(survival.model.output,imputation_count=3,output.type = 'error')
  mean.rubin.method.error<-lapply(rubin.method.error,function(x)rowMeans(x))
  
  ############################################################################################################
  ## modelselection per data transformation
  ##############################################################################################################
  
  source("output.AIC.R")
  AIC<-output.AIC(survival.model.output,imputation_count=3)
  minAIC<-lapply(AIC,function(x)order(x,decreasing=FALSE)[1])
  
  #get min AIC per data.transform (average over imputations)
  meanAIC<-lapply(AIC,rowMeans)
  A<-sapply(1:nrow(index.survival.models),function(k) {order(sapply(meanAIC,function(x){x[k]}))[1]})
  
  ##matrix, first row: model, second row data transform
  modelselection<-data.table(data.transform=1:length(run),bestmodel=A)
  modelselection[,min.Rubin.error:=
                   order(mapply(function(x,k){mean.rubin.method.error[[x]][[k]]},modelselection$bestmodel,modelselection$data.transform))[1]]
  
  ##results: best model with smallest imputation error average
  modelnames<-apply(index.survival.models,1,function(x) paste0(x,collapse="-"))
  
  print(paste0("BEST PERFORMER: DATA TRANSFORMATION ",
               run[modelselection$data.transform[1]],
               " and MODEL ",
               modelnames[modelselection$min.Rubin.error[1]]
  ))
  
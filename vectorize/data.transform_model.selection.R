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

#depending on who runs the code and where, you'll either need to look for Christian's VM, Amelia's home computer, or the cluster (the latter being a 
# passed argument for parallelized cross-validation)
dir.exists <- function(d) {
  de <- file.info(d)$isdir
  ifelse(is.na(de), FALSE, de)
}
root_dir <- ifelse(dir.exists("/home/cselinger/"), "/home/cselinger/HIV-Cascade/merge/data/", "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/")
main_dir <- ifelse(length(commandArgs()>2), commandArgs()[3], root_dir)

#load data
load(paste0(main_dir, "prepped_data.rdata"))

############################################################################################################
## loop over data transformations/imputations
##############################################################################################################

index.data.transform=expand.grid(upper_bound=3.7,
                                 debias=c(0,1),
                                 impute_type=c("with_vars", "no_vars"),
                                 impute.with.aids=c(F,T),
                                 imputation_count=10)

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
                            bins=list(0,c(seq(15,60,15),100),c(15,20,30,40,100)))

index.survival.models$spvl_method<-as.character(index.survival.models$spvl_method)

source("LinearSurvivalModel.R")
load(file=paste0(main_dir,"imputed_survival_data.rdata"))

survival.model.output<-list()
for (k in 1:length(data.for.survival)){
  survival.model.output[[k]]<-mapply(LinearSurvivalModel,
                                     spvl_method=index.survival.models$spvl_method,
                                     interaction=index.survival.models$interaction,
                                     bins=index.survival.models$bins,
                                     MoreArgs=list(data=data.for.survival[[k]]))
}

#save regression outputs for cross-validation, as well as the index values telling you what each list element means
save(survival.model.output, index.survival.models, index.data.transform, file=paste0(main_dir, "survival_model_output.rdata"))

##Rubins's method for multiple imputations
imputation_count=nrow(data.for.survival)

source("rubin.method.R")
rubin.method.error<-rubin.method(survival.model.output,imputation_count=imputation_count,output.type = 'error')
mean.rubin.method.error<-lapply(rubin.method.error,function(x)rowMeans(x))

############################################################################################################
## modelselection per data transformation
##############################################################################################################

source("output.AIC.R")
AIC<-output.AIC(survival.model.output,imputation_count=imputation_count)


#get min AIC per data.transform (average over imputations)
meanAIC<-lapply(AIC,rowMeans)
A<-sapply(1:ncol(data.for.survival),function(k) {order(sapply(meanAIC,function(x){x[k]}))[1]})

#now: for data transform with min Rubin, take the model with min AIC

delta.hats<-unlist(lapply(mean.rubin.method.error,function(x) order(x)[1]))
mu.hat<-order(mapply(function(x,k){x[k]},meanAIC,unlist(lapply(mean.rubin.method.error,function(x) order(x)[1]))))[1]
delta.hat<-delta.hats[mu.hat]

  
  ##results: best model with smallest imputation error average
  modelnames<-apply(index.survival.models,1,function(x) paste0(x,collapse="-"))

print(paste0("BEST PERFORMER: DATA TRANSFORMATION ",
             colnames(data.for.survival)[delta.hat],
             " and MODEL ",
             modelnames[mu.hat]
))


##rerun bestmodel
data=data.for.survival[,delta.hat][[1]]

bestmodel<-LinearSurvivalModel(return.modelobject=1,
                               spvl_method=index.survival.models$spvl_method[mu.hat],
                               interaction=index.survival.models$interaction[mu.hat],
                               bins=unlist(index.survival.models$bins[mu.hat]),
                               MoreArgs=list(data=data)
)


bestmodel<-list('lm'=bestmodel,'data'=data,'name'=modelnames[mu.hat],'data.name'=colnames(data.for.survival)[delta.hat])
save(bestmodel,file=paste0(main_dir,'bestmodel_in_sample.Rdata'))

source("lm2csv.R")
lm2csv(bestmodel$lm,paste0(main_dir,'table.bestmodel.coefficients.',confint=TRUE))



####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Select the best model, in-sample
## Input:  main_dir: Directory to which to save best model info
##         survival.model.outputs: list of model outputs
##                    
## Output: Saves a .txt called "bestmodel_in_sample" to main_dir, with information on the
##          in-sample best model
####################################################################################################

##Rubins's method for multiple imputations
imputation_count=nrow(data.for.survival)

source("rubin.method.R")
rubin.method.error<-rubin.method(survival.model.output,imputation_count=imputation_count,output.type = 'error')
mean.rubin.method.error<-lapply(rubin.method.error,function(x)rowMeans(x))
save(rubin.method.error, mean.rubin.method.error, file=paste0(main_dir, "rubin_method_error.rdata"))

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

# #rerun bestmodel for ALL imputations, get combined coefficients/standard errors
# bestmodel<-lapply(data.for.survival[,delta.hat], function(this_data){
#   output <- LinearSurvivalModel(
#     data=this_data,
#     return.modelobject=1,
#     spvl_method=index.survival.models$spvl_method[mu.hat],
#     interaction_type=index.survival.models$interaction_type[mu.hat],
#     bins=unlist(index.survival.models$bins[mu.hat])
#   )
#   summary <- summary(output)$coefficients
#   return(data.table(summary, covariate=rownames(summary)))
# }
# )
# 
# bestmodel <- lapply(1:10, function(imp){
#   bestmodel[[imp]][, imputation:=imp]
# })
# 
# bestmodel <- do.call("rbind", bestmodel)
# setnames(bestmodel, c("Estimate", "Std. Error"), c("beta", "se"))
# 
# #run function that calculates summary means and standard errors for a dataset of combined model results
# source("summarize_models.r")
# bestmodel_summary <- summarize_models(bestmodel)
# write.csv(bestmodel_summary, file=paste0(main_dir,"coefficients.csv", sep=''), row.names=F)
# 
# bestmodel<-list('lm'=bestmodel_summary,'data'=data.for.survival,'name'=modelnames[mu.hat],'data.name'=colnames(data.for.survival)[delta.hat])
# save(bestmodel,file=paste0(main_dir,'bestmodel_in_sample.Rdata'))

# #do the same thing, but for the out-of-sample model: data transform 3.2-1-no_vars- TRUE-10 and model specification spvl_model-1-cont_age
# oos_transform <- "3.2-1-no_vars- TRUE-10"
# data_transform_index <- which(colnames(data.for.survival)==oos_transform)
# data_specification_index <- which(modelnames=="spvl_model-1-0")
# 
# oos_bestmodel <-lapply(data.for.survival[,data_transform_index], function(this_data){
#   output <- LinearSurvivalModel(
#     data=this_data,
#     return.modelobject=1,
#     spvl_method=index.survival.models$spvl_method[data_specification_index],
#     interaction_type=index.survival.models$interaction_type[data_specification_index],
#     bins=unlist(index.survival.models$bins[data_specification_index])
#   )
#   summary <- summary(output)$coefficients
#   return(data.table(summary, covariate=rownames(summary)))
# }
# )
# 
# oos_bestmodel <- lapply(1:10, function(imp){
#   oos_bestmodel[[imp]][, imputation:=imp]
# })
# 
# oos_bestmodel <- do.call("rbind", oos_bestmodel)
# setnames(oos_bestmodel, c("Estimate", "Std. Error"), c("beta", "se"))
# 
# #run function that calculates summary means and standard errors for a dataset of combined model results
# source("summarize_models.r")
# oos_bestmodel_summary <- summarize_models(oos_bestmodel)
# oos_bestmodel<-list('lm'=oos_bestmodel_summary,'data'=data.for.survival,'name'=modelnames[data_specification_index],'data.name'=colnames(data.for.survival)[data_transform_index])
# save(oos_bestmodel,file=paste0(main_dir,'bestmodel_out_of_sample.rdata'))
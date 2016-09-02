####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Run the TransformData function on each of our index values of interest; save a 
##              matrix of imputed datasets
## Input:  main_dir: Directory to which to save imputed datasets
##         index.survival.models: data.frame with information about how to specify models
##                    
## Output: Saves a list of matrices called "survival.model.output" to main_dir, with model outpus
##          for each model specification in index.survival.models
##
####################################################################################################

load(file=paste0(main_dir,"imputed_survival_data.rdata"))
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

if (!validation){
  #save regression outputs for cross-validation, as well as the index values telling you what each list element means
  print("saving imputed survival data")
  save(survival.model.output, index.survival.models, index.data.transform, file=paste0(main_dir, "survival_model_output.rdata"))
}



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
  data=data.table(data.for.survival[[k]])
  survival.model.output[[k]]<-mapply(LinearSurvivalModel,
                                     spvl_method=index.survival.models$spvl_method,
                                     interaction_type=index.survival.models$interaction_type,
                                     bins=index.survival.models$bins,
                                     MoreArgs=list(data=data))
}

#save regression outputs for cross-validation, as well as the index values telling you what each list element means
print("saving imputed survival data")
save(survival.model.output, index.survival.models, index.data.transform, file=paste0(main_dir, "survival_model_output.rdata"))
####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Run the TransformData function on each of our index values of interest; save a 
##              matrix of imputed datasets
## Input:  main_dir: Directory to which to save imputed datasets
##         index.data.transform: data.frame with information about how to impute
##                    
## Output: Saves a matrix called "data.for.survival" to main_dir, with imputed entries for each 
##          data tranformation
##
####################################################################################################

source("TransformData.R")

run=apply(index.data.transform,1,function(x) paste(x,collapse='-'))

data.for.survival<-mapply(TransformData,
                          upper_bound=index.data.transform$upper_bound,
                          debias=index.data.transform$debias,
                          pre_1996_only=index.data.transform$pre_1996_only,
                          observed_only=index.data.transform$observed_only,
                          MoreArgs=list(surv=surv)
)

rownames(data.for.survival)<-paste0("imputation_number=",c(1:imputation_count))
colnames(data.for.survival)<-apply(index.data.transform,1,function(x)paste(x,collapse="-"))

#save(data.for.survival, index.data.transform, file=paste0(main_dir, "imputed_survival_data.rdata"))  



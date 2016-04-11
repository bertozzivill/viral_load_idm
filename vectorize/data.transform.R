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

if(test.run!=0){index.data.transform<-index.data.transform[c(test.run),]}

data.for.survival<-mapply(TransformData,
                          upper_bound=index.data.transform$upper_bound,
                          debias=index.data.transform$debias,
                          impute_type=index.data.transform$impute_type,
                          impute.with.aids=index.data.transform$impute.with.aids,
                          imputation_count=index.data.transform$imputation_count,
                          MoreArgs=list(surv=surv)
)

rownames(data.for.survival)<-paste0("imputation_number=",c(1:index.data.transform$imputation_count[1]))
colnames(data.for.survival)<-apply(index.data.transform,1,function(x)paste(x,collapse="-"))

save(data.for.survival, index.data.transform, file=paste0(main_dir, "imputed_survival_data.rdata"))  



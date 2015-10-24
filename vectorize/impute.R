####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Test different methods of imputation and their implications
## Output: TODO
####################################################################################################



#function to generate imputations based on the upper bounds and imputation count desired

impute <- function(upper_bound, imputation_count){
  #read in data
  print("importing data")
  load(paste0("vectorize/", "impute_data.rdata"))
  
    #no bounds on anything or bounded imputation
    if (upper_bound>10){
    imputed_data <- amelia(x=to_impute,cs='patient_id', m=imputation_count,p2s=0)
    }else{
      #bounded, with spvl
      col_of_interest <- which(colnames(to_impute)=="observed_survival")
      surv_len <- nrow(to_impute)
      
      surv_bds <- matrix(c(1:surv_len, #row
                           rep(col_of_interest, surv_len), #column
                           log(to_impute$event_time), #lower bound
                           rep(upper_bound, surv_len), #upper bound
                           rep(0.95, surv_len)
      ), 
      nrow=surv_len,
      ncol=5)
      imputed_data <- amelia(x=to_impute,cs='patient_id', priors=surv_bds, m=imputation_count,p2s=0)
    }



  #save imputed data
  save(imputed_data, file=paste0('vectorize/', "imputed_data.rdata"))
  
}


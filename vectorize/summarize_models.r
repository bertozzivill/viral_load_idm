summarize_models <- function(bestmodel){
  mean_bestmodel <- bestmodel[, list(beta=mean(beta)), by="covariate"] 
  
  #calculate standard errors using the method outlined on page 6 of the AMELIA documentation: https://cran.r-project.org/web/packages/Amelia/vignettes/amelia.pdf
  
  #calculate variance of the mean estimates
  bestmodel[, across_beta_var:= var(beta), by="covariate"]
  
  #calculate mean of the variances
  bestmodel[, mean_variance:= mean(se^2), by="covariate"]
  
  #calculate joint variance of imputations
  var_bestmodel<-unique(bestmodel[, list(var=mean_variance + across_beta_var*(1+1/imputation_count)), by="covariate"])
  var_bestmodel[, se:=sqrt(var)]
  
  #merge means and variances to get a full summary of results
  bestmodel_summary <- merge(mean_bestmodel, var_bestmodel, by="covariate", all=T)
  
  return(bestmodel_summary)
}

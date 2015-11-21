data.transform<-function(upper_bound=3.5,debias=1,impute_type='no_vars',impute.with.aids=T,imputation_count=2, surv){

  source('impute.R')
  print("imputing")
  
  if(debias==1) surv[,log_event_time:= log(event_timeNew)] else surv[,log_event_time:= log(event_time)]
  
  surv[, event_num:=ifelse(event_type=="mar", NA, ifelse(event_type=="aids", 1, 2))] #need to impute event type as well as time now
  surv[, observed_survival:= ifelse(is.na(event_num), NA, log_event_time)]  # equals observed time to AIDS/death if available, NA otherwise 
  surv[, Censorship:= ifelse(!is.na(event_num), "Uncensored", "Censored")]
  
  #if impute.with.aids is false, remove all art observations, and set all event nums to 2 (death)
  if (!impute.with.aids){
    surv[, event_num:= NULL] #event num will always be death
    aids_events <- surv[event_type=="aids"]
    surv <- surv[event_type!="aids"]
  }
  
  # define which columns you want to keep for imputation: this will change depending on whether you 
  # are running with or without vars (the formere will have more), and on the value of impute.with.aids
  # (if false, we can't include event_num since it's time invariant)
  if (impute.with.aids) event_num_col <- "event_num" else event_num_col <- NULL
  if (impute_type=="with_vars"){
    cols_to_keep <- c("patient_id","event_time","agesero","observed_survival","spvl_model","spvl_fraser","spvl_hybrid",event_num_col)
  } else{
    cols_to_keep <- c("patient_id","event_time","agesero","observed_survival",event_num_col)
  } 

  to_impute <- subset(surv,,cols_to_keep)
  
  ##cross-validation??
  #save(to_impute, file=paste0("main_dir","data_to_impute.rdata"))
  
  #run and save imputations
  imputed_data <- impute(to_impute, upper_bound, imputation_count)
  
  #save(imputed_data, file=paste0('main_dir', "imputed_data.rdata"))
  
  print(paste0("IMPUTATION: you have just saved imputed data", " | upper bound=", upper_bound," | debias=", debias," | impute_type=",impute_type,"| impute.with.aids=", impute.with.aids, " | imputation_count=",imputation_count))

  data<-imputed_data[['imputations']]
  
  if(impute_type=='no_vars'){
    data<- lapply(data,function(x) merge(x, surv[, list(patient_id, spvl_model, spvl_fraser, spvl_hybrid)], by="patient_id", all=T))
  }
  
  #restrict event_num to equal either one or two (also to be a factor)
  if (!impute.with.aids) data <- lapply(data, function(x) x[, event_num:=2])
  lapply(data,function(x) x[, event_num:=ifelse(event_num<1.5, 1, 2)])
  lapply(data,function(x) x[, event_num:= as.factor(event_num)])
  
  #if you've taken out the art observations: put them back
 # if (!impute.with.aids){
 #   final_cols <- c("patient_id","event_time","agesero","observed_survival","spvl_model","spvl_fraser","spvl_hybrid")
  #  aids_events <- subset(aids_events,,final_cols)
   # aids_events[, event_num:=1]
   # print(names(aids_events))
   # data <- lapply(data,function(x) rbind(x, aids_events))
 # } 
  
  out<-data
  return(out)
}

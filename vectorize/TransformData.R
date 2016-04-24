

TransformData<-function(upper_bound=3.2,debias=T,impute_type='with_vars',impute_with_aids=F, pre_1996_only=F, observed_only=F, imputation_count=10, surv){

  source('impute.R')
  print("imputing")
  
  # if running with pre-1996 data only, drop everything post-1996
  if (pre_1996_only) surv <- surv[pre_1996==T] 
  
  #Determine which event time (unbiased or not) to use; 
  if(debias) surv[,log_event_time:= log(event_timeNew)] else surv[,log_event_time:= log(event_time)]
  
  # set up event types and times that need to be imputed 
  surv[, event_num:=ifelse(event_type=="mar", NA, ifelse(event_type=="aids", 1, 2))] #need to impute event type as well as time now
  surv[, observed_survival:= ifelse(is.na(event_num), NA, log_event_time)]  # equals observed time to AIDS/death if available, NA otherwise 
  surv[, Censorship:= ifelse(!is.na(event_num), "Uncensored", "Censored")]
  
  # if running with observed events only, remove all mar events and return a list of the same format as the imputed datase
  if (observed_only){
    surv <- surv[event_type!="mar", list(patient_id, event_time, agesero ,observed_survival, event_num, spvl_model,spvl_fraser)]
    reps <- lapply(1:imputation_count, function(x) return(surv))
    names(reps) <- paste0("imp", 1:imputation_count)
    print("returning without imputation")
    return(reps)
  } 
  
  #if impute_with_aids is false, remove all art observations, and set all event nums to 2 (death)
  if (!impute_with_aids){
    surv[, event_num:= NULL] #event num will always be death
    aids_events <- surv[event_type=="aids"]
    surv <- surv[event_type!="aids"]
  }
  
  # define which columns you want to keep for imputation: this will change depending on whether you 
  # are running with or without vars (the formere will have more), and on the value of impute_with_aids
  # (if false, we can't include event_num since it's time invariant)
  if (impute_with_aids) event_num_col <- "event_num" else event_num_col <- NULL
  if (impute_type=="with_vars"){
    cols_to_keep <- c("patient_id","event_time","agesero","observed_survival","spvl_model","spvl_fraser",event_num_col)
  } else{
    cols_to_keep <- c("patient_id","event_time","agesero","observed_survival",event_num_col)
  } 

  to_impute <- subset(surv,,cols_to_keep)
  
  
  #run and return imputations
  imputed_data <- impute(to_impute, upper_bound, imputation_count)
  
  print(paste0("IMPUTATION: you have just run the imputation:", " | upper bound=", upper_bound," | debias=", debias," | pre_1996_only=",pre_1996_only))

  data<-imputed_data[['imputations']]
  
  if(impute_type=='no_vars'){
    data<- lapply(data,function(x) merge(x, surv[, list(patient_id, spvl_model, spvl_fraser)], by="patient_id", all=T))
  }
  
  #restrict event_num to equal either one or two (also to be a factor)
  if (!impute_with_aids) data <- lapply(data, function(x) x[, event_num:=2])
  lapply(data,function(x) x[, event_num:=ifelse(event_num<1.5, 1, 2)])
  lapply(data,function(x) x[, event_num:= factor(event_num, levels=c(1,2))])
  
  #if you've taken out the art observations: put them back
 if (!impute_with_aids){
   final_cols <- c("patient_id","event_time","agesero","observed_survival","spvl_model","spvl_fraser")
   aids_events <- subset(aids_events,,final_cols)
   aids_events[, event_num:=1]
   print(names(aids_events))
   data <- lapply(data,function(x) rbind(x, aids_events))
  } 
  
  out<-data
  return(out)
}

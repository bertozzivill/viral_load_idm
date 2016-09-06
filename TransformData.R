

TransformData<-function(upper_bound=3.2,debias=T, pre_1996_only=F, observed_only=F, imputation_count=10, surv){

  #print("imputing")
  
  # if running with pre-1996 data only, drop everything post-1996
  if (pre_1996_only) surv <- surv[pre_1996==T] 
  
  #Determine which event time (unbiased or not) to use; 
  if(debias) surv[,log_event_time:= log(event_time_debiased)] else surv[,log_event_time:= log(event_time)]
  
  # set up event types and times that need to be imputed 
  surv[, event_num:=ifelse(event_type=="mar", NA, ifelse(event_type=="death", 0, 1))] #need to impute event type as well as time now
  surv[, observed_survival:= ifelse(is.na(event_num), NA, log_event_time)]  # equals observed time to AIDS/death if available, NA otherwise 
  
  # if running with observed events only, remove all mar events and return a list of the same format as the imputed datase
  if (observed_only){
    surv <- surv[event_type!="mar", list(patient_id, event_time, agesero ,observed_survival, event_num, spvl_model,spvl_fraser)]
    surv[, event_num:=factor(event_num, levels=c(0,1))]
    reps <- lapply(1:imputation_count, function(x) return(surv))
    names(reps) <- paste0("imp", 1:imputation_count)
    print("returning without imputation")
    return(reps)
  } 
  
  #Remove all art observations, and set all event nums to 0 (death)
  surv[, event_num:= NULL] #event num will always be death
  aids_events <- surv[event_type=="aids"]
  surv <- surv[event_type!="aids"]

  # define which columns you want to keep for imputation:
  event_to_keep <- ifelse(debias, "event_time_debiased", "event_time")
  cols_to_keep <- c("patient_id",event_to_keep,"agesero","observed_survival","spvl_model","spvl_fraser")

  to_impute <- subset(surv,,cols_to_keep)
  aids_to_impute <- subset(aids_events,,cols_to_keep)
  
  #run imputations
  col_of_interest <- which(colnames(to_impute)=="observed_survival")
  surv_len <- nrow(to_impute)
  
  #set bounds
  surv_bds <- matrix(c(1:surv_len, #row
                       rep(col_of_interest, surv_len), #column
                       log(to_impute$event_time), #lower bound
                       rep(upper_bound, surv_len), #upper bound
                       rep(0.95, surv_len)
  ), 
  nrow=surv_len,
  ncol=5)
  imputed_data <- amelia(x=to_impute,cs='patient_id', priors=surv_bds, m=imputation_count,p2s=0)
  imputed_aids <- amelia(x=aids_to_impute,cs="patient_id", m=imputation_count,p2s=0)
  
  ## modify some columns, etc
  #print(paste0("IMPUTATION: you have just run the imputation:", " | upper bound=", upper_bound," | debias=", debias," | pre_1996_only=",pre_1996_only))

  data<-imputed_data[['imputations']]
  aids_data <- imputed_aids[['imputations']]
  
  #restrict event_num to equal zero (for imputed deaths) or one (aids deaths)
  data <- lapply(data, function(x) x[, event_num:=0])
  aids_data <- lapply(aids_data, function(x) x[, event_num:=1])
  
  #combine aids and non-aids data
   data <- lapply(1:imputation_count,function(x) rbind(data[[x]], aids_data[[x]]))
   lapply(data,function(x) x[, event_num:= factor(event_num, levels=c(0,1))])
  
  return(data)
}

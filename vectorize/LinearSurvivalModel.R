LinearSurvivalModel<-function(#data,
                              spvl_method='spvl_model',
                              interaction=0,
                              bins=c(seq(15,65,10),100),
                              return.modelobject=0){
  
  bins<-unlist(bins)
#  data<-data.table(data)
#   print(paste0("type: ",typeof(data)))
#   print(paste0("class: ",class(data)))
  ###agebins for categorical variable
  if(length(bins)>1){
    data[,agebin:=cut(data$agesero, bins, include.lowest=TRUE)]
  }
  
  ###arguments of function
  
  binning<-ifelse(length(bins)>1,
         paste0(" and age bins ", paste(bins,collapse=' ')),
         ' and continous age covariate')
         
  print(paste0('SURVIVAL: ',"log normal survival with ",spvl_method,c(' without interaction',' with interaction')[interaction+1],binning))
  
  if (interaction==1&length(bins)>1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  (agebin + ",spvl_method,")^2 + event_num")}
  if (interaction==1&length(bins)>1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  (agebin + ",spvl_method,")^2")}
  
  if (interaction==1&length(bins)==1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  (agesero + ",spvl_method,")^2 + event_num")}
  if (interaction==1&length(bins)==1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  (agesero + ",spvl_method,")^2")}
  
  if (interaction==0&length(bins)>1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  agebin + ",spvl_method," + event_num")}
  if (interaction==0&length(bins)>1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  agebin + ",spvl_method)}
  
  if (interaction==0&length(bins)==1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  agesero + ",spvl_method," + event_num")}
  if (interaction==0&length(bins)==1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  agesero + ",spvl_method)}

  ###evaluate model
  output <- lm(as.formula(formula), data=data)
  if(return.modelobject==0){
    return(list('lm'=output,'AIC'=round(AIC(output))))
  }else{
    return(output)
  }

}


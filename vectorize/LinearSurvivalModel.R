LinearSurvivalModel<-function(data,
                              spvl_method='spvl_model',
                              interaction_type="two_way",
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
         
  print(paste0("SURVIVAL: ","log normal survival with ",spvl_method, ", interaction type ", interaction_type,binning))
  
  age_str <- ifelse(length(bins)>1, "agebin", "agesero")
  
  #specify formula
  if (interaction_type=="none"){
    formula = paste("observed_survival ~ ", age_str, "+", spvl_method, "+ event_num")
  }else if (interaction_type=="two_way"){
    formula = paste("observed_survival ~ (", age_str, "+", spvl_method, ")^2 + event_num")
  }else if (interaction_type=="three_way"){
    formula = paste("observed_survival ~ (", age_str, "+", spvl_method, "+ event_num)^3")
  }
  
  print(formula)

  ###evaluate model
  output <- lm(as.formula(formula), data=data)
  if(return.modelobject==0){
    return(list('lm'=output,'AIC'=round(AIC(output))))
  }else{
    return(output)
  }

}


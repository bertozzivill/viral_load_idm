LinearSurvivalModel<-function(orig_data,
                              spvl_method='spvl_model',
                              interaction_type="none",
                              include.age=T,
                              age.type="cont",
                              return.modelobject=T){
  
  
  ###arguments of function
  
  
  ## immediately make a copy of the original dataset. 
  ## Use this instead of orig_data so we can manipulate it without interfereing with the mapply.
  data <- copy(orig_data)
         
  #print(paste0("SURVIVAL: ","log normal survival with ",spvl_method, ", interaction type ", interaction_type, " and ", ifelse(include.age, "age", "no age")))
  
  # if you have only spvl or age (but not both), make sure interaction_type is none
  if (interaction_type!="none" & (spvl_method=="none" | include.age==F)){
    stop(paste("interaction type is", interaction_type, ", but spvl method is", spvl_method, "and age is", ifelse(include.age, "included", "not included")))
  }
  
  age_str <- ifelse(include.age, "agesero +", "")
  spvl_str <- ifelse(spvl_method=="none", "", paste0(spvl_method, ""))
  
  #specify formula
  if (interaction_type=="none"){
    model_formula <- paste("observed_survival ~ ", age_str, spvl_str, "+ event_num")
  }else if (interaction_type=="two_way"){
    model_formula <- paste("observed_survival ~ (agesero +", spvl_method, ")^2 + event_num")
  }else if (interaction_type=="three_way"){
    model_formula <- paste("observed_survival ~ (agesero +", spvl_method, "+ event_num)^3")
  }
  
  #print(model_formula)

  ###evaluate model
  output <- lm(as.formula(model_formula), data=data)
  
  if(return.modelobject){
    return(list('lm'=output,'AIC'=round(AIC(output))))
  }else{
    return(output)
  }
  
  
}


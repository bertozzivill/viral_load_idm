rubin.method<-function(survival.model.output,imputation_count,output.type=c('error')){
  
  m=imputation_count
  k=length(survival.model.output)/m-1
  output<-rep(list(NULL),nrow(index.survival.models));names(output)<-apply(index.survival.models,1,function(x) paste(x,sep=".",collapse="."))
  
  print(paste0("Rubins's method for ",imputation_count," imputations...using fixed effect estimates and se...returning ",output.type, " of t-test"))
  for (model in 1:length(output)){
 
  out_p.value=out_error<-NULL
  
  for (j in 0:k){
  coeff=sapply(c(1:m)+j*m,function(x) summary(survival.model.output[[x]]["lm",model][[1]])$coefficients[,"Estimate"])
  se=sapply(c(1:m)+j*m,function(x) summary(survival.model.output[[x]]["lm",model][[1]])$coefficients[,"Std. Error"])
  
  Qbar<-apply(coeff,1,mean)
  Ubar<-apply(se,1,mean)
  B<-apply(coeff,1,var)
  TotalVar<-Ubar+(1+1/imputation_count)*B
  df<-(m-1)*(1+m*Ubar/((m+1)*B))^2
  
  error <- qt(0.975,df=df)*sqrt(TotalVar)/sqrt(imputation_count)
  p.value = round(2*pt(-abs(Qbar/sqrt(TotalVar)), df=df),3)
  out_error<-rbind(out_error,error)
  out_p.value<-rbind(out_p.value,p.value)
  }
  rownames(out_error)=rownames(out_p.value)<-run
  if (output.type=='error'){
  output[[model]]<-out_error}else
  {
  output[[model]]<-out_p.value
  }
  
  }
  
  return(output)
}


  

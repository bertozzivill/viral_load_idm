output.AIC<-function(survival.model.output,imputation_count){
  
  id=imputation_count
  k=length(survival.model.output)/id-1
  output<-rep(list(NULL),nrow(index.survival.models));names(output)<-apply(index.survival.models,1,function(x) paste(x,sep=".",collapse="."))
  
  for (model in 1:length(output))
  {
    
    out_AIC<-NULL
    
    for (j in 0:k){
      AIC=sapply(c(1:id)+j*id,function(x) survival.model.output[[x]]["AIC",model][[1]])
      out_AIC<-rbind(out_AIC,AIC)
    }
    rownames(out_AIC)<-run
    output[[model]]<-out_AIC
    
  }
  
  return(output)
}




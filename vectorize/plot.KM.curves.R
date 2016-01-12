




plot.KM.curves<-function(surv, show_censored=T){
  
  
  event_names <- data.table(event_type=c("aids", "death", "mar"), 
                            event_name=c("AIDS Only", "Death Only", "Censored"))
  
  surv <- merge(surv, event_names, by="event_type")
  
  library(survival)
  ###agebins for categorical variable
  bins<-c(seq(15,60,15),100)
  surv=data.table(surv, agebin=cut(surv$agesero, bins, include.lowest=FALSE))
  surv[, censor_indic:=ifelse(event_type=="mar", 0,1)] #binary indicator for survival model
  
  if (show_censored==F) surv <- surv[event_type!="mar"]
  
  #everybody all at once
  km <- survfit(Surv(surv$event_time, surv$censor_indic) ~ surv$agebin)
  all <- ggsurv(km, cens.col="black", main="All Events")+ylim(c(0,1))+xlim(c(0,15))
  
  km <- survfit(Surv(surv$event_timeNew, surv$censor_indic) ~ surv$agebin)
  allNew <- ggsurv(km, cens.col="black", main="All Events, Unbiased Event Time")+ylim(c(0,1))+xlim(c(0,15))
  
  #only AIDS
  subset<- surv[event_type=="aids"]
  km <- survfit(Surv(subset$event_time, subset$censor_indic) ~ subset$agebin)
  aids <- ggsurv(km, cens.col="black", main="AIDS Events")
  
  #only death (with/without aids)
  subset<- surv[event_type=="death"]
  km <- survfit(Surv(subset$event_time, subset$censor_indic) ~ subset$agebin)
  death <- ggsurv(km, cens.col="black", main="Death Events")
  
  multiplot(all,allNew, aids, death,cols=2)
   
}
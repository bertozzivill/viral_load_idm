




plot.KM.curves<-function(surv){
  
  
  event_names <- data.table(event_type=c("aids", "death", "aids_death", "mar"), 
                            event_name=c("AIDS Only", "Death Only", "AIDS and Death", "Censored"))
  
  surv <- merge(surv, event_names, by="event_type")
  
  library(survival)
  ###agebins for categorical variable
  bins<-c(seq(15,60,15),100)
  surv=data.table(surv, agebin=cut(surv$agesero, bins, include.lowest=FALSE))
  surv[, censor_indic:=ifelse(event_type=="mar", 0,1)] #binary indicator for survival model
  
  
  #everybody all at once
  km <- survfit(Surv(surv$event_time, surv$censor_indic) ~ surv$agebin)
  all <- ggsurv(km, cens.col="black", main="KM survival estimator, All Event Types \n no VL matching")+ylim(c(0,1))+xlim(c(0,15))
  
  km <- survfit(Surv(surv$event_timeNew, surv$censor_indic) ~ surv$agebin)
  allNew <- ggsurv(km, cens.col="black", main="KM survival estimator, All Event Types, bias removed \n no VL matching \n pre 1996 events")+ylim(c(0,1))+xlim(c(0,15))
  
  #only AIDS
  subset<- surv[event_type=="aids"]
  km <- survfit(Surv(subset$event_time, subset$censor_indic) ~ subset$agebin)
  aids <- ggsurv(km, cens.col="black", main="Survival, AIDS Only")
  
  #only death (with/without aids)
  subset<- surv[event_type=="death" | event_type=="aids_death"]
  km <- survfit(Surv(subset$event_time, subset$censor_indic) ~ subset$agebin)
  any_death <- ggsurv(km, cens.col="black", main="Survival, Any Death")
  
  subset<- surv[event_type=="aids" | event_type=="mar" ]
  km <- survfit(Surv(subset$event_timeNew, subset$censor_indic) ~ subset$agebin)
  death_mar <- ggsurv(km, cens.col="black", main="Survival, Any Death, AIDS, mar")
  
  

  multiplot(aids,any_death,death_mar,cols=2)
  
  
  
}
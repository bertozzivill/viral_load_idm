#################################################################################################
## prep_competing_risks_cascade.r
## Description: Prep the SPVL and agesero covariates for the viral load models, save a prepped 
##              dataset.
#################################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)

main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
#main_dir <- "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"
#main_dir <- "/home/cselinger/HIV-Cascade/merge/data/"


##################################################
## I. Define a function that runs the model, 
##    calculates random effects and spvl
##################################################

run_nonlin <- function(model_formula, data){
  
  Model<- ~ b0+ b1*time + b2*exp(-b3*time)
  ModelGradient<-deriv(Model,namevec=c("b0", "b1", "b2","b3"),function.arg=c("time","b0", "b1", "b2","b3"))
  
  print(paste("Hello I am the model formula and my name is:", model_formula))
  out<-do.call("nlmer", list( as.formula(model_formula),
                              data=quote(data),
                              start = c(b0=3,b1=0.05,b2=1,b3=4),
                              control=nlmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=200000))))
  print(out)
  
  re<- ranef(out)[[1]]
  re$patient_id<- rownames(re)
  re <- data.table(re)
  
  fe <- fixef(out)
  
  #parameters for each study subject
  b=cbind(fe['b0']+re$b0,fe['b1']+re$b1,fe['b2'],fe['b3']);colnames(b)=paste('b',c(0:3),sep="");rownames(b)=re$patient_id;b=as.data.frame(b)
  
  #setpoint viral load defined as global minimum of the viral load function, see also the pdf I had sent before
  
  timetosetpoint<-function(b){ -(log(b[1]/(b[2]*b[3])))/b[3] };
  b$tts<-as.vector(apply(b[,c('b1','b2','b3')],1,timetosetpoint))
  global_tts <- timetosetpoint(c(fe[['b1']], fe[['b2']], fe[['b3']]))
  nl_vl_model<-function(x){x[1]+x[2]*x[5]+x[3]*exp(-x[4]*x[5])}
  b$spvl<-as.vector(apply(b[,c('b0','b1','b2','b3','tts')],1,nl_vl_model))
  
  return(list(out, b, global_tts))
}

#bonus: function for calculating geometric means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x>0]))
}

################################
## II. Read in data
#################################
print("importing data")
load(paste0(main_dir, "alldata.rdata"))

#generate viral load dataset
vl <- alldata[,list(patient_id, time=visit_time, vl, assay_ll, assay_type)]

#generate survival dataset
surv <- alldata[, list(patient_id, event_type, event_time, event_timeNew, agesero=serocon_age,event_date,serocon_date,enroll_date,seroconv_before_enroll,aids_before_enroll,bias)]
setkeyv(surv, NULL)
surv <- unique(surv)

orig_surv <- copy(surv) #keep a copy with original AD and ADtime values for plotting
rm(alldata)

#################################################################
## III. Format Viral Load Data:
##   -- add age at seroconversion
##   -- rename columns to show marker (_M) loyalty
#################################################################
print("formatting viral load data")
vl[, vl_unlog:=vl]
vl[, vl:=log10(vl+2)] #convert viral load to log10 space
orig_vl <- copy(vl)
age <- orig_surv[, list(patient_id, agesero)]   #get a dataset with age at seroconversion to merge onto vl
vl <- merge(vl, age, by="patient_id", all=T)

## compute fraser viral load
print("computing fraser viral load")
fraser_spvl <- vl[time>0.5, list(patient_id, vl)]
fraser_spvl[, spvl_fraser:=gm_mean(vl), by="patient_id"]
fraser_spvl$vl <- NULL
fraser_spvl <- fraser_spvl[!duplicated(fraser_spvl)]

#############################################
## IV. Run nonlinear vl model
#############################################
print("running viral load only model")
#source("vectorize/lm2csv.R")
formula_vl <- "vl~ModelGradient(time=time,b0,b1,b2,b3)~ (b0|patient_id) + (b1|patient_id)"
output_vl <- run_nonlin(formula_vl, data=vl)
re_vl <- output_vl[[2]]
lm.spvl_model <- output_vl[[1]]
missing <- rownames(re_vl[is.na(re_vl$spvl),])

save(missing,re_vl,vl,file=paste0(main_dir,"missing_spvl_model.Rdata"))
#lm2csv(lm.spvl_model,"../data/table.spvl_model.",confint=FALSE)

new_vl <- vl[!patient_id %in% missing]
print("rerunning viral load model")
formula_vl <- "vl~ModelGradient(time=time,b0,b1,b2,b3)~ (b0|patient_id) + (b1|patient_id)"
output_vl <- run_nonlin(formula_vl, data=new_vl)
re_vl <- output_vl[[2]]
re_vl$patient_id <- rownames(re_vl)
re_vl<-data.table(re_vl)[, list(spvl, patient_id)]

#estimate spvl be geometric mean from global tts
print("computing hybrid viral load")
global_tts <- output_vl[[3]]
hybrid_spvl <- vl[time>global_tts, list(patient_id, vl)]
hybrid_spvl[, spvl_hybrid:=gm_mean(vl), by="patient_id"]
hybrid_spvl$vl <- NULL
hybrid_spvl <- hybrid_spvl[!duplicated(hybrid_spvl)]

surv <- merge(surv, re_vl[, list(patient_id, spvl_model=spvl)], by="patient_id", all=T)
surv <- merge(surv, fraser_spvl, by="patient_id", all=T)
surv <- merge(surv, hybrid_spvl, by="patient_id", all=T)

save(surv, file=paste0(main_dir, "prepped_data.rdata"))


####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Survival model coupling nonlinear viral load progression with informative censoring
##              based events.
## Output: TODO
####################################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)

#setwd("C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/code/")
setwd("C:/Users/cselinger/Dropbox (IDM)/viral_loadNew/cascade/code")

parent_dir <- "../data/"



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
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

################################
## II. Read in data
#################################
print("importing viral load data")
vl <- fread(paste0(parent_dir, "vl.csv")) #column 'time' corresponds to 'visit_time' in  alldata.rdata 
surv <- fread(paste0(parent_dir, "surv.csv")) #column 'AD' corresponds to 'aids_death_indic', 'ADtime' to 'event_time' in alldata.rdata
orig_surv <- copy(surv) #keep a copy with original AD and ADtime values for plotting

# ggplot(surv, aes(x=event_timeNew)) +
#   geom_histogram() +
#   facet_wrap(~event_type, scales="free_y")

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
formula_vl <- "vl~ModelGradient(time=time,b0,b1,b2,b3)~ (b0|patient_id) + (b1|patient_id)"
output_vl <- run_nonlin(formula_vl, data=vl)
re_vl <- output_vl[[2]]

missing <- rownames(re_vl[is.na(re_vl$spvl),])

new_vl <- vl[!patient_id %in% missing]
print("rerunnig viral load model")
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

############################################################################################################
## loop over data transformations/imputations
##############################################################################################################


index.data.transform=expand.grid(upper_bound=c(seq(3.2,4,0.5),999),
                                 debias=c(0,1),
                                 impute_type=c("with_vars", "no_vars"),
                                 impute.with.aids=c(F,T),
                                 imputation_count=3)

run=apply(index.data.transform,1,function(x) paste(x,collapse='-'))

source("vectorize/data.transform.R")


data.for.survival<-mapply(data.transform,
                          upper_bound=index.data.transform$upper_bound,
                          debias=index.data.transform$debias,
                          impute_type=index.data.transform$impute_type,
                          impute.with.aids=index.data.transform$impute.with.aids,
                          imputation_count=index.data.transform$imputation_count
                          )
rownames(data.for.survival)<-paste0("imputation_number=",c(1:index.data.transform$imputation_count[1]))
colnames(data.for.survival)<-apply(index.data.transform,1,function(x)paste(x,collapse="-"))

#names(out)<-paste0('imputed_data: ',"upper bound=", upper_bound," | debias=", debias," | impute_type=",impute_type,"| impute.with.aids=", impute.with.aids," || imputation_number=",1:imputation_count)



############################################################################################################
## loop over survival models, per data transform and multiple imputations
##############################################################################################################


index.survival.models<-expand.grid(
                            spvl_method=paste0('spvl_',c('model','fraser','hybrid')),
                            interaction=c(0,1),
                            bins=list(c(0),c(seq(15,60,15),100)))

index.survival.models$spvl_method<-as.character(index.survival.models$spvl_method)


source("vectorize/LinearSurvivalModel.R")

survival.model.output<-list()
for (k in 1:length(data.for.survival)){
          data=data.for.survival[[k]]
          survival.model.output[[k]]<-mapply(LinearSurvivalModel,
          spvl_method=index.survival.models$spvl_method,
          interaction=index.survival.models$interaction,
          bins=index.survival.models$bins)
}

##Rubins's method for multiple imputations
source("vectorize/rubin.method.R")
rubin.method.error<-rubin.method(survival.model.output,imputation_count=3,output.type = 'error')
mean.rubin.method.error<-lapply(rubin.method.error,function(x)rowMeans(x))

############################################################################################################
## modelselection per data transformation
##############################################################################################################

source("vectorize/output.AIC.R")
AIC<-output.AIC(survival.model.output,imputation_count=3)
minAIC<-lapply(AIC,function(x)order(x,decreasing=FALSE)[1])


#get min AIC per data.transform (average over imputations)
meanAIC<-lapply(AIC,rowMeans)
A<-sapply(1:nrow(index.survival.models),function(k) {order(sapply(meanAIC,function(x){x[k]}))[1]})

##matrix, first row: model, second row data transform
modelselection<-data.table(data.transform=1:length(run),bestmodel=A)
modelselection[,min.Rubin.error:=
                 order(mapply(function(x,k){mean.rubin.method.error[[x]][[k]]},modelselection$bestmodel,modelselection$data.transform))[1]]


##results: best model with smallest imputation error average
modelnames<-apply(index.survival.models,1,function(x) paste0(x,collapse="-"))

print(paste0("BEST PERFORMER: DATA TRANSFORMATION ",
             run[modelselection$data.transform[1]],
             " and MODEL ",
             modelnames[modelselection$min.Rubin.error[1]]
))



# relativeAIC<-lapply(survival.model.output,function(x) exp((min(unlist(x[2,]))-unlist(x[2,]))/2))

#########################
####PLOTS

data=data.for.survival[modelselection$bestmodel.min.Rubin.error[1]][[1]]
bestmodel<-LinearSurvivalModel(return.modelobject=1,
                    spvl_method=index.survival.models$spvl_method[modelselection$min.Rubin.error[1]],
                    interaction=index.survival.models$interaction[modelselection$min.Rubin.error[1]],
                    bins=index.survival.models$bins[modelselection$min.Rubin.error[1]])


###FIGURE 3b
source("vectorize/plot.survival.R")

plot.survival(bestmodel,modelnames[modelselection$min.Rubin.error[1]],'title')

###FIGURE 3c
source("vectorize/plot.modelcurve.R")
plot.modelcurve(bestmodel,modelnames[modelselection$min.Rubin.error[1]],'title')



source('c:/Users/cselinger/TOOLS/multiplot.R')
plots<-lapply(1:length(bestmodel.names),function(x) plot.modelcurve(bestmodel[[x]],bestmodel.names[[x]],names(bestmodel.index)[x]))

png('vectorize/modelcurve.png',width=45, height=45, units='cm', res=400)
multiplot(plotlist=plots,cols=4)
dev.off()


###FIGURE 2a
source("vectorize/plot.spvlagesero.R")
plot.spvlagesero(bestmodel,modelnames[modelselection$min.Rubin.error[1]],'title')


plots<-lapply(1:length(bestmodel.names),function(x) plot.spvlagesero(bestmodel[[x]],bestmodel.names[[x]],'Set point viral load and age at seroconversion'))

png('vectorize/spvlagesero.png',width=45, height=45, units='cm', res=400)
multiplot(plotlist=plots,cols=4)
dev.off()


###FIGURE 2b
prop.na.spvl_model<-as.data.frame(table(substring(data_orig[[16]][is.na(spvl_model)]$patient_id, 1, 3))/table(substring(data_orig[[16]]$patient_id, 1, 3)))

ggplot(prop.na.spvl_model)+
  geom_bar(aes(x=Var1,y=Freq,group=Var1,fill=Var1),stat='identity')+
  theme(legend.position = "none")+
  labs(x="cohort", y="fraction",title="Negative set point viral load slope")


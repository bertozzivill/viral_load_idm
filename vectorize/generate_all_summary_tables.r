####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Find summary measures (coefficients, SE, AIC, Rubin error, RMSE[oos]) for all models tested

## Input:  main_dir: Directory in which imputed model data live
##                    
## Output: 
##
## Run from within the "vectorize" folder!
####################################################################################################


rm(list=ls())

library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)
library(stringr)
library(Rmisc)

#automatically finds correct directory name so we can stop commenting and uncommenting lines
dir.exists <- function(d) {
  de <- file.info(d)$isdir
  ifelse(is.na(de), FALSE, de)
}
root_dir <- ifelse(dir.exists("/home/cselinger/"), "/home/cselinger/HIV-Cascade/merge/data/", "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/")
main_dir <- ifelse(length(commandArgs())>2, commandArgs()[3], root_dir)


##---------------------------
## Load data
##---------------------------

#load imputed data and indices for data transforms and survial models
print("loading imputed data")
load(paste0(main_dir, "imputed_survival_data.rdata"))
load(paste0(main_dir, "indexes.rdata"))
load(paste0(main_dir, "survival_model_output.rdata"))
load(paste0(main_dir, "prepped_data.rdata"))

# fill a named list with model outputs
imputation_count <- 10
model_spec_names <- apply(index.survival.models,1,function(x)paste(x, collapse="-"))

compute_summary <- T

if (compute_summary){
  survival.model.summaries<- NULL
  list_idx <- 1
  
  for (transform_imp in names(survival.model.output)){
    print(transform_imp)
    
    #find data transform and imputation county
    pattern <- "(.*)-imp_count_([0-9]*)"
    pattern_match <- str_match(transform_imp, pattern)
    data_transform <- pattern_match[1,2]
    imp_number <- as.numeric(pattern_match[1,3])
    
    this_transform<-survival.model.output[[transform_imp]]
    
    for (model_spec in colnames(this_transform)){
      this_model <- this_transform[, model_spec]
      model_summary <- summary(this_model$lm)$coefficients
      
      model_summary <- data.table(model_summary, 
                                  covariate=rownames(model_summary),
                                  imputation=imp_number,
                                  model_spec=model_spec,
                                  data_transform=data_transform)
      survival.model.summaries[[list_idx]] <- model_summary
      list_idx<- list_idx+1
      
    }
    
  }
  
  survival.model.summaries <- rbindlist(survival.model.summaries)



  ##-------------------------------
  ## Compute summary coefficients
  ##-------------------------------
  
  #get summary coefficients for model outputs (across imputations)
  by_names <- c("data_transform", "model_spec","covariate")
  
  setnames(survival.model.summaries, c("Estimate", "Std. Error"), c("beta", "se"))
  
  # means 
  mean.survival.model.summaries <- survival.model.summaries[, list(beta=mean(beta)), by=by_names] 
  
  #calculate standard errors using the method outlined on page 6 of the AMELIA documentation: https://cran.r-project.org/web/packages/Amelia/vignettes/amelia.pdf
  
  #calculate variance of the mean estimates
  survival.model.summaries[, across_beta_var:= var(beta), by=by_names]
  
  #calculate mean of the variances
  survival.model.summaries[, mean_variance:= mean(se^2), by=by_names]
  
  #calculate joint variance of imputations
  variance.survival.model.summaries<-unique(survival.model.summaries[, list(var=mean_variance + across_beta_var*(1+1/imputation_count)), by=by_names])
  variance.survival.model.summaries[, se:=sqrt(var)]
  
  #merge means and variances to get a full summary of results
  summary.survival.model.summaries <- merge(mean.survival.model.summaries, variance.survival.model.summaries, by=by_names, all=T)
  
  ##-------------------------------
  ## Load and merge error metrics
  ##-------------------------------
  
  # # get rubin error
  # load(paste0(main_dir, "rubin_method_error.rdata"))
  # mean_rubin <- data.frame(mean.rubin.method.error)
  # setnames(mean_rubin, names(mean.rubin.method.error))
  # mean_rubin$data_transform <- rownames(mean_rubin)
  # mean_rubin <- data.table(melt(mean_rubin, id.vars="data_transform", variable.name="model_spec", value.name="mean_rubin_error"))
  # summary.survival.model.summaries <- merge(summary.survival.model.summaries, mean_rubin, by=c("data_transform", "model_spec"), all=T)
  # 
  # # get AIC
  # load(paste0(main_dir, "AIC.rdata"))
  # mean_aic <- data.frame(meanAIC)
  # setnames(mean_aic, names(meanAIC))
  # mean_aic$data_transform <- rownames(mean_aic)
  # mean_aic <- data.table(melt(mean_aic, id.vars="data_transform", variable.name="model_spec", value.name="mean_aic"))
  # summary.survival.model.summaries <- merge(summary.survival.model.summaries, mean_aic, by=c("data_transform", "model_spec"), all=T)
  # 
  
  # get RMSE from OOS file
  load(paste0(main_dir, "compiled_rmse.rdata"))
  oos_rmse <- data.frame(compiled_rmse)
  setnames(oos_rmse, colnames(compiled_rmse))
  oos_rmse$data_transform <- rownames(oos_rmse)
  oos_rmse <- data.table(melt(oos_rmse, id.vars="data_transform", variable.name="model_spec", value.name="oos_rmse"))
  summary.survival.model.summaries <- merge(summary.survival.model.summaries, oos_rmse, by=c("data_transform", "model_spec"), all=T)
  
  
  #get ranking as per RMSE
  ranking <- unique(summary.survival.model.summaries[order(oos_rmse), list(data_transform, model_spec, oos_rmse)])
  ranking[, ranking:= as.numeric(rownames(ranking))]
  summary.survival.model.summaries <- merge(summary.survival.model.summaries, ranking, by=c("data_transform", "model_spec", "oos_rmse"), all=T)
  summary.survival.model.summaries <- summary.survival.model.summaries[order(ranking, data_transform, model_spec)]
  summary.survival.model.summaries[, lower:= beta-1.96*se]
  summary.survival.model.summaries[, upper:= beta+1.96*se]
  
  save(summary.survival.model.summaries, file=paste0(main_dir, "all_model_summaries.rdata"))

}else{
  load(paste0(main_dir, "all_model_summaries.rdata"))
  ranking <- unique(summary.survival.model.summaries[order(oos_rmse), list(data_transform, model_spec, oos_rmse)])
  ranking[, ranking:= as.numeric(rownames(ranking))]
}

summary.survival.model.summaries[,lower:= beta-1.96*se]
summary.survival.model.summaries[,upper:=beta+1.96*se]
summary.survival.model.summaries[, effect_mean:=1-exp(beta)]
summary.survival.model.summaries[, effect_lower:= 1-exp(lower)]
summary.survival.model.summaries[, effect_upper:= 1-exp(upper)]


##-----------------------------------------------
## Get means across imputations for all datasets,
## for plotting and predicion
##-----------------------------------------------

data_list <- lapply(colnames(data.for.survival), function(col){
  print(col)
  
  this_transform <- data.for.survival[, col]
  this_transform <- rbindlist(this_transform)
  if ("event_time_debiased" %in% names(this_transform)) setnames(this_transform, "event_time_debiased", "event_time")
  this_transform <- this_transform[, list(event_time=mean(event_time), agesero=mean(agesero), 
                                          observed_survival=mean(observed_survival), spvl_model=mean(spvl_model),
                                          spvl_fraser=mean(spvl_fraser), event_num=mean(event_num)), by="patient_id"]
  
})
names(data_list) <- colnames(data.for.survival)


##-------------------------------
## Visualize
##-------------------------------

##generate data from which to predict
to_predict_template <- data.table(expand.grid(event_num=0,
                                     agesero=seq(15, 100,5),
                                     spvl=seq(-0.5, 10, 0.5)))

ages_to_keep <- c(15, 30, 45, 60)

#prediction function
predict_results <- function(ranking_list, title="Predictions by Model", show_legend=T){

  to_plot <- NULL
  idx <-1
  
  for (rank in ranking_list){
    #find betas 
    this_model <- summary.survival.model.summaries[ranking==rank]
    betas <- this_model$beta
    names(betas) <- this_model$covariate
    names(betas)[[1]] <- "Intercept"
    betas <- as.list(betas)
    
    spvl_type <- ifelse("spvl_model" %in% names(betas), "spvl_model",
                        ifelse("spvl_fraser" %in% names(betas), "spvl_fraser",
                               NA))
    
    for (val in c("agesero", "spvl_model", "agesero:spvl_model", "agesero:event_num1", "spvl_model:event_num1", "agesero:spvl_model:event_num1")){
      if (is.null(betas[[val]])) betas[[val]] =0
    }
    
    to_predict <- copy(to_predict_template)
    
    #predicted values
    to_predict[, time_to_event:= exp(betas$Intercept + event_num*betas$event_num1 +
                                  agesero*betas$agesero +
                                  spvl*betas[[spvl_type]]+
                                  agesero*spvl*betas[[paste0("agesero:", spvl_type)]] +
                                  agesero*event_num*betas[["agesero:event_num1"]] +
                                  spvl*event_num*betas[[paste0(spvl_type, ":event_num1")]] +
                                  agesero*spvl*event_num*betas[[paste0("agesero:", spvl_type, ":event_num1")]])
         ]
    to_predict[, ranking:=rank]
    to_predict <- to_predict[agesero %in% ages_to_keep]
    to_predict[, binned_age:= agesero]
    to_predict[, category:="modeled"]
    
    #load and prep real data
    this_transform <- ranking[ranking==rank, data_transform]
    data <- data_list[[this_transform]]
    data[, binned_age:= cut(agesero, breaks=c(0,30, 45, 60, 100), labels=c(15,30,45,60))]
    data <- data[, list(category="data", event_num=as.numeric(as.character(event_num)), agesero, binned_age, spvl=spvl_model, ranking=rank, time_to_event=exp(observed_survival))]
  
    to_predict <- rbind(to_predict, data, use.names=T)
    to_predict[, agesero:= as.factor(agesero)]
    to_predict[, data_transform := this_transform]
    
    to_plot[[idx]] <- to_predict
    idx <- idx+1
  }

  to_plot <- rbindlist(to_plot)
  to_plot[, ranking:=factor(ranking, levels=ranking_list)]
  to_plot <- to_plot[spvl>=2 & spvl<=6.5 & event_num==0]
  
  result <- ggplot(data=to_plot, aes(x=spvl, y=time_to_event)) +
      geom_point(data=to_plot[category=="data"], aes(color=binned_age), alpha=0.5) +
      geom_line(data=to_plot[category=="modeled"], aes(color=agesero), size=1) +
      facet_grid(.~ranking) +
      scale_color_discrete(name="Age at\nSeroconversion")+
      labs(title=title,
           x="SPVL (log10 units/mL)",
           y="Time to Death (Years)")
  
  if (!show_legend) result <- result+theme(legend.position="none")
  
  return(result)
}

#pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/results_1.pdf", width=14, height=8)

one <- predict_results(33:34, title="Non-Imputed Full")
two <- predict_results(13:16, title="Non-Imputed Pre-96")
three <- predict_results(20, title="Imputed Full")
four <- predict_results(1:3, title="Imputed Pre-96")

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/results_1.pdf", width=10, height=6)
  print(one)
graphics.off()

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/results_2.pdf", width=14, height=8)
print(two)
graphics.off()

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/results_3.pdf", width=5, height=5)
print(three)
graphics.off()

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/results_4.pdf", width=11, height=7)
print(four)
graphics.off()


print(two)
print(three)
print(four)


#multiplot(one,two,three,four, layout=matrix(c(1,2,3,4,1,2,3,4,0,2,3,0, 0,2,0,0), nrow=4))

#plot age at seroconversion in the three datasets
reported_transforms <- unique(ranking[ranking %in% c(1,13,20,33), data_transform])
reported_data <- lapply(reported_transforms, function(this_transform){
                    print(this_transform)
                    this_data <- data_list[[this_transform]]
                    this_data[, data_transform:=this_transform]
                    return(this_data)
})

reported_data <- rbindlist(reported_data)
reported_data <- reported_data[, list(patient_id, agesero, observed_survival, spvl_model, event_num=as.numeric(as.character(event_num)), data_transform)]

#give data_transforms better names
reported_data[, data_transform:= factor(data_transform, labels=c("Non-Imputed Full", "Non-Imputed Pre-96", "Imputed Full", "Imputed Pre-96"))]
reported_data[, data_transform:= factor(data_transform, levels=c("Non-Imputed Pre-96", "Non-Imputed Full", "Imputed Pre-96","Imputed Full"))]

agesero <- ggplot(reported_data[event_num==0], aes(x=agesero, y=observed_survival)) +
            geom_point(alpha=0.5) +
            facet_grid(~data_transform) +
            labs(x="Age at Seroconversion", y="log(Survival Time)")

spvl <- ggplot(reported_data[event_num==0], aes(x=spvl_model, y=observed_survival)) +
          geom_point(alpha=0.5) +
          facet_grid(~data_transform)+
          labs(x="SPVL", y="log(Survival Time)")
        

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/results_data.pdf", width=11, height=7)
multiplot(spvl, agesero)
graphics.off()

#multiplot(agesero, spvl, layout=matrix(c(1,2,0,0, 1,2,0,0, 1,2,0,0), nrow=4))


#plot secular trend of spvl
surv[, serocon_year := as.numeric(format(serocon_date, "%Y"))]
surv[, serocon_month:= as.numeric(format(serocon_date, "%m"))]
surv[, serocon_time:= as.Date(paste0(serocon_year, "-", serocon_month, "-15"), origin="1970-01-01")]
surv[, mean_spvl:=mean(spvl_model, na.rm=T), by="serocon_time"]
mean_times <- unique(surv[,list(serocon_time, mean_spvl)])
mean_times <- mean_times[order(serocon_time)]

secular_trend <- gam(spvl_model~serocon_date, data=surv)

dateseq <- seq(as.Date("1983-07-15"), as.Date("2014-01-28"), by="month")
new_data <- data.table(serocon_date=dateseq)
new_data[, predicted:= predict(secular_trend, newdata=new_data)]

ggplot(new_data, aes(x=serocon_date, y=predicted)) +
      geom_line()


secular <- ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
            geom_point(alpha=0.3) +
            geom_smooth(size=2) +
            labs(title="Secular Trend of SPVL: Conditional Mean",
                 x="Seroconversion Date",
                 y="SPVL")

secular_hex <- ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
  geom_hex(bins=50)+ 
  labs(title="Secular Trend of SPVL:Density",
       x="Seroconversion Date",
       y="SPVL") +
  theme(legend.position="none")

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/presentation_figures/secular_trend.pdf", width=14, height=8)

  print(secular)


graphics.off()


ggplot(surv[event_type=="death"], aes(x=agesero, y=event_timeNew)) +
  geom_point()

ggplot(surv[event_type=="death"], aes(x=spvl_model, y=log(event_timeNew))) +
  geom_point()



#look at age distribution of deaths among entire dataset
age_test <- data_list[["3.1-1-0-0"]]
age_test[, binned_age:= cut(agesero, breaks=c(0,20, 40, 60, 80, 100), labels=c(0,20,40,60,80))]

ggplot(age_test[as.numeric(as.character(event_num))==0], aes(x=spvl_model, y=exp(observed_survival), color=binned_age)) +
      geom_point(aes(color=binned_age), alpha=0.2)

#look at predicted results when only age is included in the model
age_only <- summary.survival.model.summaries[ranking==13]
betas <- age_only$beta
names(betas) <- age_only$covariate
names(betas)[[1]] <- "Intercept"
betas <- as.list(betas)

#predict (TODO: make this modifiable)
to_predict <- copy(to_predict_template)
to_predict[, predicted_time:= exp(betas$Intercept + event_num*betas$event_num1 + agesero*betas$agesero)]

ages_to_keep <- c(20, 40, 60, 80)
to_predict <- to_predict[agesero %in% ages_to_keep]
to_predict[, agesero:= as.factor(agesero)]

ggplot(to_predict[spvl>=2.5 & spvl<=7.5], aes(x=spvl, y=predicted_time, group=agesero)) +
  geom_line(aes(color=agesero))

#look at different mar events vs event time: linear in log space?
ggplot(surv[event_type=="mar"], aes(x=spvl_model, y=log(event_timeNew))) +
      geom_point(aes(color=event_detail), alpha=0.3)

test <- surv[!is.na(spvl_model) & event_type=="mar"]
test_lm <- lm(log(event_timeNew)~spvl_model,data=test)
test$predict <- predict(test_lm, data=test)
ggplot(test, aes(x=spvl_model)) +
      geom_point(aes(y=log(event_timeNew), color=event_detail), alpha=0.1)

#find mean within each spvl bin
test[, spvl_bin := cut(spvl_model, breaks=-1:9, labels=-1:8)]
test[, mean_surv:= mean(event_timeNew), by="spvl_bin"]

ggplot(test[pre_1996==T], aes(x=spvl_model)) +
  geom_point(aes(y=log(event_timeNew), color=event_detail), alpha=0.3)

#look at different mar events vs event time: linear in log space?
ggplot(surv[event_type=="mar"], aes(x=agesero, y=log(event_timeNew))) +
  geom_point(aes(color=event_detail), alpha=0.1)


#look at secular trend
test[, year_enrolled := as.numeric(format(enroll_date, "%Y"))]
test[, year_evented:= as.numeric(format(event_date, "%Y"))]
ggplot(test, aes(x=spvl_model)) +
  geom_point(aes(y=log(event_timeNew), color=event_detail), alpha=0.1)+
  facet_wrap(~year_evented)


#plot this all nice
pdf(paste0(main_dir, "mar_plots.pdf"), width=15, height=8)

#the whole cloud
full <- ggplot(surv[event_type=="mar"], aes(x=spvl_model, y=log(event_timeNew))) +
          geom_point(aes(color=event_detail), alpha=0.3) +
          labs(title="SPVL vs log Time to Event, by MAR Type",
               x="log10(spvl)",
               y="log(years to event)")
print(full)

secular_by_enroll <- ggplot(test, aes(x=spvl_model)) +
                    geom_point(aes(y=log(event_timeNew), color=event_detail), alpha=0.1)+
                    facet_wrap(~year_enrolled) +
                    labs(title="SPVL vs log Time to Event, by MAR Type and Year of Enrollment",
                         x="log10(spvl)",
                         y="log(years to event)")

print(secular_by_enroll)

secular_by_event <- ggplot(test, aes(x=spvl_model)) +
                    geom_point(aes(y=log(event_timeNew), color=event_detail), alpha=0.1)+
                    facet_wrap(~year_evented) +
                    labs(title="SPVL vs log Time to Event, by MAR Type and Year of Event",
                         x="log10(spvl)",
                         y="log(years to event)")

print(secular_by_event)

#just pre-1996
part <- ggplot(surv[event_type=="mar" & pre_1996==T], aes(x=spvl_model, y=log(event_timeNew))) +
        geom_point(aes(color=event_detail), alpha=0.5) +
        labs(title="SPVL vs log Time to Event, by MAR Type, Pre-1996 Only",
             x="log10(spvl)",
             y="log(years to event)")
print(part)


#also look at non-mar
full_non_mar <- ggplot(surv[event_type!="mar"], aes(x=spvl_model, y=log(event_timeNew))) +
                geom_point(aes(color=event_detail), alpha=0.5) +
                labs(title="SPVL vs log Time to Event, by Event Type",
                     x="log10(spvl)",
                     y="log(years to event)")
print(full_non_mar)

part_non_mar <- ggplot(surv[event_type!="mar" & pre_1996==T], aes(x=spvl_model, y=log(event_timeNew))) +
                geom_point(aes(color=event_detail), alpha=0.5) +
                labs(title="SPVL vs log Time to Event, by Event Type, Pre-1996 Only",
                     x="log10(spvl)",
                     y="log(years to event)")
print(part_non_mar)


graphics.off()




##for number-plugging, find a breakdown of mar patients by mar type
load(paste0(main_dir, "alldata.rdata"))
alldata[, event_detail:=""]
alldata[event_type=="mar" & art_indic==1, event_detail:="art"]
alldata[event_type=="mar" & ltfu_date<art_start_date, event_detail:="ltfu"] #happens for two people
alldata[event_type=="mar" & event_detail=="" & ltfu_indic==1, event_detail:="ltfu"]
alldata[event_type=="mar" & event_detail=="", event_detail:="admin_censor"]
alldata[event_detail=="", event_detail:= event_type]
event_detail <- unique(alldata[, list(patient_id, event_detail)])
sexes <- unique(alldata[, list(patient_id, sex, inf_mode)])


#paste on to our final dataset
surv <- merge(surv, event_detail, by="patient_id", all.x=T)
surv <- merge(surv, sexes, by="patient_id", all.x=T)
surv[, unlog_spvl_model:= 10^(spvl_model)]
surv[, unlog_spvl_fraser:= 10^(spvl_fraser)]





















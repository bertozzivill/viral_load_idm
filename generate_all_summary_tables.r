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
library(xtable)


main_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/"

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
    
    #find data transform and imputation count
    pattern <- "(.*)-imp_count_([0-9]*)"
    pattern_match <- str_match(transform_imp, pattern)
    data_transform <- pattern_match[1,2]
    # data_transform <- gsub("FALSE", "0", data_transform)
    # data_transform <- gsub(" TRUE", "1", data_transform)
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
  
  # get RMSE from OOS file
  load(paste0(main_dir, "validation/compiled_rmse.rdata"))
  oos_rmse <- data.frame(compiled_rmse)
  setnames(oos_rmse, colnames(compiled_rmse))
  oos_rmse$data_transform <- rownames(oos_rmse)
  oos_rmse <- data.table(melt(oos_rmse, id.vars="data_transform", variable.name="model_spec", value.name="oos_rmse"))
  summary.survival.model.summaries <- merge(summary.survival.model.summaries, oos_rmse, by=c("data_transform", "model_spec"), all=T)


  # get ranking as per RMSE
  ranking <- unique(summary.survival.model.summaries[order(oos_rmse), list(data_transform, model_spec, oos_rmse)])
  ranking[, ranking:= as.numeric(rownames(ranking))]
  summary.survival.model.summaries <- merge(summary.survival.model.summaries, ranking, by=c("data_transform", "model_spec", "oos_rmse"), all=T)
  summary.survival.model.summaries <- summary.survival.model.summaries[order(ranking, data_transform, model_spec)]

  ## summary stats
  
  summary.survival.model.summaries[,lower:= beta-1.96*se]
  summary.survival.model.summaries[,upper:=beta+1.96*se]
  summary.survival.model.summaries[, effect_mean:=1-exp(beta)]
  summary.survival.model.summaries[, effect_lower:= 1-exp(lower)]
  summary.survival.model.summaries[, effect_upper:= 1-exp(upper)]
  
  save(summary.survival.model.summaries, file=paste0(main_dir, "all_model_summaries.rdata"))

}else{
  load(paste0(main_dir, "all_model_summaries.rdata"))
  #ranking <- unique(summary.survival.model.summaries[order(oos_rmse), list(data_transform, model_spec, oos_rmse)])
  ranking <- unique(summary.survival.model.summaries[,list(data_transform, model_spec)])
  ranking[, ranking:= as.numeric(rownames(ranking))]
}


##-----------------------------------------------
## Get means across imputations for all datasets,
## for plotting and predicion
##-----------------------------------------------

data_list <- lapply(colnames(data.for.survival), function(col){
  print(col)
  
  this_transform <- data.for.survival[, col]
  this_transform <- rbindlist(this_transform)
  if ("event_time_debiased" %in% names(this_transform)) setnames(this_transform, "event_time_debiased", "event_time")
  this_transform <- this_transform[, list(event_time=mean(event_time), agesero=unique(agesero), 
                                          observed_survival=mean(observed_survival), spvl_model=mean(spvl_model),
                                          spvl_fraser=mean(spvl_fraser), event_num=unique(event_num)), by="patient_id"]
  
})
names(data_list) <- colnames(data.for.survival)
save(data_list, file=paste0(main_dir, "mean_imputed_data.rdata"))

##-------------------------------
## Visualize
##-------------------------------

## give data transforms and model specs beautiful names

name_transform <- function(this_transform){
  this_transform <- as.list(strsplit(this_transform, split="-")[[1]])
  names(this_transform) <- c("upper_bound", "debias", "pre_96", "no_impute")
  
  imp_type <- ifelse(this_transform$no_impute=="0", "Imputed", "Nonimputed")
  time <- ifelse(this_transform$pre_96=="1", "Pre 1996", "Full Time")
  bias_type <- ifelse(this_transform$debias=="1", "Debiased", "Nondebiased")
  imp_ub <- ifelse(this_transform$upper_bound=="2.9" & this_transform$no_impute=="1", "none", this_transform$upper_bound)
  
  describe_transform <- data.table(imp_type=imp_type, time=time, bias_type=bias_type, imp_ub=imp_ub)

return(describe_transform)
}

name_spec <- function(this_spec){
  this_spec <- as.list(strsplit(this_spec, split="-")[[1]])
  names(this_spec) <- c("spvl_method", "interaction_type", "include.age", "age.type")
  
  spvl_type <- ifelse(this_spec$spvl_method=="none", "none", ifelse(this_spec$spvl_method=="spvl_model", "Nonlinear", "Geometric"))
  interaction_type <- ifelse(this_spec$interaction_type=="none", "none",
                             ifelse(this_spec$interaction_type=="two_way", "Two Way", "Three Way"))
  include_age <- ifelse(this_spec$include.age=="FALSE", F, T)
  age_type <- ifelse(this_spec$age.type=="cont", "Continuous Age",
                     ifelse(this_spec$age.type=="bin_10", "10yr Age Bins",
                            ifelse(this_spec$age.type=="quint", "Age Quintiles", "No Age")))
  
  describe_spec <- data.table(spvl_type=spvl_type, interaction_type=interaction_type, include_age=include_age, age_type=age_type)
  
  return(describe_spec)

}

full_transform <- rbindlist(lapply(ranking$data_transform, name_transform))
full_transform[, ranking:= as.numeric(rownames(full_transform))]
full_spec <- rbindlist(lapply(as.character(ranking$model_spec), name_spec))
full_spec[, ranking:= as.numeric(rownames(full_spec))]

ranking <- merge(ranking, full_transform, by="ranking", all=T)
ranking <- merge(ranking, full_spec, by="ranking", all=T)
ranking[, imp_ub:= ifelse(imp_ub=="none", "none",
                          ifelse(imp_ub=="2.9", "18",
                                 ifelse(imp_ub=="3", "20", "22")))]

ranking[, full_transform:= paste0(time, "-", bias_type, "-", imp_type,   ifelse(imp_ub=="none", "", paste0("-UB ", imp_ub)))]

ranking[, full_spec:= ifelse(spvl_type=="none" & include_age==F, "Null",
                             ifelse(spvl_type=="none" & include_age==T, paste(age_type, "Only"),
                                    paste0(spvl_type, "-", ifelse(include_age==F, "SPVL Only", 
                                                                  ifelse(interaction_type=="none", "Central", interaction_type)))))]

###heatmap of data tranforms and model specs, with color=rmse

# for_heatmap <- copy(ranking)
# 
# setnames(for_heatmap, c("oos_rmse", "ranking"), c("RMSE", "Ranking"))
# for_heatmap[imp_ub=="none", imp_ub:="0"] #so ranking works
# for_heatmap <- for_heatmap[order(-time, -bias_type, -imp_type, imp_ub)]
# 
# for_heatmap[, full_transform:= factor(full_transform, levels=rev(unique(for_heatmap$full_transform)))]
# 
# for_heatmap[, full_spec:= factor(full_spec, levels=c( "Null", "Age Only", "Geometric-SPVL Only", "Geometric-Central", "Geometric-Two Way", "Geometric-Three Way",
#                                                       "Nonlinear-SPVL Only", "Nonlinear-Central",  "Nonlinear-Two Way", "Nonlinear-Three Way"))]
# 
# pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/heatmap.pdf", width=11, height=9)
# 
# heatmap <- ggplot(for_heatmap, aes(x=full_spec, y=full_transform)) +
#             geom_tile(aes(fill=RMSE)) +
#             geom_text(color="grey", aes(fill = RMSE, label = paste0(round(RMSE, 2), "\n", "(", Ranking, ")" ))) +
#             scale_fill_gradient(low = "blue", high = "red") +
#             theme(axis.text.x=element_text(angle = 45, hjust=1))+
#             labs(title="RMSE by Data Transform and Model Specification",
#                  y="Data Transform",
#                  x="Model Specification")
# 
# print(heatmap)
# 
# graphics.off()


# ##find mean RMSE's over each possible type of action
# summary_rmse <- melt(ranking, id.vars=c("ranking", "oos_rmse"), measure.vars=c("imp_type", "time", "bias_type", "imp_ub",
#                                                                              "spvl_type", "interaction_type", "include_age"))
# summary_rmse[, type:= ifelse(variable %in% c("imp_type", "time", "bias_type", "imp_ub"), "data_transform", "model_spec")]
# summary_rmse[, mean_rmse:= mean(oos_rmse), by=c("variable", "value")]
# unique(summary_rmse[, list(type, variable, value, mean_rmse)])


##line plots of predicted results and data


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
  to_plot <- merge(to_plot, ranking[, list(ranking, full_spec, full_transform)], by="ranking", all.x=T)
  to_plot[, full_spec:=gsub("Nonlinear-", "", full_spec)]
  to_plot[, full_transform:= gsub("Pre 1996-Nondebiased-", "", full_transform)]
  to_plot[, full_name:= paste(full_transform, full_spec, sep="-")]
  to_plot <- to_plot[order(ranking)]
  ordered_names <- unique(to_plot$full_name)
  to_plot[, full_name:=factor(full_name, levels=ordered_names)]
  to_plot <- to_plot[spvl>=2 & spvl<=6.5 & event_num==0]
  
  result <- ggplot(data=to_plot, aes(x=spvl, y=time_to_event)) +
      geom_point(data=to_plot[category=="data"], aes(color=binned_age), alpha=0.5) +
      geom_line(data=to_plot[category=="modeled"], aes(color=agesero), size=1) +
      facet_grid(.~full_name) +
      scale_color_discrete(name="Age at\nSeroconversion")+
      labs(title=title,
           x="SPVL (log10 units/mL)",
           y="Time to Death (Years)")
  
  if (!show_legend) result <- result+theme(legend.position="none")
  
  return(result)
}

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/predictions.pdf", width=11, height=4)

  print(predict_results(1:5))

graphics.off()


### better look at % change for agesero and spvl
perc_table <- summary.survival.model.summaries[ranking %in% 1:10]
perc_table <- perc_table[covariate %in% c("agesero", "spvl_model"), list(ranking, data_transform, model_spec, oos_rmse, covariate, effect_mean, effect_lower, effect_upper)]
perc_table[,effect_mean:= ifelse(covariate=="agesero", effect_mean*1000, effect_mean*100)]
perc_table[,effect_lower:= ifelse(covariate=="agesero", effect_lower*1000, effect_lower*100)]
perc_table[,effect_upper:= ifelse(covariate=="agesero", effect_upper*1000, effect_upper*100)]

### regression tables

make_reg_table <- function(rankings){
  reg_table <- summary.survival.model.summaries[ranking %in% rankings]
  reg_table <- reg_table[, list(ranking, covariate, beta, lower, upper)]
  reg_table[, uncert:= paste0("(", round(lower,2), ", ", round(upper,2), ")")]
  reg_table[, beta:= as.character(round(beta,2))]
  reg_table[, c("lower", "upper"):=NULL]
  reg_table <- merge(reg_table, ranking[, list(ranking, full_spec, full_transform)], by="ranking", all.x=T)
  reg_table <- melt(reg_table, id.vars = c("ranking", "full_transform", "full_spec", "covariate"))
  reg_table[, full_spec:= gsub("Nonlinear-", "", full_spec)]
  reg_table[, full_transform:= gsub("Pre 1996-Nondebiased-", "", full_transform)]
  
  reg_table <- data.table(dcast(reg_table, ranking+full_transform+covariate+variable~full_spec))
  
  #format rows and columns
  setnames(reg_table, c("ranking", "full_transform"), c("Ranking", "Data Transform"))
  setcolorder(reg_table, c("Ranking", "Data Transform", "covariate", "variable", "SPVL Only", "Central", "Two Way"))
  
  reg_table[, covariate:= factor(covariate, labels=c("Intercept", "Age",
                                                      "Age*SPVL", "I(AIDS)", "SPVL"))]
  
  reg_table[, covariate:= factor(covariate, levels=c("Intercept", "I(AIDS)", "Age",
                                                     "SPVL", "Age*SPVL"))]
  reg_table <- reg_table[order(Ranking, covariate)]
  
  reg_table[, covariate:=as.character(covariate)]
  reg_table[variable=="uncert", covariate:=NA]
  reg_table[, variable:=NULL]
  reg_table[covariate!="Intercept" | is.na(covariate), Ranking:=NA]
  reg_table[covariate!="Intercept" | is.na(covariate)][["Data Transform"]] <- NA
  setnames(reg_table, "covariate", "Model Term")
  
  return(reg_table)
}

top_5 <- make_reg_table(1:5)
print(xtable(top_5, caption="Top Five Regression Outputs", label="reg_table", include.rownames=F))


### what's the mean difference in rmse between geometric and nonlinear spvl?
test_diff <- ranking[model_spec %like% "spvl"]
test_diff[, spvl_method:= ifelse(full_spec %like% "Nonlinear", "Nonlinear", "Geometric")]
test_diff[, short_spec:= gsub("(Nonlinear|Geometric)-", "", full_spec)]
test_diff[, full_model:=paste(full_transform, short_spec, sep="-")]
test_diff <- data.table(dcast(test_diff, full_model~spvl_method, value.var="oos_rmse"))
test_diff[, diff:= Geometric-Nonlinear]

## between pre-96 and full timeseries for those two?
test_diff[, pre_96 := ifelse(full_model %like% "Pre 1996", T, F)]
year_diff <- test_diff[, list(Geometric=mean(Geometric), Nonlinear=mean(Nonlinear)), by="pre_96"]
year_diff[, diff:= Geometric-Nonlinear]

###plot age at seroconversion and spvl in the four datasets
reported_transforms <- unique(ranking[bias_type=="Nondebiased" & imp_ub %in% c("none", "18"), data_transform])
reported_data <- lapply(reported_transforms, function(this_transform){
                    print(this_transform)
                    this_data <- data_list[[this_transform]]
                    this_data[, data_transform:=this_transform]
                    this_data[, binned_age:=NULL]
                    return(this_data)
})

reported_data <- rbindlist(reported_data)
reported_data <- reported_data[, list(patient_id, agesero, observed_survival, spvl_model, event_num=as.numeric(as.character(event_num)), data_transform)]

#give data_transforms better names
reported_data[, data_transform:= factor(data_transform, labels=c("Imputed-Full Time", "Nonimputed-Full Time",  "Imputed-Pre 1996",   "Nonimputed-Pre 1996"))]
reported_data[, data_transform:= factor(data_transform, levels=c("Nonimputed-Pre 1996", "Nonimputed-Full Time", "Imputed-Pre 1996","Imputed-Full Time"))]

agesero <- ggplot(reported_data[event_num==0], aes(x=agesero, y=exp(observed_survival))) +
            geom_point(alpha=0.5) +
            facet_grid(~data_transform) +
            scale_y_log10(breaks=c(0,1,2,4,6,8,10,20))+
            labs(x="Age at Seroconversion", y="Survival Time(Years)")

spvl <- ggplot(reported_data[event_num==0], aes(x=spvl_model, y=exp(observed_survival))) +
          geom_point(alpha=0.5) +
          facet_grid(~data_transform)+
          scale_y_log10(breaks=c(0,1,2,4,6,8,10,20))+
          labs(x="SPVL(log10 units/mL)", y="Survival Time(Years)")
        

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/agesero_spvl_data.pdf", width=11, height=7)
  multiplot(spvl, agesero)
graphics.off()


###plot secular trend of spvl

secular <- ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
            geom_smooth(size=2, color="red") +
            labs(title="",
                 x="Seroconversion Date",
                 y="SPVL(log10 units/mL)")

secular_data <-  ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
                geom_point(alpha=0.3) +
                geom_smooth(size=2, color="red") +
                labs(title="Secular Trend of SPVL: Conditional Mean",
                     x="Seroconversion Date",
                     y="SPVL(log10 units/mL)")

secular_hex <- ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
  geom_hex(bins=50)+ 
  geom_smooth(size=2, color="red")+
  labs(title="Secular Trend of SPVL:Density",
       x="Seroconversion Date",
       y="SPVL(log10 units/mL)")

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/secular_trend.pdf", width=7, height=10)

  multiplot(secular_data, secular)

graphics.off()




















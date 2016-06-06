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

compute_summary <- F

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

## give data transforms and model specs beautiful names

this_transform <- "2.9-1-1-0"
return_full <- F

name_transform <- function(this_transform, return_full=T){
  this_transform <- as.list(strsplit(this_transform, split="-")[[1]])
  names(this_transform) <- c("upper_bound", "debias", "pre_96", "no_impute")
    
  if (return_full){
    this_name <- paste(ifelse(this_transform$no_impute=="0", "Imputed", "Nonimputed"),
                       ifelse(this_transform$pre_96=="1", "Pre 1996", "Full Time"),
                       ifelse(this_transform$debias=="1", "Debiased", "Nondebiased"),
                       ifelse(this_transform$upper_bound=="2.9" & this_transform$no_impute=="0", "", paste("UB", this_transform$upper_bound)),
                       sep="-")
  }else{
    this_name <- paste(ifelse(this_transform$no_impute=="0", "Imputed", "Nonimputed"),
                       ifelse(this_transform$pre_96=="1", "Pre 1996", "Full Time"),
                       sep="-")
  }
  
  if (str_sub(this_name, -1,-1)=="-") this_name <- substr(this_name, 1, nchar(this_name)-1)

return(this_name)
}

spec_transform <- function(this_spec, return_full=T){
  this_spec <- as.list(strsplit(this_spec, split="-")[[1]])
  names(this_spec) <- c("spvl_method", "interaction_type", "include.age")
  
  if (this_spec$spvl_method=="none"){
    this_name <- ifelse(this_spec$include.age=="FALSE", "Null", "Age Only")
  }else if (this_spec$include.age=="FALSE"){
    this_name <- "SPVL Only"
  }else{
    this_name <- ifelse(this_spec$interaction_type=="none", "Main",
                        ifelse(this_spec$interaction_type=="two_way", "Two Way", "Three Way"))
  }
  
  if (return_full & this_spec$spvl_method!="none"){
    this_name <- paste(this_name,
                       ifelse(this_spec$spvl_method=="spvl_model", "Nonlinear SPVL", "Geometric SPVL"),
                       sep="-")
  }
  
  return(this_name)
}

full_names <- unlist(lapply(ranking$data_transform, name_transform))
full_specs <- unlist(lapply(as.character(ranking$model_spec), spec_transform))

ranking[, full_transform:=full_names]
ranking[, full_spec:= full_specs]

###heatmap of data tranforms and model specs, with color=rmse

for_heatmap <- ranking[(full_spec %like% "Nonlinear" | !full_spec %like% "SPVL") & full_transform %like% "Debiased" &
                         (full_transform %like% "UB 3.1" | full_transform %like% "Nonimputed")]

short_names <- unlist(lapply(for_heatmap$data_transform, function(x){name_transform(x, return_full=F)}))
short_spec <- unlist(lapply(as.character(for_heatmap$model_spec), function(x){spec_transform(x, return_full=F)}))
for_heatmap[, short_transform:=short_names]
for_heatmap[, short_specs:= short_spec]

for_heatmap[, short_specs:= factor(short_specs, levels=c( "Null", "Age Only", "SPVL Only","Main", "Two Way", "Three Way"))]
for_heatmap[, short_transform:= factor(short_transform, levels=rev(c("Imputed-Pre 1996", "Nonimputed-Pre 1996", "Imputed-Full Time", "Nonimputed-Full Time")))]
setnames(for_heatmap, c("oos_rmse", "ranking"), c("RMSE", "Ranking"))

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/heatmap.pdf", width=7, height=9)

heatmap <- ggplot(for_heatmap, aes(x=short_specs, y=short_transform)) +
            geom_tile(aes(fill=RMSE)) +
            geom_text(aes(fill = RMSE, label = round(RMSE, 2))) +
            scale_fill_gradient(low = "blue", high = "red") +
            labs(title="RMSE by Data Transform and Model Specification",
                 y="Data Transform",
                 x="Model Specification")

ranking_heatmap <- ggplot(for_heatmap, aes(x=short_specs, y=short_transform)) +
                    geom_tile(aes(fill=Ranking)) +
                    geom_text(aes(fill = Ranking, label=Ranking)) +
                    scale_fill_gradient(low = "blue", high = "red") +
                    labs(title="Ranking by Data Transform and Model Specification",
                         y="Data Transform",
                         x="Model Specification")

multiplot(heatmap, ranking_heatmap)


graphics.off()


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

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/predictions.pdf", width=9, height=12)

imp_96 <- c(1:3, 10)
nonimp_96 <- 13:16
imp <- c(23, 28:30)
nonimp <- c(37:38, 42:43)
all_models <- c(imp_96, nonimp_96, imp, nonimp)


one <- predict_results(imp_96, title="A: Imputed-Pre 96")
two <- predict_results(nonimp_96, title="B: Nonimputed-Pre 96")
three <- predict_results(imp, title="C: Imputed-Full Time")
four <- predict_results(nonimp, title="D: Nonimputed-Full Time")

#multiplot(one,two,three,four, layout=matrix(c(1,2,3,4,1,2,3,4,0,2,3,0, 0,2,0,0), nrow=4))
multiplot(one,two,three,four)

graphics.off()


### better look at % change for agesero and spvl
perc_table <- summary.survival.model.summaries[ranking %in% nonimp]
perc_table <- perc_table[covariate %in% c("agesero", "spvl_model"), list(ranking, data_transform, model_spec, oos_rmse, covariate, effect_mean, effect_lower, effect_upper)]
perc_table[,effect_mean:= ifelse(covariate=="agesero", effect_mean*1000, effect_mean*100)]
perc_table[,effect_lower:= ifelse(covariate=="agesero", effect_lower*1000, effect_lower*100)]
perc_table[,effect_upper:= ifelse(covariate=="agesero", effect_upper*1000, effect_upper*100)]

### regression tables

make_reg_table <- function(this_transform, this_caption="", this_label=""){
  reg_table <- summary.survival.model.summaries[data_transform==this_transform & !model_spec %like% "fraser"]
  reg_table <- reg_table[, list(model_spec, covariate, beta)]
  spec_names <- unlist(lapply(as.character(reg_table$model_spec), function(x){spec_transform(x, return_full = F)}))
  reg_table[, model_spec:= spec_names]
  reg_table <- data.table(dcast(reg_table, covariate~model_spec))
  
  #format rows and columns
  setnames(reg_table, "covariate", "Covariate")
  setcolorder(reg_table, c("Covariate", "Null", "Age Only", "SPVL Only", "Main", "Two Way", "Three Way"))
  
  reg_table[, Covariate:= factor(Covariate, labels=c("Intercept", "Age at Seroconversion", "Age*I(AIDS)",
                                                      "Age*SPVL", "Age*SPVL*I(AIDS)", "I(AIDS)", "SPVL",
                                                      "SPVL*I(AIDS)"))]
  
  reg_table[, Covariate:= factor(Covariate, levels=c("Intercept", "I(AIDS)", "Age at Seroconversion",
                                                     "SPVL", "Age*SPVL", "Age*I(AIDS)", "SPVL*I(AIDS)",
                                                     "Age*SPVL*I(AIDS)"))]
  reg_table <- reg_table[order(Covariate)]
  this_table <- xtable(reg_table, caption=this_caption, label=this_label)
  return(this_table)
}

make_reg_table("3.1-1-1-0", this_caption="Imputed-Pre 1996 Regression Outputs", this_label="reg_imp_96")
make_reg_table("2.9-1-1-0", this_caption="Nonimputed-Pre 1996 Regression Outputs", this_label="reg_nonimp_96")
make_reg_table("3.1-1-0-0", this_caption="Imputed-Full Time Regression Outputs", this_label="imp")
make_reg_table("2.9-1-0-1", this_caption="Nonimputed-Full Time Regression Outputs", this_label="nonimp")



### what's the mean difference in rmse between geometric and nonlinear spvl?
test_diff <- ranking[model_spec %like% "spvl"]
test_diff[, spvl_method:= ifelse(full_spec %like% "Nonlinear", "Nonlinear", "Geometric")]
test_diff[, short_spec:= gsub("-(Nonlinear|Geometric) SPVL", "", full_spec)]
test_diff[, full_model:=paste(full_transform, short_spec, sep="-")]
test_diff <- data.table(dcast(test_diff, full_model~spvl_method, value.var="oos_rmse"))
test_diff[, diff:= Geometric-Nonlinear]



###plot age at seroconversion and spvl in the four datasets
reported_transforms <- unique(ranking[ranking %in% all_models, data_transform])
reported_data <- lapply(reported_transforms, function(this_transform){
                    print(this_transform)
                    this_data <- data_list[[this_transform]]
                    this_data[, data_transform:=this_transform]
                    return(this_data)
})

reported_data <- rbindlist(reported_data)
reported_data <- reported_data[, list(patient_id, agesero, observed_survival, spvl_model, event_num=as.numeric(as.character(event_num)), data_transform)]

#give data_transforms better names
reported_data[, data_transform:= factor(data_transform, labels=c("Nonimputed-Full Time", "Nonimputed-Pre 1996", "Imputed-Full Time", "Imputed-Pre 1996"))]
reported_data[, data_transform:= factor(data_transform, levels=c("Nonimputed-Pre 1996", "Nonimputed-Full Time", "Imputed-Pre 1996","Imputed-Full Time"))]

agesero <- ggplot(reported_data[event_num==0], aes(x=agesero, y=observed_survival)) +
            geom_point(alpha=0.5) +
            facet_grid(~data_transform) +
            labs(x="Age at Seroconversion", y="log(Survival Time)")

spvl <- ggplot(reported_data[event_num==0], aes(x=spvl_model, y=observed_survival)) +
          geom_point(alpha=0.5) +
          facet_grid(~data_transform)+
          labs(x="SPVL", y="log(Survival Time)")
        

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/agesero_spvl_data.pdf", width=11, height=7)
  multiplot(spvl, agesero)
graphics.off()

#multiplot(agesero, spvl, layout=matrix(c(1,2,0,0, 1,2,0,0, 1,2,0,0), nrow=4))


###plot secular trend of spvl
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



secular <- ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
            #geom_point(alpha=0.3) +
            geom_smooth(size=2, color="red") +
            labs(title="Secular Trend of SPVL: Conditional Mean",
                 x="Seroconversion Date",
                 y="SPVL")

secular_hex <- ggplot(surv, aes(x=serocon_date, y=spvl_model)) +
  geom_hex(bins=50)+ 
  geom_smooth(size=2, color="red")+
  labs(title="Secular Trend of SPVL:Density",
       x="Seroconversion Date",
       y="SPVL") +
  theme(legend.position="none")

pdf("C:/Users/abertozz/Documents/work/classes/thesis_spring2016/paper_figures/secular_trend.pdf", width=8, height=10)

  multiplot(secular, secular_hex)

graphics.off()




















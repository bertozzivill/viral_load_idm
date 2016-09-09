####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Make basic plots of raw data: age vs survival, as we know it

## Input: 
##                    
## Output: 
##
####################################################################################################

rm(list=ls())

library(ggplot2)
library(data.table)

main_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/"

load(paste0(main_dir, "mean_imputed_data.rdata")) #imputed data
load(paste0(main_dir, "prepped_data.rdata")) #non-imputed data

pdf(paste0(main_dir, "age_plots.pdf"), width=12, height=7)

censored <- surv[event_type=="mar", patient_id]

name_transform <- function(transform_name){
  transform_name <- as.list(strsplit(transform_name, split="-")[[1]])
  names(transform_name) <- c("upper_bound", "debias", "pre_96", "no_impute")
  
  imp_type <- ifelse(transform_name$no_impute=="0", "Imputed", "Nonimputed")
  time <- ifelse(transform_name$pre_96=="1", "Pre 1996", "Full Time")
  bias_type <- ifelse(transform_name$debias=="1", "Debias", "Nondebias")
  imp_ub <- ifelse(transform_name$upper_bound=="2.9" & transform_name$no_impute=="1", "Nonimputed", paste("UB", transform_name$upper_bound))
  
  describe_transform <- paste(imp_ub, bias_type, time, sep="-")
  
  return(describe_transform)
}

idx <- 1

for (idx in 1:length(data_list)){
  this_transform <- data_list[[idx]]
  transform_name <- names(data_list)[[idx]]
  
  #make transform_name understandable
  transform_name <- name_transform(transform_name)

  print(transform_name)
  
  this_transform[, event_time:= exp(observed_survival)]
  this_transform[, imputed_event:= patient_id %in% censored]
  this_transform[, agesero_bin10:= cut(agesero, breaks=c(15, 25, 35, 45, Inf), labels=c("15-25", "25-35", "35-45", "45+"))]
  this_transform[, agesero_quint:=cut(agesero, breaks=5)]
  this_transform[, event_type:= ifelse(event_num==0, "Death", "AIDS")]
  
  
  scatter <- ggplot(this_transform, aes(x=agesero, y=observed_survival)) +
    geom_point(aes(color=event_type), alpha=0.7)+
    guides(color=guide_legend(title="Event Type")) +
    theme(legend.position="left") +
    labs(title=transform_name,
         x="Age (Years)",
         y="log(Years to Event)")
  
  
  box_10_yr <- ggplot(this_transform, aes(x=agesero_bin10, y=observed_survival)) +
    geom_boxplot(aes(fill=event_type), alpha=0.7)+
    theme(legend.position="none") +
    labs(title="10 Year Age Bins",
         x="Age Bin",
         y="log(Years to Event)")
  
  box_quints <- ggplot(this_transform, aes(x=agesero_quint, y=observed_survival)) +
    geom_boxplot(aes(fill=event_type), alpha=0.7)+
    theme(legend.position="none") +
    labs(title="Age Quintiles",
         x="Age Quintile",
         y="log(Years to Event)")
  
  multiplot(scatter, box_10_yr, box_quints, cols=3)
  
  idx <- idx +1
}

graphics.off()


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

pdf(paste0(main_dir, "age_plots.pdf"), width=12, height=7)

name_transform <- function(transform_name){
  transform_name <- as.list(strsplit(transform_name, split="-")[[1]])
  names(transform_name) <- c("upper_bound", "debias", "pre_96", "no_impute", "age_type")
  
  imp_type <- ifelse(transform_name$no_impute=="0", "Imputed", "Nonimputed")
  time <- ifelse(transform_name$pre_96=="1", "Pre 1996", "Full Time")
  bias_type <- ifelse(transform_name$debias=="1", "Debias", "Nondebias")
  imp_ub <- ifelse(transform_name$upper_bound=="2.9" & transform_name$no_impute=="1", "Nonimputed", paste("UB", transform_name$upper_bound))
  age_type <- ifelse(transform_name$age_type=="cont", "Continuous Age",
                     ifelse(transform_name$age_type=="bin_10", "10yr Age Bins",
                            ifelse(transform_name$age_type=="quint", "Age Quintiles", "No Age")))
  
  describe_transform <- paste(imp_ub, bias_type, time, sep="-")
  
  return(describe_transform)
}

# Change the names of data_list to be nicely formatted, and *not* to include age
# so there are three datasets with the same name, repeated.
names(data_list) <- unlist(lapply(names(data_list), name_transform))

for (transform_name in unique(names(data_list))){
  print(transform_name)
  
  this_transform_idxs <- which(names(data_list)==transform_name)
  this_transform <- data_list[this_transform_idxs]
  
  this_transform <- lapply(this_transform, function(x) x[, event_type:= ifelse(event_num==0, "Death", "AIDS")])
  
  scatter <- ggplot(this_transform[[1]], aes(x=agesero, y=observed_survival)) +
    geom_point(aes(color=event_type), alpha=0.7)+
    guides(color=guide_legend(title="Event Type")) +
    theme(legend.position="left") +
    labs(title=transform_name,
         x="Age (Years)",
         y="log(Years to Event)")

  box_10_yr <- ggplot(this_transform[[2]], aes(x=agesero, y=observed_survival)) +
    geom_boxplot(aes(fill=event_type), alpha=0.7)+
    theme(legend.position="none") +
    labs(title= "10yr Age Bins",
         x="Age Bin",
         y="log(Years to Event)")

  # box_quints <- ggplot(this_transform[[3]], aes(x=agesero, y=observed_survival)) +
  #   geom_boxplot(aes(fill=event_type), alpha=0.7)+
  #   theme(legend.position="none") +
  #   labs(title="Age Quintiles",
  #        x="Age Quintile",
  #        y="log(Years to Event)")

  multiplot(scatter, box_10_yr, cols=2)

}

graphics.off()


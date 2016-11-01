## Test a cox model on the CASCADE survival data

library(data.table)
library(ggplot2)
library(survival)

main_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/cox_model/"

load(paste0(main_dir, "clean_data.rdata"))

#keep pre-96 only
pre_96 <- data[pre_1996==1]

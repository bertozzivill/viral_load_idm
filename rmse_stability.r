####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Analyze stability of RMSE across dataset
####################################################################################################

rm(list=ls())

library(data.table)
library(ggplot2)

main_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/"

load(paste0(main_dir, "validation/1/1/rmse.rdata"))


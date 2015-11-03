####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: See how well each model did by calculating rmse between model and testing dataset
## Input:  main_dir: Directory in which the file "alldata.rdata" (either original or the training data)
##                    lives. If not running cross-validation, this will be one of:
##                    "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
##                    "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"
##                    
## Output: Saves a dataset called "survival.model.output" to main_dir, which contains a matrix of 
##          "lm" elements from which we can later predict and compare to the testing dataset if 
##            we're validating. Also generates a number of plots for whatever model is deemed best
##            by the AIC method of model selection. 
##
## Run from within the "vectorize" folder!
####################################################################################################
rm(list=ls())
library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
library(Amelia)

#main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
#main_dir <- "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"

main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/cross_validation/1/1/"

#load models  and data
load(paste0(main_dir, "survival_model_output.rdata"))
load(paste0(main_dir, "test_data.rdata"))


##build rmse framework
this_lm <- survival.model.output[[1]]["lm", "spvl_model"][[1]]







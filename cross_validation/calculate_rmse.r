####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: See how well each model did by calculating rmse between model and testing dataset
## Input:  main_dir: Directory in which the file "alldata.rdata" (either original or the training data)
##                    lives. If not running cross-validation, this will be one of:
##                    "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
##                    "C:/Users/cselinger/Dropbox (IDM)/viral_load (1)/cascade/data"
##                    
## Output: Saves a dataset called "rmse.rdata" to main_dir, which contains a matrix of 
##            rmse values which we can later combine across all testing splits to find a best model.
##
## Run from within the "vectorize" folder!
####################################################################################################
rm(list=ls())
library(data.table)
library(ggplot2)
library(lme4)
library(reshape2)
#library(Metrics)

main_dir <- commandArgs()[3]
#main_dir <- "/ihme/scratch/users/abertozz/vl_cross_validation/1/1/"
 
 #load models  and data
  print("loading data")
  load(paste0(main_dir, "survival_model_output.rdata"))
  load(paste0(main_dir, "test_data.rdata"))
  
  #only need the events from test data
  test_data <- test_data[event_type!="mar"]
  test_data[, event_num:= factor(ifelse(event_type=="aids", 1, 2))]
  
  ##build rmse framework
  
  # survival.model.output is a list of length nrow(index.data.transform)*n_i, where n_i is the imputation count.
  # the first n_i elements of survival.model.output will correspond to the first element of index.data.transform, and so on.
  # Each element of survival.model.output is a matrix with width length(index.survival.model). Each element of that matrix 
  # corresponds to the model specified in index.survival.model, in order. The "lm" row of surival.model.output contains the lm
  # object we use to predict.
  
  imputation_count <- length(survival.model.output)/nrow(index.data.transform)
  
  #fill this array with the RMSE's you find.
  all_rmse <- array(dim=c(nrow(index.data.transform), nrow(index.survival.models), imputation_count))
  
  #name the array
  rownames(all_rmse)<-apply(index.data.transform,1,function(x) paste0(x,collapse="-"))
  colnames(all_rmse) <- apply(index.survival.models,1,function(x) paste0(x,collapse="-"))
  
  #loop through data transforms
  print("predicting and calculating rmse")
  for (transform_index in 1:nrow(index.data.transform)){
    start_index <- (transform_index-1)*imputation_count +1
    transform_info <- index.data.transform[transform_index,]
    #print(transform_info)
    
    this_transform <- survival.model.output[start_index:(start_index+imputation_count-1)]
    
    #loop through imputations, find RMSE for each model
    for (impute_index in 1:imputation_count){
      #print(paste("impute count", impute_index))
      
      this_imputation <- this_transform[[impute_index]]
      
      #loop through models
      for (model_index in 1:nrow(index.survival.models)){
        model_info <- index.survival.models[model_index,]
        #print(model_info)
        
        this_model <- this_imputation["lm", model_index][[1]]
        
        #bin the test data into the appropriate age groups for predictions.
        bins <- model_info$bins[[1]]
        if (bins!=0) test_data[,agebin:=cut(test_data$agesero, bins, include.lowest=TRUE)]
        
        #predict
        test_data$predicted <- predict(this_model, newdata=test_data)
        
        #remove na's
        rmse_data <- test_data[!is.na(predicted)]
        
        # find RMSE: if using a transform type with debias==1, use event_timeNew for the event time (it removed bias).
        # otherwise, use event_time.
        
        if (transform_info$debias ==1) event_colname <- "event_timeNew" else event_colname <- "event_time"
        #rmse <- rmse(rmse_data[[event_colname]], rmse_data$predicted)
        #until the cluster gets the Metrics package loaded, run rmse by hand:
        rmse <- sqrt(mean( (rmse_data[[event_colname]]-rmse_data$predicted)^2, na.rm=T) )
        
        #put rmse into the proper matrix: row==data transform, column=model specification, array==imputation
        all_rmse[transform_index, model_index, impute_index] <- rmse
        
      } #end model loop
    } #end imputation loop
  } #end transform loop
  
  #save this data!
  print("saving")
  save(all_rmse, file=paste0(main_dir, "rmse.rdata"))


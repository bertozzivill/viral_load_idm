##loop through directories, running models for that cross-validation each time

library(data.table)

rm(list=ls())
main_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/validation/"
validation <- T

iteration_count <- 10
split_count <- 5

#age_quints <- c(15.4, 28.2, 41, 53.8, 66.6, Inf)
age_quints <- c(0, 28.4, 36.3, 44.3, 52.2, Inf)

# prep cross validation
source("prep_cross_validation.r")

#run cross validation
print("running imputation and models, finding rmse")

for (iteration in 1:iteration_count){
  for (split in 1:split_count){
    print(paste("analyzing for iteration", iteration, "and split", split))
    source("_run_all_data_transform_model_selection.r")
  }
}

#once all these have run, calculate the best model of them all
print("finding best model")
source("find_best_model.r")


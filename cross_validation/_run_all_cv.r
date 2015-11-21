##loop through directories, running models for that cross-validation each time
# run from "cross_validation" folder!

library(data.table)

source("calculate_rmse.r")
setwd("../vectorize")
source("data.transform_model.selection.R")
main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/cross_validation/"

for (iteration in 1:8){
  for (split in 1:10){
    new_dir <-paste0(main_dir, iteration, "/", split, "/")
    print(paste("iteration", iteration, "split", split))
    #data.transform_model.selection(new_dir)
    calc_rmse(new_dir)
  }
}


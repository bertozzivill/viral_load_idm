##compile all rmse values; take the mean; find best model

library(data.table)

main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/cross_validation/"


#load rmse results
print("loading rmse")
for (iteration in 1:8){
  for (split in 1:10){
    new_dir <-paste0(main_dir, iteration, "/", split, "/")
    print(paste("iteration", iteration, "split", split))
    load(paste0(new_dir, "rmse.rdata"))
    apply(all_rmse, c(1,2), mean)
  }
}

compiled_rmse <- lapply(1:8, function(iteration){
      iteration_rmse <- lapply(1:10, function(split){
        new_dir <-paste0(main_dir, iteration, "/", split, "/")
        print(paste("iteration", iteration, "split", split))
        load(paste0(new_dir, "rmse.rdata"))
        split_rmse <- apply(all_rmse, c(1,2), mean)
        return(split_rmse)
      })
      #take the mean of this list
      iteration_count <- length(iteration_rmse)
      iteration_rmse <- Reduce("+", iteration_rmse) / iteration_count
      return(iteration_rmse)
})

#take the mean of THIS list
compiled_count <- length(compiled_rmse)
compiled_rmse <- Reduce("+", compiled_rmse) / compiled_count

#find the indices of the minimum RMSE
best_model <- which(compiled_rmse==min(compiled_rmse), arr.ind=T)
best_row <- rownames(compiled_rmse)[[best_model[,"row"]]]
best_col <- colnames(compiled_rmse)[[best_model[,"col"]]]

print(paste("best model has data transform", best_row, "and model specification", best_col, "!"))


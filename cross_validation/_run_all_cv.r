##loop through directories, running models for that cross-validation each time
# run from "cross_validation" folder!

library(data.table)

rm(list=ls())
main_dir <- ("/ihme/scratch/users/abertozz/vl_cross_validation/")
setwd("../vectorize/")

## define qsub function
qsub <- function(code, name="abertozz", arguments=NULL, hold=NULL, shell="r_shell.sh", slots=1,
                 error_dir="/share/temp/sgeoutput/abertozz/errors",
                 output_dir="/share/temp/sgeoutput/abertozz/output") {
  
  # format arguments and hold IDs
  if(!is.null(arguments)) arguments <- paste(paste("\"", arguments, "\"", sep=""), collapse=" ")
  if(!is.null(hold) & length(hold)>1) hold <- paste(hold, collapse=",")
  
  # construct and submit qsub command and return the job ID
  x <- paste("/usr/local/bin/SGE/bin/lx-amd64/qsub -cwd",
             "-N", name,
             "-e", error_dir,
             "-o", output_dir,
             if (slots>1) paste("-pe multi_slot", slots),
             if (!is.null(hold)) paste("-hold_jid", hold),
             shell, code,
             if (!is.null(arguments)) arguments)
  id <- system(x, intern = T)
  return(as.numeric(strsplit(id, " ")[[1]][3]))
}

## ONLY uncomment this if you intend to overwrite the imputed dataset that was used for paper publication (not recommended!)
# # prep cross validation
# prep_jid <- qsub(code="../cross_validation/prep_cross_validation.r",
#                  name="prep_cv",
#                  arguments=c(main_dir),
#                  slots=10)

#run cross validation
cv_jid <- lapply(1:10, function(iteration){
        iter_jid <- lapply(1:10, function(split){
          new_dir <-paste0(main_dir, iteration, "/", split, "/")
          
          print(paste("iteration", iteration, "split", split))
          
          #submit job/get id for main modeling step
          split_jid <- qsub(code = "_run_all_data_transform_model_selection.r",
                       name= paste0("cv_", iteration, "_", split),
                       arguments = c(new_dir),
                       slots=10,
                       hold=prep_jid)
          
          #submit rmse job once the main modeling job is done
          rmse_jid <- qsub(code = "../cross_validation/calculate_rmse.r",
                           name= paste0("rmse_", iteration, "_", split),
                           arguments = c(new_dir),
                           slots=10,
                           hold=split_jid)
        return(rmse_jid)
          
        })
        iter_jid <- unlist(iter_jid)
        return(iter_jid)
})

#once all these have run, calculate the best model of them all
print("finding best model")
qsub(code = "../cross_validation/find_best_model.r",
     name= "find_best_model",
     arguments = c(main_dir),
     slots=5,
     hold=unlist(cv_jid))


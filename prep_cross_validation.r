#################################################################################################
## prep_cross_validation.r

## DESCRIPTION: Split the dataset into 10 parts. In 10 different folders, save each part
##              (testing dataset), plus the remaining 9 parts (training dataset). Do this ten 
##              times, for a total of 100 different validation datasets.
## 
#################################################################################################

#setup
set.seed(strtoi("0xCAFE"))

library(foreign)
library(data.table)
library(reshape2)
require(bit64)

#delete extant files
lapply(1:iteration_count, function(x) unlink(paste0(main_dir, x), recursive=T))

#make new files
lapply(1:iteration_count, function(x){
  lapply(1:split_count, function(y){
    dir.create(paste0(main_dir, x, "/", y), recursive = T)
  })
})


#load all data 
load(paste0(main_dir, "/../prepped_data.rdata"))
patients <- unique(surv$patient_id)
split_size <- ceiling(length(patients)/split_count)

summary <- data.table(iteration=numeric(), split=numeric(), type=character(), mar=numeric(), aids=numeric(), death=numeric()) 

#rename data because the subsets must be named "surv"
fulldata <- copy(surv)
setkeyv(fulldata, "patient_id")

for (iteration in 1:iteration_count){
  print(paste("iteration", iteration))
  #split into ten groups of patients
  shuffled <- sample(patients)
  startval <- 1
  
  for (split in 1:split_count){
    
    print(paste("split", split))
    endval <- split_size*split
    test_patients <- shuffled[startval:endval]
    train_patients <- shuffled[!shuffled %in% test_patients]
    test_data <- fulldata[J(test_patients)]
    surv <- fulldata[J(train_patients)]
    
    #spot checks: tell me how many mar vs AIDS vs death events there are
    for (type in c("test", "train")){
      if (type=="test") thisdata <- copy(test_data) else thisdata <- copy(surv)
      mar <- length(unique(thisdata[event_type=="mar"]$patient_id))
      aids <- length(unique(thisdata[event_type=="aids"]$patient_id))
      death <- length(unique(thisdata[event_type=="death"]$patient_id))
      summary <- rbind(summary, list(iteration, split, type, mar, aids, death))
    }
    
    save(surv, file=paste0(main_dir, iteration, "/", split, "/prepped_data.rdata"))
    save(test_data, file=paste0(main_dir, iteration, "/", split, "/test_data.rdata"))
    startval <- startval+split_size
  }
}

rm(thisdata); gc

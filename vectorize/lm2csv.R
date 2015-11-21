##write summary of lm to csv
lm2csv <- function(lm, fname='lm',confint=TRUE){

# grab the coefficients and append confint
lm.tb<-as.data.frame(summary(lm)$coefficients)
if(confint==TRUE){
lmt.tb<-cbind(lm.tb,as.data.frame(confint(lm)))
}

 
# format the table
summary = format(lm.tb, autoformat = 1)
 
# write it as a csv file 
write.csv(summary, paste(fname,"coefficients.csv", sep=''))
}

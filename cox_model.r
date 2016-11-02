## Test a cox model on the CASCADE survival data

library(data.table)
library(ggplot2)
library(survival)

main_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/cox_model/"

load(paste0(main_dir, "clean_data.rdata"))
source("ggsurv.R")

#keep pre-96 only
pre_96 <- data[pre_1996==1]
pre_96[, event:= ifelse(event_type=="death", 1, 0)] # pretend uninformative censoring
pre_96[, binned_age:= cut(serocon_age, breaks=c(15, 25, 35, 45, Inf), labels=c("15", "25", "35", "45"))]

## Plot age and event time/type
ggplot(pre_96, aes(x=serocon_age, y=log(event_time))) +
    geom_point(aes(color=as.factor(event)))

# Bivariate model
model <- coxph(Surv(event_time, event)~serocon_age, data=pre_96)

## coef: 0.02177, exp(coef) = 1.022. Interpretation: an extra year of age increases the hazard of an event by a factor of 1.022, or 2.2%
## --> An extra ten years of age increases hazard of an event by 22%

## Plot using binned ages 
newdata <- copy(pre_96)
newdata[, serocon_age:= as.numeric(binned_age)]


# Including sex and infection type
multivar_model <- coxph(Surv(event_time, event)~serocon_age+as.factor(sex)+as.factor(inf_mode), data=pre_96)

# Binning age: 15, 25, 35, 45
binned_model <- coxph(Surv(event_time, event)~binned_age, data=pre_96, model=T)
newdata <- unique(pre_96[, list(binned_age)])
fit <- survfit(binned_model, newdata=newdata)
plot(fit, lty=1:4)
ggsurv(fit)



####################################################################################################
## Author: Amelia Bertozzi-Villa
## Description: Plot time-to-event for ART initiators vs AIDS and death events, see how long people
##              wait before initiation
####################################################################################################


rm(list=ls())

library(data.table)
library(ggplot2)

#automatically finds correct directory name so we can stop commenting and uncommenting lines
dir.exists <- function(d) {
  de <- file.info(d)$isdir
  ifelse(is.na(de), FALSE, de)
}
root_dir <- ifelse(dir.exists("/home/cselinger/"), "/home/cselinger/HIV-Cascade/merge/data/", "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/")
main_dir <- ifelse(length(commandArgs())>2, commandArgs()[3], root_dir)

#load data
load(paste0(main_dir, "prepped_data.rdata"))
load(paste0(main_dir, "alldata.rdata"))

#determine more detailed event type
detailed_event <- unique(alldata[, list(patient_id, event_type, enroll_date, mar_date, serocon_date, aids_date, death_date, art_start_date, ltfu_date, admin_censor_date)])
detailed_event[event_type=="mar", event_type:= ifelse((!is.na(art_start_date) & art_start_date==mar_date), "art_start",
                                                     ifelse((!is.na(ltfu_date) & ltfu_date==mar_date), "ltfu", "admin_censor"))]

for (event_label in c("death", "aids", "art_start", "ltfu", "admin_censor")){
  
  detailed_event[event_type==event_label, event_time:=detailed_event[event_type==event_label][[paste0(event_label, "_date")]] - serocon_date]
}
detailed_event[,event_time:= as.numeric(event_time/365)]


density_plot <- ggplot(detailed_event[event_type %in% c("aids", "death", "art_start")], aes(x=event_time, group=event_type)) +
                geom_density(aes(color=event_type, fill=event_type), alpha=0.3, size=1)

print(density_plot)


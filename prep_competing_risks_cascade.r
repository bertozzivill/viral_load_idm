#################################################################################################
## prep_competing_risks_cascade.r
## using CASCADE data, prep competing risk model data (at the moment use only age as a covariate). 
## you will need:
## patient_id: 
## time: 
## event: possible events: AIDS death, other death, mono/combo therapy initiation, ltfu, censorship 
## serocon_age: age at diagnosis/seroconversion
## spvl: set point viral load

## tblART= "list(patient_id=PATIENT, art_drug=ART_ID, art_start=ART_SD, art_last=ART_ED)",
#################################################################################################

#setup
rm(list=ls())
set.seed(strtoi("0xCAFE"))

library(foreign)
library(data.table)
library(reshape2)
library(ggplot2)
require(bit64)

main_dir <- ("C:/Users/abertozzivilla/Dropbox (IDM)/viral_load/cascade/data/")

#################################################################
# I. import data
#################################################################

files <- list(tblBAS="list(patient_id=PATIENT, birth_date=BIRTH_D, first_visit_date=FRSVIS_D, enroll_date=ENROL_D, sex=GENDER, inf_mode=MODE, first_pos_date=HIV_POS_D, serocon_date=SEROCO_D, serohow=SEROHOW, aids_indic=AIDS_Y, aids_date=AIDS_D, art_indic=RECART_Y, art_start_date=RECART_D, last_artcheck_date=LTART_D, admin_censor_date=CENS_D)", 
              tblLAB_RNA= "list(patient_id=PATIENT, visit_date=RNA_D, vl=RNA_V, assay_ll=RNA_L, assay_ul=RNA_UL, assay_type=RNA_T)", 
              tblLAB_CD4="list(patient_id=PATIENT, visit_date=CD4_D, cd4_count=CD4_V, cd4_meas_type=CD4_U)",
              tblLTFU= "list(patient_id=PATIENT, ltfu_indic=DROP_Y, ltfu_date=DROP_D, death_indic=DEATH_Y, death_date=DEATH_D, cod_1=DEATH_R1, cod_2=DEATH_R2, cod_1_causal=DEATH_RC1, cod_2_causal=DEATH_RC2, last_alive_date=L_ALIVE)")

#generate a large list containing all data
data <- lapply(names(files), function(fname){
  print(fname)
  cols<-files[[fname]]
  #import<-fread(paste0("CASCADE/", fname, ".csv"))
  import<-fread(paste0(main_dir, "raw_data/", fname, ".csv"))
  import[, eval(parse(text=cols))] #keep only the columns we want, and change the names
})

baseline <- data[[1]]
viral <- data[[2]]
cd4 <- data[[3]]
ltfu <- data[[4]]

#merge into a single dataset

#baseline and ltfu should have the same number of entries (one/patient), merge those first
baseline <- baseline[, list(patient_id, birth_date, sex, inf_mode, serocon_date, aids_indic, aids_date, art_indic, art_start_date, admin_censor_date,enroll_date)]
ltfu <- ltfu[, list(patient_id, ltfu_date, death_indic, death_date)] #not using cause of death FOR NOW
data <- merge(baseline, ltfu, by=intersect(names(baseline), names(ltfu)), all=T)

############################################################
#viral load:
# - get rid of duplicate measurements
# -drop if vl is below detectable limit and no assay lower limit was given
# -sometimes the assay lower limit is set to a negative number.  assume this is supposed to be positive and switch it (?????)
# - if vl is below detectable limit and that limit is given, set vl to half of lower limit
# -keep VERY low vl? 
#############################################################

#drop null vl values
viral <- viral[!is.na(vl) & vl!=0]

#drop duplicated rows
viral <- viral[!duplicated(viral)]

#drop if vl==-1 (below detectable) and there's no assay lower limit
viral <- viral[vl>0 | (vl==-1 & !is.na(assay_ll))]

#change assay ll from negative to positive (check with someone about this)
viral[, assay_ll:= ifelse(assay_ll<0, assay_ll*-1, assay_ll)]

#set undetectable vl's to half of the assay ll
viral[, vl := ifelse(vl==-1, assay_ll/2, vl)]
viral[, visit_type:="vl"]

#convert visit_date to date object
viral[, visit_date:= as.Date(as.POSIXct(strptime(visit_date, format="%Y-%m-%d"), tz="GMT"))]

#########################################################################################
# II. Individual data:
#     --age at seroconversion (years)
#     -- event type (AIDS/death or censoring due to AIDS initiation/study termination)
#     -- time from seroconversion to event (years) 
# this requires removing some nonsensical/unknown results:
#  -- aids status unknown 
#  -- patient on art, but no art initiation date
#  --age under 15
#  -- status other than 1 (homosexual) or 2=6 (heterosexua) transmission
##########################################################################################


#convert strings into 'Date' objects
# sometimes, for no reason I know of, dates get coded to 1911. The earliest date in this 
# series should always be the patient's birth date; drop anything where this is not the case.
print(paste("converting birth_date to date object"))
data[['birth_date']] <- as.Date(as.POSIXct(strptime(data[['birth_date']], format="%Y-%m-%d"), tz="GMT"))

for (colname in names(data)){
  if (grepl("_date", colname) & colname!="birth_date"){
    print(paste("converting", colname, "to date object"))
    data[[colname]] <- as.Date(as.POSIXct(strptime(data[[colname]], format="%Y-%m-%d"), tz="GMT"))
    data<- data[is.na(data[[colname]]) | data[[colname]]>=data[['birth_date']]]
  }
}

#calculate age at seroconversion
data[, serocon_age:= as.numeric((serocon_date - birth_date)/365)]

#drop those younger than 15
data<- data[serocon_age>=15]

#drop those with an infection mode other than 1 or 6
data<- data[inf_mode==1 | inf_mode==6]

### Set up indicators for each event type; those for art, death, and AIDS are already in place. We need to remove impossible combinations
### and create indicators for loss to follow-up and administrative censoring
data[, ltfu_indic:= ifelse(is.na(ltfu_date), 0, 1)]
# administrative censoring: if a patient did not die, get AIDS, or go on ART, 
# we can use their administrative censoring date as their 'event' date
data[, admin_censor_indic:= ifelse(aids_indic==1| death_indic==1|art_indic==1|ltfu_indic==1, 0, 1)]

#drop those with an indicators of 1 but no event date
data <- data[art_indic==0 | (art_indic==1 & !is.na(art_start_date))]
data <- data[death_indic==0 | (death_indic==1 & !is.na(death_date))]
data <- data[aids_indic==0 | (aids_indic==1 & !is.na(aids_date))] #also drops those with unknown aids status
data <- data[ltfu_indic==0 | (ltfu_indic==1 & ltfu_date<=admin_censor_date)] #ltfu should never occur after administrative censoring

#now, to sort out the timing of the dates. we start by determining the earlies mar date (ltfu, or admin_censor)
data[, mar_date:=as.Date(pmin(ltfu_date, admin_censor_date, na.rm=T), origin="1970-01-01")]

# think about determining event types as a decision tree: the first split is on whether or not the patient 
# initiated ART

## PATIENT INITIATED ART
# If patient initiated ART and an AIDS event preceded ART, event is AIDS
data[art_indic==1 & aids_indic==1 & aids_date<art_start_date, event_type:="aids"]
# If ART precedes AIDS and/or death, or there is no logged AIDS/death event, event is ART
data[art_indic==1 & (art_start_date<=pmin(aids_date, death_date, na.rm=T) | (aids_indic==0 & death_indic==0)), event_type:="art"]

## PATIENT DID NOT INITIATE ART
# In the absence of ART, the second split is on whether the patient had AIDS

##    PATIENT HAD AIDS
# If patient had AIDS and died, event type is death
data[art_indic==0 & aids_indic==1 & death_indic==1, event_type:="death"]
# If patient had AIDS and did not die, event type is AIDS
data[art_indic==0 & aids_indic==1 & death_indic==0, event_type:="aids"]

##    PATIENT DID NOT HAVE AIDS
# If patient died, event type is death
data[art_indic==0 & aids_indic==0 & death_indic==1, event_type:="death"]
# If patient did not die, event type is MAR
data[art_indic==0 & aids_indic==0 & death_indic==0, event_type:="mar"]

data[, event_date:= ifelse(event_type=="mar", mar_date,
                                  ifelse(event_type=="art", art_start_date,
                                         ifelse(event_type=="aids", aids_date, death_date)))]
data[, event_date:=as.Date(event_date, origin="1970-01-01")]


# calculate "event time" and debiased versions (to correct for survival bias)
data[, event_time:= as.numeric((event_date - serocon_date)/365)] #original "event time": time from recorded seroconversion to event
data[, bias:= as.numeric((enroll_date - serocon_date)/365)] # "bias" factor: time from recorded seroconversion to enrollment (if it's positive, it's "bias")
data[, bias:= max(bias, 0), by=rownames(data)]
data[, debiased_event_time:= event_time-bias] # "corrected" event time: time from enrollment to event

## DROP INDIVIDUALS HERE: exclude patients who had an event before, or when they enrolled (e.g. enrolling in ART the same day as enrolment)
data <- data[event_date>enroll_date]
# Also exclude individuals who record an event prior to seroconversion (event_date<0)
data <- data[event_time>0]

## Create a pre-1996 indicator
data[, pre_1996:= ifelse(event_date<as.Date("1996-01-01"), 1, 0)]

# Don't worry about viral load for now, just save this
save(data, file=paste0(main_dir, "cox_model/surv_data.rdata"))

# Also save a cleaned version
data <- data[, list(patient_id, sex, inf_mode, serocon_age, event_type, event_time, debiased_event_time, pre_1996)]
save(data, file=paste0(main_dir, "cox_model/clean_data.rdata"))


# #########################################################################################
# # III. Merge vl and individual together
# ##########################################################################################
# 
# alldata <- merge(data, viral, by="patient_id", all.x=T) #keep even those with no vl measurement
# alldata[, visit_time:= as.numeric((visit_date - serocon_date)/365)]
# 
# #drop visits prior to seroconversion
# alldata <- alldata[visit_time>=0 | is.na(visit_time)]
# 
# # drop visits after event
# alldata <- alldata[visit_date<event_date | is.na(visit_date)]
# 
# #save full dataset
# save(alldata, file=paste0(main_dir,"alldata.rdata"))
# 
# #################################################################
# # V. Split into relevant datasets, save
# #################################################################
# 
# ## viral load
# vl <- alldata[,list(patient_id, time=visit_time, vl, assay_ll, assay_type)]
# write.csv(vl, file=paste0(main_dir, "vl.csv"), row.names=F)
# 
# ## survival
# surv <- alldata[, list(patient_id, event_type, event_time, event_timeNew, agesero=serocon_age,event_date,serocon_date,enroll_date,seroconv_before_enroll,aids_before_enroll,bias)]
# setkeyv(surv, NULL)
# surv <- unique(surv)
# write.csv(surv, file=paste0(main_dir, "surv.csv"), row.names=F)
# 
# #################################################################
# # VI. Source code to prep covariates
# #################################################################
# 
# source("prep_covariates.r")
##plotting figures
library(data.table)
#main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/"
#main_dir <- "C:/Users/cselinger/Dropbox (IDM)/viral_loadNew/cascade/data/"
#main_dir <- "C:/Users/abertozz/Dropbox (IDM)/viral_load/cascade/data/cross_validation/1/1/"
main_dir <-"/home/cselinger/HIV-Cascade/merge/data/"

#load plot scripts and data
sapply(c("multiplot.R","ggsurv.R","plot.KM.curves.R","plot.spvl.agesero.distribution.R","plot.negative.slope_vl.measurements.R","plot.surface.R","plot.survival.R","plot.modelcurve.R"),source)
surv <- fread(paste0(main_dir, "surv.csv")) #column 'AD' corresponds to 'aids_death_indic', 'ADtime' to 'event_time' in alldata.rdata
load(paste0(main_dir,"bestmodel.Rdata"))
load(paste0(main_dir,"missing_spvl_model.Rdata"))

#Figure 1
#png('../figures/Figure1a.png',width=20, height=20, units='cm', res=400)
#plot.KM.curves(surv)
#dev.off()


#Figure 2a
png('../figures/Figure2a.png',width=20, height=15, units='cm', res=400)
plot.spvl.agesero.distribution(bestmodel$data)
dev.off()

#Figure 2b
png('../figures/Figure2b.png',width=35, height=20, units='cm', res=400)
plot.negative.slope_vl.measurements(vl,missing,re_vl)
dev.off()

#Figure 3b
png('../figures/Figure3b.png',width=20, height=15, units='cm', res=400)
plot.surface(bestmodel$lm,bestmodel$name)
dev.off()

#Figure 3bis
png('../figures/Figure3bis.png',width=45, height=35, units='cm', res=400)
plot.survival(bestmodel$lm,bestmodel$name)
dev.off()

#Figure 3d
png('../figures/Figure3d.png',width=35, height=30, units='cm', res=400)
plot.modelcurve(bestmodel$lm,bestmodel$name)
dev.off()


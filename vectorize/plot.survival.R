plot.survival<-function(model,bestmodel.name,title='title'){

##PLOTS####
  ##survival surface

  spvl_method<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[1]
  interaction<-as.numeric(unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[2])
  bins<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[3]
  
  z1 <- seq(20,60,0.1);z2<-seq(2,6,0.1);z3=as.factor(c(1,2));#z4=as.factor(levels(data$agebin))
  newdf <- expand.grid(agesero=z1,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
  
  data.surface<-predict(model, newdf,se.fit=TRUE)
  surface=data.table(transform(newdf, survival=data.surface$fit));surface$se=data.surface$se.fit
  surface[,survival:=exp(survival)];surface[,se:=100*(exp(se)-1)];surface[,se:=round(se,4)]; setnames(surface,spvl_method,'spvl')
  levels(surface$event_num)<-c("AIDS","Death")

  library(directlabels)
  library(ggplot2)

  
  plot.surface <- ggplot(surface, aes(y=agesero, x=spvl, z = survival))+
    geom_tile(aes(fill=survival))+
    scale_fill_gradient(low="blue", high="orange")+
    stat_contour(aes(colour=..level..),binwidth = 1,colour="black",linetype=3)+
    labs(y="age at seroconversion (years)", x=paste0("setpoint viral load (log10): ",spvl_method), title=paste0('survival surface \n ',c('without interaction', 'with interaction ')[interaction+1])) + 
    theme_bw()+
    facet_grid(. ~event_num)
  plot.surface <-direct.label(plot.surface,"top.pieces")




z1 <- seq(20,60,10);z2<-seq(2,6,1);z3=as.factor(c(1,2))
newdf <- expand.grid(agesero=z1,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
data.agesero<-data.table(transform(newdf, survival=predict(model, newdf,se.fit=TRUE)));setnames(data.agesero,spvl_method,'spvl')
levels(data.agesero$event_num)<-c("AIDS","Death")
dodge <- position_dodge(width=0.5)
limits <- aes(x=agesero,ymax = exp(survival.fit + survival.se.fit), ymin=exp(survival.fit - survival.se.fit))



  plot.agesero<-ggplot(data=data.agesero)+
    geom_point(aes(y=exp(survival.fit), x=agesero, color=factor(spvl)))+
    stat_smooth(aes(y=exp(survival.fit), x=agesero,color=factor(spvl),group=as.factor(spvl)),se=FALSE,method='lm')+
    scale_colour_discrete(name=paste(spvl_method,'(log10)')) + 
    geom_errorbar(limits, position=dodge, width=0.25)+
    labs(x="age at seroconversion (years)", y="average predicted survival (years)",title=paste0('survival \n ',c('without interaction ', 'with interaction ')[interaction+1])) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12))+
    facet_grid(. ~event_num)+labs(spvl='custom title')
  
  source('multiplot.R')

  return(multiplot(plot.surface,plot.agesero,cols=2,maintitle = title))
}



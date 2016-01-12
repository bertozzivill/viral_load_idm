
plot.surface<-function(bestmodel.lm,bestmodel.name,title='survival surface'){
  
  ##PLOTS####
  coeff<-bestmodel.lm$coefficients
  moderation.agesero<-function(spvl=5,agesero){coeff[1]+(coeff[2]+coeff[4]*spvl)*agesero}
  moderation.spvl<-function(agesero=30,spvl){coeff[1]+(coeff[3]+coeff[4]*agesero)*spvl}
  data<-na.omit(bestmodel$data)
  
  ##PLOTS####
  ##survival surface
  
  spvl_method<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[1]
  interaction<-as.numeric(unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[2])
  bins<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[3]
  
  z1 <- seq(20,60,0.1);z2<-seq(2,6,0.1);z3=as.factor(c(1,2));#z4=as.factor(levels(data$agebin))
  newdf <- expand.grid(agesero=z1,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
  
  data.surface<-predict(bestmodel.lm, newdf,se.fit=TRUE)
  surface=data.table(transform(newdf, survival=data.surface$fit));surface$se=data.surface$se.fit
  surface[,survival:=exp(survival)];surface[,se:=100*(exp(se)-1)];surface[,se:=round(se,4)]; setnames(surface,spvl_method,'spvl')
  levels(surface$event_num)<-c("AIDS","Death")[as.numeric(levels(surface$event_num))]
  
  library(directlabels)
  library(ggplot2)
  

plot.surface <- ggplot(surface, aes(y=agesero, x=spvl, z = survival))+
  geom_tile(aes(fill=survival))+
  scale_fill_gradient(low="blue", high="orange")+
  stat_contour(aes(colour=..level..),binwidth = 1,colour="black",linetype=1,linesize=5)+
  labs(y="age at seroconversion (years)", x=paste0("setpoint viral load (log10): ",spvl_method), title=paste0('survival surface \n ',c('without interaction', 'with interaction ')[interaction+1])) + 
  theme_bw()+
  facet_grid(. ~event_num)
plot.surface <-direct.label(plot.surface,"top.pieces")

return(plot.surface)
}



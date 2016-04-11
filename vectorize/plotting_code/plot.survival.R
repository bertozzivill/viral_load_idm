

plot.survival<-function(bestmodel.lm=bestmodel$lm,bestmodel.name=bestmodel$name,title='moderation and predicted survival'){

  coeff<-bestmodel.lm$coefficients
  moderation.agesero<-function(spvl=5,agesero){coeff[1]+(coeff[2]+coeff[4]*spvl)*agesero}
  moderation.spvl<-function(agesero=30,spvl){coeff[1]+(coeff[3]+coeff[4]*agesero)*spvl}
  data<-na.omit(bestmodel$data)
 
##PLOTS####
  ##survival surface

  spvl_method<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[1]
  interaction<-as.numeric(unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[2])
  bins<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[3]
  
  z1 <- seq(20,60,0.1);z2<-seq(2,6,0.1);z3=as.factor(levels(data$event_num));#z4=as.factor(levels(data$agebin))
  newdf <- expand.grid(agesero=z1,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
  
  data.surface<-predict(bestmodel.lm, newdf,se.fit=TRUE)
  surface=data.table(transform(newdf, survival=data.surface$fit));surface$se=data.surface$se.fit
  surface[,survival:=exp(survival)];surface[,se:=100*(exp(se)-1)];surface[,se:=round(se,4)]; setnames(surface,spvl_method,'spvl')
  levels(surface$event_num)<-c("AIDS","Death")[as.numeric(levels(data$event_num))]


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




z1 <- seq(20,60,10);z2<-seq(2,6,1);z3=as.factor(levels(data$event_num))
newdf <- expand.grid(agesero=z1,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
data.agesero<-data.table(transform(newdf, survival=predict(bestmodel.lm, newdf,se.fit=TRUE)));setnames(data.agesero,spvl_method,'spvl')
levels(data.agesero$event_num)<-c("AIDS","Death")[as.numeric(levels(data$event_num))]


dodge <- position_dodge(width=0.5)
limits <- aes(x=agesero,ymax = exp(survival.fit + survival.se.fit), ymin=exp(survival.fit - survival.se.fit))



  plot.agesero<-ggplot(data=data.agesero)+
    geom_point(aes(y=exp(survival.fit), x=agesero, color=factor(spvl)))+
    geom_line(aes(y=exp(survival.fit), x=agesero,color=factor(spvl),group=as.factor(spvl)))+
    #stat_smooth(aes(y=exp(survival.fit), x=agesero,color=factor(spvl),group=as.factor(spvl)),se=FALSE,method='lm')+
    scale_colour_discrete(name=paste(spvl_method,'(log10)')) + 
    geom_errorbar(limits, position=dodge, width=2.5)+
    labs(x="age at seroconversion (years)", y="average predicted survival (years)",title=paste0('survival \n ',c('without interaction ', 'with interaction ')[interaction+1])) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0,size=12))+
    labs(spvl='custom title')+
    facet_grid(. ~event_num)+
    scale_y_continuous(breaks=seq(5,20,1))
    
 
  
  
  limits <- aes(x=spvl,ymax = exp(survival.fit + survival.se.fit), ymin=exp(survival.fit - survival.se.fit))
  dodge <- position_dodge(width=5)
  
  plot.spvl<-ggplot(data=data.agesero)+
    geom_point(aes(y=exp(survival.fit), x=spvl, color=factor(agesero)))+
    geom_line(aes(y=exp(survival.fit), x=spvl,color=factor(agesero),group=as.factor(agesero)))+
    #stat_smooth(aes(y=exp(survival.fit), x=spvl,color=factor(agesero),group=as.factor(agesero)),se=FALSE,method='lm')+
    scale_colour_discrete(name='agesero (years)') + 
    geom_errorbar(limits, position=dodge, width=0.25)+
    labs(x="set point viral load (log10)", y="average predicted survival (years)",title=paste0('survival \n ',c('without interaction ', 'with interaction ')[interaction+1])) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0,size=12))+
    labs(spvl='custom title')+
    facet_grid(. ~event_num)+
    scale_y_continuous(breaks=seq(5,20,1))
  

  mod.agesero<-data.frame(
  agesero=seq(20,60,10),
  min=moderation.agesero(spvl=min(data$spvl_model),agesero=seq(20,60,10)),
  max=moderation.agesero(spvl=max(data$spvl_model),agesero=seq(20,60,10))
  )
  #levels(mod.agesero$event_num)<-c("AIDS","Death")[as.numeric(levels(data$event_num))]
  

  mod.spvl<-data.frame(
  spvl=c(2:7),
  min=moderation.spvl(agesero=min(data$agesero),spvl=c(2:7)),
  max=moderation.spvl(agesero=max(data$agesero),spvl=c(2:7))
  )
  #levels(mod.spvl$event_num)<-c("AIDS","Death")[as.numeric(levels(data$event_num))]


  mod.agesero<-melt(mod.agesero,id.vars='agesero')
  mod.spvl<-melt(mod.spvl,id.vars='spvl')
 
  plot.mod.agesero<-
    ggplot(data=mod.agesero,aes(x=agesero,y=value,group=as.factor(variable),color=as.factor(variable)))+
    geom_point()+geom_line()+
    scale_colour_discrete(name='spvl (log10)') + 
    ylab("change in survival (log years)")+
    ggtitle("conditional fixed effect of \nage at seroconversion (by spvl) on survival")+
    theme_bw()
  
  plot.mod.spvl<-
  ggplot(data=mod.spvl,aes(x=spvl,y=value,group=as.factor(variable),color=as.factor(variable)))+
    geom_point()+geom_line()+
    scale_colour_discrete(name='agesero (years)') + 
    ylab("change in survival (log years)")+
    ggtitle("conditional fixed effect of \nset point viral load (by agesero) on survival")+
    theme_bw()
  
  
  
  index=expand.grid(spvl=c(2:7),agesero=seq(20,60,10))
  index$moderation.spvl<-mapply(moderation.spvl,agesero=index$agesero,spvl=index$spvl)
  

  plot.mod.spvl<-
    ggplot(data=index,aes(x=spvl,y=moderation.spvl,color=as.factor(agesero),group=as.factor(agesero)))+
    geom_point()+geom_line()+
    scale_colour_discrete(name='agesero (years)') + 
    ylab("change in survival (log years)")+
    ggtitle("conditional fixed effect of set point viral load (by agesero) on survival")+
    theme_bw()
  
  #sjp.int(bestmodel.lm,axisLimits.x =c(2,8),showCI=T)

  
  library(interplot)
  plot.mod.agesero<-interplot(bestmodel.lm,var1='spvl_model',var2='agesero')+
  xlab("age at seroconversion (years)") +
  ylab("estimated coefficent for set point viral load") +
  theme_bw() +
  ggtitle("Estimated Coefficient of set point viral load \non survival by age at seroconversion") +
  theme(plot.title = element_text(face="bold")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = median(data$agesero), color='red')+
  annotate("text",x=median(data$agesero)+8,y=-0.1,label=paste0("median=",round(median(data$agesero),2)),color="red",size=5)
  
  
  plot.mod.spvl<-interplot(bestmodel.lm,var1='agesero',var2='spvl_model')+
    ylab("estimated coefficent for age at seroconversion") +
    xlab("set point viral load (log10)") +
    theme_bw() +
    ggtitle("Estimated Coefficient of age at seroconversion \non survival by set point viral load") +
    theme(plot.title = element_text(face="bold")) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    xlim(c(2,7.5))+
    geom_vline(xintercept = median(data$spvl_model), color='red')+
    annotate("text",x=median(data$spvl_model)+0.8,y=-0.008,label=paste0("median=",round(median(data$spvl_model),2)),color="red",size=5)
  
  
  
  source('multiplot.R')

  return(multiplot(plot.mod.spvl,plot.mod.agesero,plot.spvl,plot.agesero,layout=matrix(c(1,3,2,4),nrow=2,byrow=TRUE),maintitle = title))
}



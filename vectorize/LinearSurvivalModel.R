LinearSurvivalModel<-function(#data,
                              spvl_method='spvl_model',
                              interaction=0,
                              bins=c(seq(15,65,10),100),
                              return.modelobject=0){
  
  bins<-unlist(bins)
  data<-data.table(data)
#   print(paste0("type: ",typeof(data)))
#   print(paste0("class: ",class(data)))
  ###agebins for categorical variable
  if(length(bins)>1){
    data[,agebin:=cut(data$agesero, bins, include.lowest=TRUE)]
  }
  
  
  ###arguments of function
  
  binning<-ifelse(length(bins)>1,
         paste0(" and age bins ", paste(bins,collapse=' ')),
         ' and continous age covariate')
         
  print(paste0('SURVIVAL: ',"log normal survival with ",spvl_method,c(' without interaction',' with interaction')[interaction+1],binning))
  
  
  if (interaction==1&length(bins)>1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  (agebin + ",spvl_method,")^2 + event_num")}
  if (interaction==1&length(bins)>1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  (agebin + ",spvl_method,")^2")}
  
  if (interaction==1&length(bins)==1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  (agesero + ",spvl_method,")^2 + event_num")}
  if (interaction==1&length(bins)==1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  (agesero + ",spvl_method,")^2")}
  
  if (interaction==0&length(bins)>1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  agebin + ",spvl_method," + event_num")}
  if (interaction==0&length(bins)>1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  agebin + ",spvl_method)}
  
  if (interaction==0&length(bins)==1&length(levels(data$event_num))>1){formula <- paste0("observed_survival~  agesero + ",spvl_method," + event_num")}
  if (interaction==0&length(bins)==1&length(levels(data$event_num))==1){formula <- paste0("observed_survival~  agesero + ",spvl_method)}
  
  

  
  
  ###evaluate model
  output <- lm(as.formula(formula), data=data)
  if(return.modelobject==0){
    return(list('summary'=summary(output),'AIC'=round(AIC(output))))
  }else{
    return(output)
  }

}

#   
#   #smalldata <- data[which(!is.na(data[[spvl_method]])),] #use only valid spvl estimates for survival model
# 
#   
#   ###PLOTS####
#   ##survival surface
#   z1 <- seq(20,60,0.1);z2<-z2<-seq(2,6,0.1);z3=as.factor(c(1,2));#z4=as.factor(levels(data$agebin))
#   newdf <- expand.grid(agesero=z1,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
#   
#   if (length(bins)>1){
#   data.surface<-predict(output, newdf,se.fit=TRUE)
#   surface=data.table(transform(newdf, survival=data.surface$fit));surface$se=data.surface$se.fit
#   surface[,survival:=exp(survival)];surface[,se:=100*(exp(se)-1)];surface[,se:=round(se,4)]; setnames(surface,spvl_method,'spvl')
#   
#   
#   library(directlabels)
#   library(ggplot2)
#   
#   plot.surface <- ggplot(surface[event_num==1], aes(y=agesero, x=spvl, z = survival))+
#     geom_tile(aes(fill=survival))+
#     scale_fill_gradient(low="blue", high="orange",limits=c(1,20))+
#     stat_contour(aes(colour=..level..),binwidth = 1,colour="black",linetype=3)+
#     labs(y="age at seroconversion (years)", x=paste0("setpoint viral load (log10): ",spvl_method), title=paste0('survival surface \n ',c('without interaction', 'with interaction ')[interaction+1],'\n event: AIDS')) + theme_bw()
#   plot.surface.aids <-direct.label(plot.surface,"top.pieces")
#   
#   plot.surface <- ggplot(surface[event_num==2], aes(y=agesero, x=spvl, z = survival))+
#     geom_tile(aes(fill=survival))+
#     scale_fill_gradient(low="blue", high="orange",limits=c(1,20))+
#     stat_contour(aes(colour=..level..),binwidth = 1,colour="black",linetype=3)+
#     labs(y="age at seroconversion (years)", x=paste0("setpoint viral load (log10): ",spvl_method), title=paste0('survival surface \n ',c('without interaction', 'with interaction ')[interaction+1],'\n event: Death')) + theme_bw()
#   plot.surface.death <-direct.label(plot.surface,"top.pieces")
#   }
# 
#   
#   
#   
#   
#   ##survival curves for agebins and event_num
#   z1 <- seq(20,60,0.1);z2<-seq(2,6,1);z3=as.factor(c(1,2));z4=as.factor(levels(data$agebin))
#   newdf <- expand.grid(agebin=z4,spvl_model=z2,event_num=z3);colnames(newdf)[2]<-spvl_method
#   
#   
#   data.agebin<-data.table(transform(newdf, survival=predict(output_2, newdf,se.fit=TRUE)));setnames(data.agebin,spvl_method,'spvl')
#   data.agebin$agebin = factor(data.agebin$agebin,levels(data.agebin$agebin)[c(6,1:5)])
#   dodge <- position_dodge(width=0.5)
#   
#   limits.aids <- aes(x=agebin,ymax = exp(data.agebin[event_num==1]$survival.fit + data.agebin[event_num==1]$survival.se.fit), ymin=exp(data.agebin[event_num==1]$survival.fit - data.agebin[event_num==1]$survival.se.fit))
#   
#   plot.agebin.aids<-ggplot(data=data.agebin[event_num==1])+
#     geom_point(aes(y=exp(survival.fit), x=agebin, color=factor(spvl),shape=factor(event_num)))+
#     stat_smooth(aes(y=exp(survival.fit), x=agebin,color=factor(spvl),group=as.factor(spvl)),se=FALSE)+
#     scale_colour_discrete(name=spvl_method) + 
#     scale_x_discrete(limits=levels(data.agebin$agebin))+
#     geom_errorbar(limits.aids, position=dodge, width=0.25)+
#     labs(x="agesero", y="average predicted survival",title=paste0('survival \n ',c('without interaction ', 'with interaction ')[interaction+1],'\n event: AIDS \n ','AIC=',round(AIC(output_2)))) + theme_bw()+ylim(c(1,20))+
#     theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12))
#   
#   
#   limits.death <- aes(x=agebin,ymax = exp(data.agebin[event_num==2]$survival.fit + data.agebin[event_num==2]$survival.se.fit), ymin=exp(data.agebin[event_num==2]$survival.fit - data.agebin[event_num==2]$survival.se.fit))
#   
#   plot.agebin.death<-ggplot(data=data.agebin[event_num==2])+
#     geom_point(aes(y=exp(survival.fit), x=agebin, color=factor(spvl),shape=factor(event_num)))+
#     stat_smooth(aes(y=exp(survival.fit), x=agebin,color=factor(spvl),group=as.factor(spvl)),se=FALSE)+
#     scale_colour_discrete(name=spvl_method) + 
#     scale_x_discrete(limits=levels(data.agebin$agebin))+
#     geom_errorbar(limits.death, position=dodge, width=0.25)+
#     labs(x="agesero", y="average predicted survival",title=paste0('survival \n ',c('without interaction ', 'with interaction ')[interaction+1],'\n event: Death \n ','AIC=',round(AIC(output_2)))) + theme_bw()+ylim(c(1,20))+
#     theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12))
#   
#   
#   #   ###aggreagte survival curves for even_num confounded
#   #   limits <- aes(x=agebin,ymax = exp(data.agebin$survival.fit + data.agebin$survival.se.fit), ymin=exp(data.agebin$survival.fit - data.agebin$survival.se.fit))
#   #   
#   #   plot.agebin<-ggplot(data=data.agebin)+
#   #     geom_point(aes(y=exp(survival.fit), x=agebin, color=factor(spvl),shape=factor(event_num)))+
#   #     stat_smooth(aes(y=exp(survival.fit), x=agebin,color=factor(spvl),group=factor(spvl)),se=FALSE)+
#   #     scale_colour_discrete(name=spvl_method) + 
#   #     scale_x_discrete(limits=levels(data.agebin$agebin))+
#   #     #geom_errorbar(limits, position=dodge, width=0.25)+
#   #     labs(x="agesero", y="average predicted survival",title=paste0('survival \n ',c('without interaction ', 'with interaction ')[interaction+1],'\n event: Death \n ','IC=',round(AIC(output_2)))) + theme_bw()+ylim(c(1,20))+
#   #     theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12))
#   #   
#   
#   #   
#   #   png(paste0('log_normal_survival_',c("biased","debiased")[debiased+1],'_',spvl_method,'.png'),units='cm',width=32,height=32,res=400)
#   #   multiplot(plot.agebin,plot.surface,plot.spvl,plot.surface.se,cols=2)
#   #   dev.off()
#   
#   ##prepare for output
#   PLOTS.summary<-list(filename=paste0('output\\PLOTS.summary_',runname,'_',spvl_method,'_interaction.',interaction,'.RData'),plot.surface.aids=plot.surface.aids,plot.surface.death=plot.surface.death,plot.agebin.aids=plot.agebin.aids,plot.agebin.death=plot.agebin.death,model.agesero=output_1,summary.agesero=summary(output_1),summary.agebin=summary(output_2),model.agebin=output_2,agesero.AIC=round(AIC(output_1)),agebin.AIC=round(AIC(output_2)),data.agebin=data.agebin)
#   return(PLOTS.summary)
#   save(PLOTS.summary,file=paste0('output\\PLOTS.summary_',runname,'_',spvl_method,'_interaction.',interaction,'_agesero.AIC.',round(AIC(output_1)),'_agebin.AIC.',round(AIC(output_2)),'.RData'))
# }

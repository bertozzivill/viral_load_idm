plot.modelcurve<-function(model,bestmodel.name,title='title'){

  ##PLOTS####
  ##model curve
  spvl_method<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[1]
  interaction<-as.numeric(unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[2])
  bins<-unlist(strsplit(as.character(bestmodel.name),"-",fixed=TRUE))[3]
  
  data=eval(getCall(model)$data,environment(formula(model)))
  data$predicted_survival<-predict(model)
  
  df.plot<-melt(data[,list(agesero,event_num,observed_survival,predicted_survival)],id=c('agesero','event_num'))
  levels(df.plot$event_num)<-c('AIDS','Death')
  
  ggplot(df.plot)+
    geom_point(aes(agesero,exp(value),group=variable,color=variable),alpha=0.3)+
    geom_smooth(data=df.plot[variable=='predicted_survival'],aes(agesero,exp(value),color='black'),method='lm')+
    labs(x="age at seroconversion (years)", y="survival until event (years)", title=paste0('survival  \n ',c('without interaction', 'with interaction ')[interaction+1],'\n',c('continuous agesero ', 'binned agesero ')[ifelse(length(bins)==1,1,2)])) + 
    theme_bw()+
    scale_x_continuous(breaks=seq(20,80,10))+
    scale_y_continuous(breaks=seq(0,40,10))+
    facet_grid(. ~event_num)
 
}
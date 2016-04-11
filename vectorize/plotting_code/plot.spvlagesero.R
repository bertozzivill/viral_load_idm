plot.spvlagesero<-function(model,bestmodel.name,title='title'){
  
  
  data=eval(getCall(model)$data,environment(formula(model)))

  df.plot<-melt(data[,list(agesero,event_num,spvl_model,spvl_fraser,spvl_hybrid)],id=c('agesero','event_num'))
  levels(df.plot$event_num)<-c('AIDS','Death')
  
  
  ggplot(df.plot)+
    geom_point(aes(agesero,value,group=variable,color=variable),alpha=0.3)+
    geom_smooth(aes(agesero,value,group=variable,color=variable),method='lm')+
    labs(x="age at seroconversion (years)", y="set point viral load (log10 RNA copies per mm^3)", title=title) +
    scale_x_continuous(breaks=seq(20,80,10))+
    scale_y_continuous(breaks=seq(0,8,1))+
    facet_grid(. ~event_num)
}
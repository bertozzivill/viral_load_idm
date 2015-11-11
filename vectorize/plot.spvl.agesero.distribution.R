library(ggplot2)
library(gridExtra)
plot.spvl.agesero.distribution<-function(data){

df=melt(data,id.vars=c("patient_id","agesero"),measure.vars=c("spvl_model","spvl_fraser"))  

hist_right <- ggplot(df)+geom_density(aes(value,group=as.factor(variable),fill=as.factor(variable)),alpha=0.4)+
  coord_flip()+theme_bw()+theme(legend.position = "none")

empty <- ggplot()+geom_point(aes(1,1), colour="white")+theme(
  axis.line = element_blank(), 
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())

scatter <- ggplot(data=df,aes(x=agesero,y=value,group=as.factor(variable),color=as.factor(variable)))+
  geom_point(alpha=0.1)+
  stat_smooth(size=4,method="lm")+
  xlab("age at seroconversion")+
  ylab("set point viral load")+
  theme_bw()+
  theme(legend.position = c(0.8,0.2))

hist_top <- ggplot(df)+geom_density(aes(agesero))+theme_bw()+
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
  

grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

}